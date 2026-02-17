"""
GO Tools FastAPI Application.

Provides REST API for GO Term Finder and GO Slim Mapper analysis.
Replaces the previous Flask + Perl implementation with pure Python.
"""
from __future__ import annotations

import hashlib
import os
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Optional

import boto3
from fastapi import FastAPI, Form, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse

from go_parser import GOOntology, GeneAssociations, load_slim_terms
from go_term_finder_service import (
    GOTermFinderService,
    CorrectionMethod,
    format_pvalue,
)
from go_slim_mapper_service import (
    GOSlimMapperService,
    parse_slim_terms_from_string,
)
from schemas import (
    EnrichedTermResponse,
    EnrichmentTermResponse,
    GeneHitResponse,
    GOSlimMapperResponse,
    GOTermFinderResponse,
    MappedGeneResponse,
    MappedSlimTermResponse,
)
from html_generator import (
    generate_term_finder_html,
    generate_slim_mapper_html,
    generate_term_finder_tab,
    generate_slim_mapper_tab,
)


# Configuration
DATA_DIR = os.environ.get("DATA_DIR", "/var/www/data/")
TMP_DIR = os.environ.get("TMP_DIR", "/var/www/tmp/")
S3_BUCKET = os.environ.get("S3_BUCKET", "")
CACHE_DIR = os.environ.get("CACHE_DIR", "/var/www/cache/")

# Global service instances
ontology: Optional[GOOntology] = None
associations: Optional[GeneAssociations] = None
term_finder_service: Optional[GOTermFinderService] = None
slim_mapper_service: Optional[GOSlimMapperService] = None


def load_data():
    """Load ontology and gene associations data."""
    global ontology, associations, term_finder_service, slim_mapper_service

    # Ensure cache directory exists
    os.makedirs(CACHE_DIR, exist_ok=True)
    os.makedirs(TMP_DIR, exist_ok=True)

    obo_path = os.path.join(DATA_DIR, "gene_ontology.obo")
    gaf_path = os.path.join(DATA_DIR, "gene_association.sgd")
    cache_path = os.path.join(CACHE_DIR, "ontology.pkl")

    # Load ontology
    ontology = GOOntology()
    ontology.load_obo(obo_path, cache_path)

    # Load gene associations
    associations = GeneAssociations(ontology)
    associations.load_gaf(gaf_path)

    # Initialize services
    term_finder_service = GOTermFinderService(ontology, associations)
    slim_mapper_service = GOSlimMapperService(ontology, associations)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan handler for loading data on startup."""
    load_data()
    yield


# Create FastAPI app
app = FastAPI(
    title="GO Tools API",
    description="GO Term Finder and GO Slim Mapper analysis tools",
    version="2.0.0",
    lifespan=lifespan,
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# =============================================================================
# Utility Functions
# =============================================================================

def get_s3_url(filename: str) -> str:
    """Get S3 URL for a file."""
    if S3_BUCKET:
        return f"https://{S3_BUCKET}.s3.amazonaws.com/gotermfinder/{filename}"
    return f"/download/{filename}"


def upload_to_s3(filepath: str, filename: str) -> str:
    """Upload a file to S3 and return the URL."""
    if not S3_BUCKET:
        return f"/download/{filename}"

    s3 = boto3.client("s3")
    s3_key = f"gotermfinder/{filename}"

    with open(filepath, "rb") as f:
        s3.upload_fileobj(f, S3_BUCKET, s3_key, ExtraArgs={"ACL": "public-read"})

    return f"https://{S3_BUCKET}.s3.amazonaws.com/{s3_key}"


def save_and_upload(content: str, base_name: str, extension: str) -> str:
    """Save content to file, compute hash, rename, and upload to S3."""
    # Save to temp file
    temp_path = os.path.join(TMP_DIR, f"{base_name}.{extension}")
    with open(temp_path, "w", encoding="utf-8") as f:
        f.write(content)

    # Compute MD5 hash
    with open(temp_path, "rb") as f:
        md5_hash = hashlib.md5(f.read()).hexdigest()

    # Rename with hash
    final_name = f"{md5_hash}.{extension}"
    final_path = os.path.join(TMP_DIR, final_name)
    os.rename(temp_path, final_path)

    # Upload to S3
    return upload_to_s3(final_path, final_name)


def get_param(request: Request, name: str, form_data: dict = None) -> Optional[str]:
    """Get parameter from query string or form data."""
    # Check query params first
    value = request.query_params.get(name)
    if value:
        return value

    # Check form data
    if form_data and name in form_data:
        return form_data[name]

    return None


# =============================================================================
# Routes
# =============================================================================

@app.get("/")
async def root():
    """Health check endpoint."""
    return {"message": "GO Tools API - Hello, we all love SGD!!"}


@app.api_route("/gotermfinder", methods=["GET", "POST"])
async def go_term_finder(
    request: Request,
    genes: Optional[str] = Query(None),
    aspect: Optional[str] = Query(None),
    genes4bg: Optional[str] = Query(None),
    evidence: Optional[str] = Query(None),
    pvalue: Optional[float] = Query(None),
    FDR: Optional[bool] = Query(None),
    file: Optional[str] = Query(None),
):
    """
    GO Term Finder endpoint - find enriched GO terms in a gene list.

    Returns HTML table with enriched terms and download URLs for various formats.
    """
    # Handle file download
    if file:
        filepath = os.path.join(TMP_DIR, file)
        if os.path.exists(filepath):
            return FileResponse(filepath, filename=file)
        return JSONResponse({"error": "File not found"}, status_code=404)

    # Get form data if POST
    form_data = {}
    if request.method == "POST":
        form_data = dict(await request.form())

    # Get parameters
    genes = genes or form_data.get("genes")
    aspect = aspect or form_data.get("aspect") or "F"
    genes4bg = genes4bg or form_data.get("genes4bg")
    evidence = evidence or form_data.get("evidence")
    pvalue = pvalue or float(form_data.get("pvalue", 0.01))
    FDR = FDR or form_data.get("FDR")

    if not genes:
        return JSONResponse({"error": "NO GENE NAME PASSED IN"})

    # Parse genes
    genes = genes.upper().replace("SGD:", "")
    gene_list = [g.strip() for g in genes.split("|") if g.strip()]

    # Parse background genes
    background_list = None
    if genes4bg:
        background_list = [g.strip() for g in genes4bg.split("|") if g.strip()]

    # Parse evidence codes
    evidence_codes = None
    if evidence:
        evidence = evidence.strip("|")
        evidence_codes = [e.strip() for e in evidence.replace("|", ",").split(",") if e.strip()]

    # Determine correction method
    correction = CorrectionMethod.BENJAMINI_HOCHBERG if FDR else CorrectionMethod.BONFERRONI

    # Run analysis
    result = term_finder_service.run_enrichment(
        query_genes=gene_list,
        aspect=aspect.upper(),
        background_genes=background_list,
        evidence_codes=evidence_codes,
        p_value_cutoff=pvalue,
        correction_method=correction,
    )

    # Check if no results
    if result.total_enriched_terms == 0:
        return JSONResponse({
            "output": "No significant GO terms were found for your input list of genes."
        })

    # Generate session ID for file names
    import random
    session_id = str(random.randint(1, 10000000))

    # Generate HTML
    html_content = generate_term_finder_html(result, aspect)

    # Generate tab-delimited file
    tab_content = generate_term_finder_tab(result)

    # Generate terms text file
    terms_content = generate_terms_text(result)

    # Save input genes
    input_content = "\n".join(gene_list)

    # Save files and get URLs
    html_url = save_and_upload(html_content, f"{session_id}", "html")
    tab_url = save_and_upload(tab_content, f"{session_id}_tab", "txt")
    term_url = save_and_upload(terms_content, f"{session_id}_terms", "txt")
    input_url = save_and_upload(input_content, session_id, "txt")

    # Clean HTML for embedding (remove html/body tags)
    clean_html = html_content.replace("<html><body>", "").replace("</body></html>", "")
    clean_html = clean_html.replace("color=red", "color=maroon")

    return JSONResponse({
        "html": clean_html,
        "image_html": "",  # No image generation in Python version
        "image_page": "",
        "tab_page": tab_url,
        "term_page": term_url,
        "table_page": html_url,
        "png_page": "",
        "svg_page": "",
        "ps_page": "",
        "input_page": input_url,
    })


@app.api_route("/termfinder", methods=["GET", "POST"])
async def term_finder_simple(
    request: Request,
    genes: Optional[str] = Query(None),
    aspect: Optional[str] = Query(None),
):
    """
    Simplified GO Term Finder endpoint - returns enriched terms as JSON array.

    Used for programmatic access to enrichment results.
    """
    # Get form data if POST
    form_data = {}
    if request.method == "POST":
        form_data = dict(await request.form())

    genes = genes or form_data.get("genes")
    aspect = aspect or form_data.get("aspect") or "P"

    if not genes:
        return JSONResponse([])

    # Parse genes
    genes = genes.upper().replace("SGD:", "")
    gene_list = [g.strip() for g in genes.split("|") if g.strip()]

    # Run analysis with FDR
    result = term_finder_service.run_enrichment(
        query_genes=gene_list,
        aspect=aspect.upper(),
        correction_method=CorrectionMethod.BENJAMINI_HOCHBERG,
    )

    # Format response
    data = []
    for term in result.enriched_terms:
        pvalue_str = format_pvalue(term.p_value)
        data.append({
            "goid": term.go_id,
            "term": term.go_term,
            "pvalue": pvalue_str,
            "num_gene_annotated": str(term.query_count),
        })

    return JSONResponse(data)


@app.api_route("/goslimmapper", methods=["GET", "POST"])
async def go_slim_mapper(
    request: Request,
    genes: Optional[str] = Query(None),
    aspect: Optional[str] = Query(None),
    terms: Optional[str] = Query(None),
    file: Optional[str] = Query(None),
):
    """
    GO Slim Mapper endpoint - map genes to GO Slim categories.

    Returns HTML table with mapped slim terms and download URLs.
    """
    # Handle file download
    if file:
        filepath = os.path.join(TMP_DIR, file)
        if os.path.exists(filepath):
            return FileResponse(filepath, filename=file)
        return JSONResponse({"error": "File not found"}, status_code=404)

    # Get form data if POST
    form_data = {}
    if request.method == "POST":
        form_data = dict(await request.form())

    genes = genes or form_data.get("genes")
    aspect = aspect or form_data.get("aspect")
    terms = terms or form_data.get("terms")

    if not genes:
        return JSONResponse({"error": "NO GENE NAME PASSED IN"})

    if not aspect:
        return JSONResponse({"error": "NO GO ASPECT PASSED IN"})

    if not terms:
        return JSONResponse({"error": "NO SLIM TERMS PASSED IN"})

    # Parse genes
    genes = genes.upper().replace("SGD:", "")
    gene_list = [g.strip() for g in genes.split("|") if g.strip()]

    # Parse slim terms
    slim_terms = parse_slim_terms_from_string(terms, ontology)

    # Run analysis
    result = slim_mapper_service.run_mapping(
        query_genes=gene_list,
        slim_terms=slim_terms,
        aspect=aspect.upper(),
    )

    # Generate session ID
    import random
    session_id = str(random.randint(1, 10000000))
    base_name = f"mapper_genes_{session_id}"

    # Generate HTML
    html_content = generate_slim_mapper_html(result, aspect)

    # Generate tab-delimited file
    tab_content = generate_slim_mapper_tab(result)

    # Generate terms text file
    terms_content = generate_slim_terms_text(result)

    # Save input genes
    input_content = "\n".join(gene_list)

    # Save slim input
    slim_input = terms

    # Save files and get URLs
    html_url = save_and_upload(html_content, f"{base_name}_slimTerms", "html")
    tab_url = save_and_upload(tab_content, f"{base_name}_slimTab", "txt")
    term_url = save_and_upload(terms_content, f"{base_name}_slimTerms", "txt")
    gene_input_url = save_and_upload(input_content, base_name, "txt")
    slim_input_url = save_and_upload(slim_input, f"mapper_terms_{session_id}", "txt")

    # Clean HTML for embedding
    clean_html = html_content.replace("<html><body>", "").replace("</body></html>", "")
    clean_html = clean_html.replace("<br><b>", "").replace("</b><br><br><center>", "<center>")
    clean_html = clean_html.replace("color=red", "color=maroon")
    clean_html = clean_html.replace("<td colspan=5>", "<td colspan=5 bgcolor='#FFCC99'>")
    clean_html = clean_html.replace("<font color=#FFFFFF>", "").replace("</font>", "")
    clean_html = clean_html.replace("<th align=center nowrap>", "<th bgcolor='#CCCCFF' align=center nowrap>")
    clean_html = clean_html.replace("<th align=center>", "<th bgcolor='#CCCCFF' align=center nowrap>")
    clean_html = clean_html.replace('<tr bgcolor="FFE4C4">', '')
    clean_html = clean_html.replace(' nowrap=""', '')
    clean_html = clean_html.replace("( ", "(").replace(" )", ")")
    clean_html = clean_html.replace("infowin", "_extwin")

    return JSONResponse({
        "html": clean_html,
        "table_page": html_url,
        "tab_page": tab_url,
        "term_page": term_url,
        "gene_input_page": gene_input_url,
        "slim_input_page": slim_input_url,
    })


@app.get("/download/{filename}")
async def download_file(filename: str):
    """Download a generated file."""
    filepath = os.path.join(TMP_DIR, filename)
    if os.path.exists(filepath):
        return FileResponse(filepath, filename=filename)
    return JSONResponse({"error": "File not found"}, status_code=404)


# =============================================================================
# Helper Functions for Text Output
# =============================================================================

def generate_terms_text(result) -> str:
    """Generate plain text terms output."""
    lines = [
        f"GO Term Finder Results",
        f"Query genes submitted: {result.query_genes_submitted}",
        f"Query genes found: {result.query_genes_found}",
        f"Query genes with GO: {result.query_genes_with_go}",
        f"Background size: {result.background_size}",
        f"Total enriched terms: {result.total_enriched_terms}",
        "",
        "Enriched Terms:",
        "-" * 80,
    ]

    for term in result.enriched_terms:
        lines.append(f"\n{term.go_id}\t{term.go_term}")
        lines.append(f"  Aspect: {term.aspect_name}")
        lines.append(f"  P-value: {term.p_value:.2e}")
        if term.fdr is not None:
            lines.append(f"  FDR: {term.fdr:.2e}")
        lines.append(f"  Query: {term.query_count}/{term.query_total} ({term.query_frequency}%)")
        lines.append(f"  Background: {term.background_count}/{term.background_total} ({term.background_frequency}%)")
        lines.append(f"  Fold enrichment: {term.fold_enrichment}")
        lines.append(f"  Genes: {', '.join(g.gene_symbol for g in term.genes)}")

    return "\n".join(lines)


def generate_slim_terms_text(result) -> str:
    """Generate plain text slim mapper output."""
    lines = [
        f"GO Slim Mapper Results",
        f"Query genes submitted: {result.query_genes_submitted}",
        f"Query genes found: {result.query_genes_found}",
        f"Query genes with GO: {result.query_genes_with_go}",
        "",
        "Mapped Slim Terms:",
        "-" * 80,
    ]

    for term in result.mapped_terms:
        lines.append(f"\n{term.go_id}\t{term.go_term}")
        lines.append(f"  Gene count: {term.gene_count}/{term.total_genes} ({term.frequency_percent}%)")
        lines.append(f"  Genes: {', '.join(g.gene_symbol for g in term.genes)}")

    if result.other_genes:
        lines.append("\n\nOther genes (have GO but not mapped to slim):")
        lines.append(", ".join(g.gene_symbol for g in result.other_genes))

    if result.not_annotated_genes:
        lines.append("\n\nNot annotated genes (no GO annotations):")
        lines.append(", ".join(g.gene_symbol for g in result.not_annotated_genes))

    return "\n".join(lines)


# =============================================================================
# Main Entry Point
# =============================================================================

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
