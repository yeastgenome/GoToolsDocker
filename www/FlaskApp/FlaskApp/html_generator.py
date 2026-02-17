"""
HTML Generator - Generate HTML output for GO Term Finder and GO Slim Mapper.

Produces HTML tables compatible with the original Perl output format.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from go_term_finder_service import EnrichmentResult
    from go_slim_mapper_service import SlimMapperResult

# SGD base URL for gene links
SGD_LOCUS_URL = "https://www.yeastgenome.org/locus/"
GO_TERM_URL = "https://www.yeastgenome.org/go/"


def generate_term_finder_html(result: "EnrichmentResult", aspect: str) -> str:
    """
    Generate HTML table for GO Term Finder results.

    Args:
        result: EnrichmentResult from the analysis
        aspect: GO aspect (P, F, or C)

    Returns:
        HTML string
    """
    aspect_names = {
        "P": "Biological Process",
        "F": "Molecular Function",
        "C": "Cellular Component",
    }
    aspect_name = aspect_names.get(aspect.upper(), aspect)

    html_parts = ["<html><body>"]

    # Summary header
    html_parts.append(f"<b>GO Term Finder Results - {aspect_name}</b><br><br>")
    html_parts.append("<table border=0 cellpadding=3>")
    html_parts.append(f"<tr><td>Query genes submitted:</td><td><b>{result.query_genes_submitted}</b></td></tr>")
    html_parts.append(f"<tr><td>Query genes found:</td><td><b>{result.query_genes_found}</b></td></tr>")
    html_parts.append(f"<tr><td>Query genes with GO annotation:</td><td><b>{result.query_genes_with_go}</b></td></tr>")
    html_parts.append(f"<tr><td>Background genes:</td><td><b>{result.background_size}</b></td></tr>")
    html_parts.append(f"<tr><td>Total enriched terms:</td><td><b>{result.total_enriched_terms}</b></td></tr>")
    html_parts.append("</table><br><br>")

    if not result.enriched_terms:
        html_parts.append("<p>No significantly enriched GO terms were found.</p>")
        html_parts.append("</body></html>")
        return "".join(html_parts)

    # Not found genes
    if result.query_genes_not_found:
        html_parts.append("<b>Genes not found in database:</b><br>")
        html_parts.append(", ".join(result.query_genes_not_found))
        html_parts.append("<br><br>")

    # Results table
    html_parts.append('<a name="table" />')
    html_parts.append('<center><table border=1 cellpadding=5>')

    # Header row
    html_parts.append('<tr bgcolor="#CCCCFF">')
    html_parts.append('<th align=center nowrap>GO ID</th>')
    html_parts.append('<th align=center>Term</th>')
    html_parts.append('<th align=center nowrap>Cluster<br>frequency</th>')
    html_parts.append('<th align=center nowrap>Background<br>frequency</th>')
    html_parts.append('<th align=center nowrap>P-value</th>')
    if result.correction_method == "benjamini_hochberg":
        html_parts.append('<th align=center nowrap>FDR</th>')
    html_parts.append('<th align=center>Genes</th>')
    html_parts.append('</tr>')

    # Data rows
    for term in result.enriched_terms:
        html_parts.append('<tr>')

        # GO ID with link
        go_url = f"{GO_TERM_URL}{term.go_id}"
        html_parts.append(f'<td nowrap><a href="{go_url}" target="_extwin">{term.go_id}</a></td>')

        # Term name
        html_parts.append(f'<td>{term.go_term}</td>')

        # Cluster frequency
        html_parts.append(f'<td align=center nowrap>{term.query_count} of {term.query_total}<br>({term.query_frequency}%)</td>')

        # Background frequency
        html_parts.append(f'<td align=center nowrap>{term.background_count} of {term.background_total}<br>({term.background_frequency}%)</td>')

        # P-value
        pval_str = f"{term.p_value:.2e}" if term.p_value < 0.01 else f"{term.p_value:.4f}"
        html_parts.append(f'<td align=center nowrap>{pval_str}</td>')

        # FDR
        if result.correction_method == "benjamini_hochberg" and term.fdr is not None:
            fdr_str = f"{term.fdr:.2e}" if term.fdr < 0.01 else f"{term.fdr:.4f}"
            html_parts.append(f'<td align=center nowrap>{fdr_str}</td>')

        # Genes with links
        gene_links = []
        for gene in term.genes:
            gene_url = f"{SGD_LOCUS_URL}{gene.db_object_id}"
            gene_links.append(f'<a href="{gene_url}" target="_extwin">{gene.gene_symbol}</a>')
        html_parts.append(f'<td>{", ".join(gene_links)}</td>')

        html_parts.append('</tr>')

    html_parts.append('</table></center>')
    html_parts.append("</body></html>")

    return "".join(html_parts)


def generate_term_finder_tab(result: "EnrichmentResult") -> str:
    """
    Generate tab-delimited output for GO Term Finder results.

    Args:
        result: EnrichmentResult from the analysis

    Returns:
        Tab-delimited string
    """
    lines = []

    # Header
    headers = [
        "GOID",
        "Term",
        "P-value",
        "FDR",
        "Cluster frequency",
        "Background frequency",
        "Genes"
    ]
    lines.append("\t".join(headers))

    # Data rows
    for term in result.enriched_terms:
        fdr_str = f"{term.fdr:.6e}" if term.fdr is not None else ""
        genes_str = ", ".join(g.gene_symbol for g in term.genes)

        row = [
            term.go_id,
            term.go_term,
            f"{term.p_value:.6e}",
            fdr_str,
            f"{term.query_count} of {term.query_total}",
            f"{term.background_count} of {term.background_total}",
            genes_str,
        ]
        lines.append("\t".join(row))

    return "\n".join(lines)


def generate_slim_mapper_html(result: "SlimMapperResult", aspect: str) -> str:
    """
    Generate HTML table for GO Slim Mapper results.

    Args:
        result: SlimMapperResult from the analysis
        aspect: GO aspect (P, F, or C)

    Returns:
        HTML string
    """
    aspect_names = {
        "P": "Biological Process",
        "F": "Molecular Function",
        "C": "Cellular Component",
    }
    aspect_name = aspect_names.get(aspect.upper(), aspect)

    html_parts = ["<html><body>"]

    # Summary header
    html_parts.append(f"<b>GO Slim Mapper Results - {aspect_name}</b><br><br>")

    # Summary table
    html_parts.append("<table border=0 cellpadding=3>")
    html_parts.append(f"<tr><td>Query genes submitted:</td><td><b>{result.query_genes_submitted}</b></td></tr>")
    html_parts.append(f"<tr><td>Query genes found:</td><td><b>{result.query_genes_found}</b></td></tr>")
    html_parts.append(f"<tr><td>Query genes with GO annotation:</td><td><b>{result.query_genes_with_go}</b></td></tr>")
    html_parts.append("</table><br><br>")

    if not result.mapped_terms:
        html_parts.append("<p>No genes were mapped to the specified GO Slim terms.</p>")
        html_parts.append("</body></html>")
        return "".join(html_parts)

    # Not found genes
    if result.query_genes_not_found:
        html_parts.append("<b>Genes not found in database:</b><br>")
        html_parts.append(", ".join(result.query_genes_not_found))
        html_parts.append("<br><br>")

    # Results table
    html_parts.append('<center><table border=1 cellpadding=5>')

    # Header row
    html_parts.append('<tr bgcolor="#CCCCFF">')
    html_parts.append('<th align=center nowrap>GO Slim Term</th>')
    html_parts.append('<th align=center nowrap>GO ID</th>')
    html_parts.append('<th align=center nowrap>Gene Count</th>')
    html_parts.append('<th align=center nowrap>Frequency</th>')
    html_parts.append('<th align=center>Genes</th>')
    html_parts.append('</tr>')

    # Data rows
    for term in result.mapped_terms:
        html_parts.append('<tr>')

        # Term name
        html_parts.append(f'<td>{term.go_term}</td>')

        # GO ID with link
        go_url = f"{GO_TERM_URL}{term.go_id}"
        html_parts.append(f'<td nowrap><a href="{go_url}" target="_extwin">{term.go_id}</a></td>')

        # Gene count
        html_parts.append(f'<td align=center>{term.gene_count}</td>')

        # Frequency
        html_parts.append(f'<td align=center>{term.frequency_percent}%</td>')

        # Genes with links
        gene_links = []
        for gene in term.genes:
            gene_url = f"{SGD_LOCUS_URL}{gene.db_object_id}"
            gene_links.append(f'<a href="{gene_url}" target="_extwin">{gene.gene_symbol}</a>')
        html_parts.append(f'<td>{", ".join(gene_links)}</td>')

        html_parts.append('</tr>')

    html_parts.append('</table></center>')

    # Other genes section
    if result.other_genes:
        html_parts.append('<br><br><td colspan=5 bgcolor="#FFCC99"><b>Other genes (have GO annotations but not mapped to slim):</b></td><br>')
        gene_links = []
        for gene in result.other_genes:
            gene_url = f"{SGD_LOCUS_URL}{gene.db_object_id}"
            gene_links.append(f'<a href="{gene_url}" target="_extwin">{gene.gene_symbol}</a>')
        html_parts.append(", ".join(gene_links))

    # Not annotated genes section
    if result.not_annotated_genes:
        html_parts.append('<br><br><td colspan=5 bgcolor="#FFCC99"><b>Genes without GO annotations:</b></td><br>')
        gene_links = []
        for gene in result.not_annotated_genes:
            gene_url = f"{SGD_LOCUS_URL}{gene.db_object_id}"
            gene_links.append(f'<a href="{gene_url}" target="_extwin">{gene.gene_symbol}</a>')
        html_parts.append(", ".join(gene_links))

    html_parts.append("</body></html>")

    return "".join(html_parts)


def generate_slim_mapper_tab(result: "SlimMapperResult") -> str:
    """
    Generate tab-delimited output for GO Slim Mapper results.

    Args:
        result: SlimMapperResult from the analysis

    Returns:
        Tab-delimited string
    """
    lines = []

    # Header
    headers = [
        "GO Slim Term",
        "GOID",
        "Gene Count",
        "Frequency (%)",
        "Genes"
    ]
    lines.append("\t".join(headers))

    # Data rows
    for term in result.mapped_terms:
        genes_str = ", ".join(g.gene_symbol for g in term.genes)

        row = [
            term.go_term,
            term.go_id,
            str(term.gene_count),
            f"{term.frequency_percent}",
            genes_str,
        ]
        lines.append("\t".join(row))

    return "\n".join(lines)
