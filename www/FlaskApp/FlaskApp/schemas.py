"""
Pydantic schemas for GO Tools API requests and responses.
"""
from __future__ import annotations

from enum import Enum
from typing import List, Optional

from pydantic import BaseModel, Field


class CorrectionMethod(str, Enum):
    """Multiple testing correction methods."""
    NONE = "none"
    BONFERRONI = "bonferroni"
    BENJAMINI_HOCHBERG = "benjamini_hochberg"


class GOAspect(str, Enum):
    """GO ontology aspects."""
    PROCESS = "P"
    FUNCTION = "F"
    COMPONENT = "C"


# =============================================================================
# GO Term Finder Schemas
# =============================================================================

class GOTermFinderRequest(BaseModel):
    """Request for GO Term Finder analysis."""
    genes: str = Field(..., description="Pipe-separated list of gene names")
    aspect: GOAspect = Field(default=GOAspect.FUNCTION, description="GO aspect (P, F, or C)")
    genes4bg: Optional[str] = Field(None, description="Pipe-separated background genes")
    evidence: Optional[str] = Field(None, description="Pipe-separated evidence codes")
    pvalue: float = Field(default=0.01, description="P-value cutoff")
    FDR: Optional[bool] = Field(None, description="Calculate FDR")


class GeneHitResponse(BaseModel):
    """A gene contributing to an enriched term."""
    gene_id: str
    gene_symbol: str
    gene_name: str
    evidence_codes: List[str]


class EnrichedTermResponse(BaseModel):
    """An enriched GO term."""
    goid: str
    term: str
    aspect: str
    aspect_name: str
    query_count: int
    query_total: int
    background_count: int
    background_total: int
    query_frequency: float
    background_frequency: float
    fold_enrichment: float
    pvalue: str
    fdr: Optional[str]
    genes: List[GeneHitResponse]


class GOTermFinderResponse(BaseModel):
    """Response for GO Term Finder analysis."""
    html: Optional[str] = None
    enriched_terms: Optional[List[EnrichedTermResponse]] = None
    query_genes_submitted: int = 0
    query_genes_found: int = 0
    query_genes_with_go: int = 0
    query_genes_not_found: List[str] = []
    background_size: int = 0
    total_enriched_terms: int = 0
    output: Optional[str] = None
    error: Optional[str] = None


class EnrichmentTermResponse(BaseModel):
    """Simplified enrichment term for /termfinder endpoint."""
    goid: str
    term: str
    pvalue: str
    num_gene_annotated: str


# =============================================================================
# GO Slim Mapper Schemas
# =============================================================================

class GOSlimMapperRequest(BaseModel):
    """Request for GO Slim Mapper analysis."""
    genes: str = Field(..., description="Pipe-separated list of gene names")
    aspect: GOAspect = Field(..., description="GO aspect (P, F, or C)")
    terms: str = Field(..., description="Pipe-separated slim terms")


class MappedGeneResponse(BaseModel):
    """A gene mapped to a slim term."""
    gene_id: str
    gene_symbol: str
    gene_name: str


class MappedSlimTermResponse(BaseModel):
    """A GO Slim term with mapped genes."""
    goid: str
    term: str
    aspect: str
    gene_count: int
    total_genes: int
    frequency_percent: float
    genes: List[MappedGeneResponse]


class GOSlimMapperResponse(BaseModel):
    """Response for GO Slim Mapper analysis."""
    html: Optional[str] = None
    mapped_terms: Optional[List[MappedSlimTermResponse]] = None
    query_genes_submitted: int = 0
    query_genes_found: int = 0
    query_genes_with_go: int = 0
    query_genes_not_found: List[str] = []
    other_genes: List[MappedGeneResponse] = []
    not_annotated_genes: List[MappedGeneResponse] = []
    error: Optional[str] = None
