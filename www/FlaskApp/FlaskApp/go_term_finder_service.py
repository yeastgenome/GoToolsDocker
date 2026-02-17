"""
GO Term Finder Service - Gene Ontology enrichment analysis.

Performs hypergeometric test for GO term enrichment with optional
multiple testing correction (Bonferroni or Benjamini-Hochberg FDR).
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Set, Tuple

from scipy.stats import hypergeom

from go_parser import (
    GOOntology,
    GeneAssociations,
    ASPECT_NAMES,
)


class CorrectionMethod(str, Enum):
    NONE = "none"
    BONFERRONI = "bonferroni"
    BENJAMINI_HOCHBERG = "benjamini_hochberg"


@dataclass
class GeneHit:
    """A gene that contributes to an enriched term."""
    db_object_id: str
    gene_symbol: str
    gene_name: str
    evidence_codes: List[str]


@dataclass
class EnrichedTerm:
    """An enriched GO term from the analysis."""
    go_id: str
    go_term: str
    go_aspect: str
    aspect_name: str
    query_count: int      # k - genes in query with this term
    query_total: int      # n - total query genes
    background_count: int  # K - genes in background with this term
    background_total: int  # N - total background genes
    query_frequency: float
    background_frequency: float
    fold_enrichment: float
    p_value: float
    fdr: Optional[float]
    genes: List[GeneHit]


@dataclass
class EnrichmentResult:
    """Complete result of GO Term Finder analysis."""
    query_genes_submitted: int
    query_genes_found: int
    query_genes_with_go: int
    query_genes_not_found: List[str]
    background_size: int
    background_type: str
    aspect: str
    evidence_codes_used: List[str]
    p_value_cutoff: float
    correction_method: str
    enriched_terms: List[EnrichedTerm]
    total_enriched_terms: int


class GOTermFinderService:
    """
    Service for GO term enrichment analysis using hypergeometric test.
    """

    def __init__(self, ontology: GOOntology, associations: GeneAssociations):
        self.ontology = ontology
        self.associations = associations

    def run_enrichment(
        self,
        query_genes: List[str],
        aspect: str = "P",
        background_genes: Optional[List[str]] = None,
        evidence_codes: Optional[List[str]] = None,
        p_value_cutoff: float = 0.01,
        correction_method: CorrectionMethod = CorrectionMethod.BENJAMINI_HOCHBERG,
        min_genes_in_term: int = 1,
    ) -> EnrichmentResult:
        """
        Run GO term enrichment analysis.

        Args:
            query_genes: List of gene names/IDs to analyze
            aspect: GO aspect (P, F, or C)
            background_genes: Optional custom background gene list
            evidence_codes: Optional list of evidence codes to include
            p_value_cutoff: P-value cutoff for significance
            correction_method: Multiple testing correction method
            min_genes_in_term: Minimum genes required per term

        Returns:
            EnrichmentResult with enriched terms
        """
        evidence_set = set(evidence_codes) if evidence_codes else None

        # Step 1: Validate query genes
        query_gene_ids, not_found = self._validate_genes(query_genes)

        if not query_gene_ids:
            return EnrichmentResult(
                query_genes_submitted=len(query_genes),
                query_genes_found=0,
                query_genes_with_go=0,
                query_genes_not_found=not_found,
                background_size=0,
                background_type="default",
                aspect=aspect,
                evidence_codes_used=evidence_codes or [],
                p_value_cutoff=p_value_cutoff,
                correction_method=correction_method.value,
                enriched_terms=[],
                total_enriched_terms=0,
            )

        # Step 2: Build background set
        if background_genes:
            background_type = "custom"
            background_gene_ids, _ = self._validate_genes(background_genes)
        else:
            background_type = "default"
            background_gene_ids = list(self.associations.get_all_annotated_genes(
                aspect=aspect,
                evidence_codes=evidence_set,
            ))

        # Ensure query genes are in background
        background_set = set(background_gene_ids)
        query_gene_ids = [g for g in query_gene_ids if g in background_set]

        if not query_gene_ids:
            return EnrichmentResult(
                query_genes_submitted=len(query_genes),
                query_genes_found=len(query_gene_ids),
                query_genes_with_go=0,
                query_genes_not_found=not_found,
                background_size=len(background_gene_ids),
                background_type=background_type,
                aspect=aspect,
                evidence_codes_used=evidence_codes or [],
                p_value_cutoff=p_value_cutoff,
                correction_method=correction_method.value,
                enriched_terms=[],
                total_enriched_terms=0,
            )

        # Step 3: Get GO annotations with ancestors
        query_annotations = self._get_annotations_with_ancestors(
            query_gene_ids, aspect, evidence_set
        )
        background_annotations = self._get_annotations_with_ancestors(
            background_gene_ids, aspect, evidence_set
        )

        # Filter to genes with GO annotations
        query_genes_with_go = [g for g in query_gene_ids if g in query_annotations]

        if not query_genes_with_go:
            return EnrichmentResult(
                query_genes_submitted=len(query_genes),
                query_genes_found=len(query_gene_ids),
                query_genes_with_go=0,
                query_genes_not_found=not_found,
                background_size=len(background_annotations),
                background_type=background_type,
                aspect=aspect,
                evidence_codes_used=evidence_codes or [],
                p_value_cutoff=p_value_cutoff,
                correction_method=correction_method.value,
                enriched_terms=[],
                total_enriched_terms=0,
            )

        # Step 4: Calculate enrichment
        enrichment_results = self._calculate_enrichment(
            {g: query_annotations[g] for g in query_genes_with_go},
            background_annotations,
            p_value_cutoff,
            min_genes_in_term,
        )

        # Step 5: Apply multiple testing correction
        corrected_results = self._apply_correction(
            enrichment_results,
            correction_method,
            p_value_cutoff,
        )

        # Step 6: Build enriched term objects
        enriched_terms = self._build_enriched_terms(
            corrected_results,
            query_annotations,
            query_genes_with_go,
            aspect,
            evidence_set,
        )

        # Sort by p-value
        enriched_terms.sort(key=lambda x: x.p_value)

        return EnrichmentResult(
            query_genes_submitted=len(query_genes),
            query_genes_found=len(query_gene_ids) + len([g for g in query_genes if self.associations.lookup_gene(g)]),
            query_genes_with_go=len(query_genes_with_go),
            query_genes_not_found=not_found,
            background_size=len(background_annotations),
            background_type=background_type,
            aspect=aspect,
            evidence_codes_used=evidence_codes or [],
            p_value_cutoff=p_value_cutoff,
            correction_method=correction_method.value,
            enriched_terms=enriched_terms,
            total_enriched_terms=len(enriched_terms),
        )

    def _validate_genes(self, genes: List[str]) -> Tuple[List[str], List[str]]:
        """Validate gene names and return (found_ids, not_found_names)."""
        found_ids = []
        not_found = []
        seen = set()

        for gene in genes:
            gene = gene.strip()
            if not gene:
                continue

            gene_id = self.associations.lookup_gene(gene)
            if gene_id and gene_id not in seen:
                found_ids.append(gene_id)
                seen.add(gene_id)
            elif gene_id is None:
                not_found.append(gene)

        return found_ids, not_found

    def _get_annotations_with_ancestors(
        self,
        gene_ids: List[str],
        aspect: str,
        evidence_codes: Optional[Set[str]],
    ) -> Dict[str, Set[str]]:
        """Get GO annotations with ancestor terms for each gene."""
        result = {}
        for gene_id in gene_ids:
            go_ids = self.associations.get_gene_go_annotations_with_ancestors(
                gene_id, aspect, evidence_codes
            )
            if go_ids:
                result[gene_id] = go_ids
        return result

    def _calculate_enrichment(
        self,
        query_annotations: Dict[str, Set[str]],
        background_annotations: Dict[str, Set[str]],
        p_value_cutoff: float,
        min_genes_in_term: int,
    ) -> List[Tuple[str, int, int, int, int, float]]:
        """
        Calculate enrichment using hypergeometric test.

        P(X >= k) = hypergeom.sf(k-1, N, K, n)

        Where:
        - N = background set size
        - K = genes in background annotated to term
        - n = query set size
        - k = genes in query annotated to term

        Returns list of (go_id, k, n, K, N, p_value) tuples.
        """
        N = len(background_annotations)
        n = len(query_annotations)

        if N == 0 or n == 0:
            return []

        # Count genes per GO term
        query_term_counts: Dict[str, Set[str]] = defaultdict(set)
        background_term_counts: Dict[str, Set[str]] = defaultdict(set)

        for gene_id, go_ids in query_annotations.items():
            for go_id in go_ids:
                query_term_counts[go_id].add(gene_id)

        for gene_id, go_ids in background_annotations.items():
            for go_id in go_ids:
                background_term_counts[go_id].add(gene_id)

        # Calculate p-values
        results = []
        for go_id, query_genes in query_term_counts.items():
            k = len(query_genes)
            K = len(background_term_counts.get(go_id, set()))

            if k < min_genes_in_term or K == 0:
                continue

            # Hypergeometric test: P(X >= k)
            p_value = hypergeom.sf(k - 1, N, K, n)

            if p_value <= p_value_cutoff:
                results.append((go_id, k, n, K, N, p_value))

        return results

    def _apply_correction(
        self,
        results: List[Tuple[str, int, int, int, int, float]],
        method: CorrectionMethod,
        p_value_cutoff: float,
    ) -> List[Tuple[str, int, int, int, int, float, Optional[float]]]:
        """Apply multiple testing correction."""
        if not results:
            return []

        if method == CorrectionMethod.NONE:
            return [(go_id, k, n, K, N, p_val, None)
                    for go_id, k, n, K, N, p_val in results]

        n_tests = len(results)

        if method == CorrectionMethod.BONFERRONI:
            corrected = []
            for go_id, k, n, K, N, p_val in results:
                corrected_p = min(p_val * n_tests, 1.0)
                if corrected_p <= p_value_cutoff:
                    corrected.append((go_id, k, n, K, N, p_val, corrected_p))
            return corrected

        elif method == CorrectionMethod.BENJAMINI_HOCHBERG:
            # Sort by p-value
            sorted_results = sorted(results, key=lambda x: x[5])

            # Calculate FDR for each rank
            fdr_values = []
            for i, (go_id, k, n, K, N, p_val) in enumerate(sorted_results):
                rank = i + 1
                fdr = (p_val * n_tests) / rank
                fdr_values.append((go_id, k, n, K, N, p_val, fdr))

            # Enforce monotonicity (FDR can only decrease as rank increases)
            for i in range(len(fdr_values) - 2, -1, -1):
                go_id, k, n, K, N, p_val, fdr = fdr_values[i]
                next_fdr = fdr_values[i + 1][6]
                if fdr > next_fdr:
                    fdr_values[i] = (go_id, k, n, K, N, p_val, next_fdr)

            # Cap FDR at 1.0 and filter by cutoff
            corrected = []
            for go_id, k, n, K, N, p_val, fdr in fdr_values:
                fdr = min(fdr, 1.0)
                if fdr <= p_value_cutoff:
                    corrected.append((go_id, k, n, K, N, p_val, fdr))

            return corrected

        return [(go_id, k, n, K, N, p_val, None)
                for go_id, k, n, K, N, p_val in results]

    def _build_enriched_terms(
        self,
        corrected_results: List[Tuple[str, int, int, int, int, float, Optional[float]]],
        query_annotations: Dict[str, Set[str]],
        query_gene_ids: List[str],
        aspect: str,
        evidence_codes: Optional[Set[str]],
    ) -> List[EnrichedTerm]:
        """Build EnrichedTerm objects from results."""
        enriched_terms = []

        for go_id, k, n, K, N, p_val, fdr in corrected_results:
            term = self.ontology.get_term(go_id)
            if not term:
                continue

            # Build gene hits
            gene_hits = []
            for gene_id in query_gene_ids:
                if go_id in query_annotations.get(gene_id, set()):
                    # Get evidence codes for this gene-term pair
                    ev_codes = []
                    for ann_go_id, ev_code, ann_aspect in self.associations.id_to_annotations.get(gene_id, set()):
                        if ann_go_id == go_id:
                            ev_codes.append(ev_code)

                    gene_hits.append(GeneHit(
                        db_object_id=gene_id,
                        gene_symbol=self.associations.get_gene_symbol(gene_id),
                        gene_name=self.associations.get_gene_name(gene_id),
                        evidence_codes=ev_codes,
                    ))

            # Calculate frequencies
            query_frequency = (k / n) * 100 if n > 0 else 0.0
            background_frequency = (K / N) * 100 if N > 0 else 0.0
            fold_enrichment = (k / n) / (K / N) if K > 0 and N > 0 and n > 0 else 0.0

            enriched_terms.append(EnrichedTerm(
                go_id=self.ontology.format_goid(go_id),
                go_term=term.name,
                go_aspect=term.aspect or aspect,
                aspect_name=ASPECT_NAMES.get(term.aspect or aspect, ""),
                query_count=k,
                query_total=n,
                background_count=K,
                background_total=N,
                query_frequency=round(query_frequency, 2),
                background_frequency=round(background_frequency, 4),
                fold_enrichment=round(fold_enrichment, 2),
                p_value=p_val,
                fdr=fdr,
                genes=gene_hits,
            ))

        return enriched_terms


def format_pvalue(p_value: float) -> str:
    """Format p-value for display."""
    if p_value == 0:
        return "0"
    if "e-" in f"{p_value}":
        # Scientific notation
        parts = f"{p_value}".split(".")
        if len(parts) == 2:
            decimal_exp = parts[1].split("e-")
            if len(decimal_exp) == 2:
                return f"{parts[0]}.{decimal_exp[0][:2]}e-{decimal_exp[1]}"
        return f"{p_value:.2e}"
    else:
        return f"{p_value:.5f}"
