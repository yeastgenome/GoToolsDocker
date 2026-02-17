"""
GO Slim Mapper Service - Map genes to predefined GO Slim categories.

Maps a list of genes to broader GO Slim terms via direct or indirect
(ancestor) annotations, without statistical enrichment analysis.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Set

from go_parser import (
    GOOntology,
    GeneAssociations,
    load_slim_terms,
    ASPECT_NAMES,
)


@dataclass
class MappedGene:
    """A gene that maps to a slim term."""
    db_object_id: str
    gene_symbol: str
    gene_name: str


@dataclass
class MappedSlimTerm:
    """A GO Slim term with its mapped genes."""
    go_id: str
    go_term: str
    go_aspect: str
    gene_count: int
    total_genes: int
    frequency_percent: float
    genes: List[MappedGene]


@dataclass
class SlimMapperResult:
    """Complete result of GO Slim Mapper analysis."""
    query_genes_submitted: int
    query_genes_found: int
    query_genes_with_go: int
    query_genes_not_found: List[str]
    aspect: str
    mapped_terms: List[MappedSlimTerm]
    other_genes: List[MappedGene]  # Have GO but not mapped to slim
    not_annotated_genes: List[MappedGene]  # No GO annotations


class GOSlimMapperService:
    """
    Service for mapping genes to GO Slim categories.
    """

    def __init__(self, ontology: GOOntology, associations: GeneAssociations):
        self.ontology = ontology
        self.associations = associations

    def run_mapping(
        self,
        query_genes: List[str],
        slim_terms: Set[str],
        aspect: str = "P",
        evidence_codes: Optional[List[str]] = None,
    ) -> SlimMapperResult:
        """
        Map genes to GO Slim terms.

        Args:
            query_genes: List of gene names/IDs to map
            slim_terms: Set of GO IDs representing the slim
            aspect: GO aspect (P, F, or C)
            evidence_codes: Optional list of evidence codes to include

        Returns:
            SlimMapperResult with mapped terms and genes
        """
        evidence_set = set(evidence_codes) if evidence_codes else None

        # Step 1: Validate genes
        gene_map, not_found = self._validate_genes(query_genes)

        if not gene_map:
            return SlimMapperResult(
                query_genes_submitted=len(query_genes),
                query_genes_found=0,
                query_genes_with_go=0,
                query_genes_not_found=not_found,
                aspect=aspect,
                mapped_terms=[],
                other_genes=[],
                not_annotated_genes=[],
            )

        # Step 2: Normalize slim terms
        normalized_slim = set()
        for term in slim_terms:
            normalized = self.ontology.normalize_go_id(term)
            if normalized:
                normalized_slim.add(normalized)

        if not normalized_slim:
            return SlimMapperResult(
                query_genes_submitted=len(query_genes),
                query_genes_found=len(gene_map),
                query_genes_with_go=0,
                query_genes_not_found=not_found,
                aspect=aspect,
                mapped_terms=[],
                other_genes=[],
                not_annotated_genes=list(gene_map.values()),
            )

        # Step 3: Get GO annotations for genes
        feature_to_go_ids: Dict[str, Set[str]] = {}
        for gene_id in gene_map.keys():
            go_ids = self.associations.get_gene_go_annotations(
                gene_id, aspect, evidence_set
            )
            if go_ids:
                feature_to_go_ids[gene_id] = go_ids

        # Step 4: Build ancestor map for all annotations
        all_go_ids = set()
        for go_ids in feature_to_go_ids.values():
            all_go_ids.update(go_ids)

        go_to_ancestors: Dict[str, Set[str]] = {}
        for go_id in all_go_ids:
            go_to_ancestors[go_id] = self.ontology.get_ancestors(go_id)

        # Step 5: Map genes to slim terms
        slim_to_genes: Dict[str, Set[str]] = defaultdict(set)
        genes_with_go = set()
        genes_mapped_to_slim = set()

        for gene_id, go_ids in feature_to_go_ids.items():
            genes_with_go.add(gene_id)
            mapped = False

            for go_id in go_ids:
                # Check if direct annotation is a slim term
                if go_id in normalized_slim:
                    slim_to_genes[go_id].add(gene_id)
                    mapped = True

                # Check if any ancestor is a slim term
                for ancestor in go_to_ancestors.get(go_id, set()):
                    if ancestor in normalized_slim:
                        slim_to_genes[ancestor].add(gene_id)
                        mapped = True

            if mapped:
                genes_mapped_to_slim.add(gene_id)

        # Step 6: Build result
        total_genes_with_go = len(genes_with_go)

        # Build mapped terms
        mapped_terms = []
        for slim_go_id in sorted(slim_to_genes.keys()):
            term = self.ontology.get_term(slim_go_id)
            if not term:
                continue

            gene_ids = slim_to_genes[slim_go_id]
            gene_count = len(gene_ids)
            frequency = (gene_count / total_genes_with_go * 100) if total_genes_with_go > 0 else 0.0

            genes = [gene_map[g] for g in sorted(gene_ids) if g in gene_map]

            mapped_terms.append(MappedSlimTerm(
                go_id=self.ontology.format_goid(slim_go_id),
                go_term=term.name,
                go_aspect=term.aspect or aspect,
                gene_count=gene_count,
                total_genes=total_genes_with_go,
                frequency_percent=round(frequency, 2),
                genes=genes,
            ))

        # Sort by gene count descending
        mapped_terms.sort(key=lambda x: -x.gene_count)

        # Build "other" genes (have GO but not mapped to slim)
        other_gene_ids = genes_with_go - genes_mapped_to_slim
        other_genes = [gene_map[g] for g in sorted(other_gene_ids) if g in gene_map]

        # Build "not annotated" genes (no GO annotations)
        not_annotated_ids = set(gene_map.keys()) - genes_with_go
        not_annotated_genes = [gene_map[g] for g in sorted(not_annotated_ids) if g in gene_map]

        return SlimMapperResult(
            query_genes_submitted=len(query_genes),
            query_genes_found=len(gene_map),
            query_genes_with_go=total_genes_with_go,
            query_genes_not_found=not_found,
            aspect=aspect,
            mapped_terms=mapped_terms,
            other_genes=other_genes,
            not_annotated_genes=not_annotated_genes,
        )

    def _validate_genes(self, genes: List[str]) -> tuple[Dict[str, MappedGene], List[str]]:
        """Validate gene names and return (gene_map, not_found)."""
        gene_map: Dict[str, MappedGene] = {}
        not_found = []

        for gene in genes:
            gene = gene.strip()
            if not gene:
                continue

            gene_id = self.associations.lookup_gene(gene)
            if gene_id and gene_id not in gene_map:
                gene_map[gene_id] = MappedGene(
                    db_object_id=gene_id,
                    gene_symbol=self.associations.get_gene_symbol(gene_id),
                    gene_name=self.associations.get_gene_name(gene_id),
                )
            elif gene_id is None:
                not_found.append(gene)

        return gene_map, not_found


def parse_slim_terms_from_string(terms_string: str, ontology: GOOntology) -> Set[str]:
    """
    Parse slim terms from a pipe-separated string.

    Format: "term name ; GO:XXXXXXX|term name ; GO:XXXXXXX|..."

    Args:
        terms_string: Pipe-separated slim terms
        ontology: GOOntology for normalization

    Returns:
        Set of normalized GO IDs
    """
    slim_terms = set()

    for term in terms_string.split("|"):
        term = term.strip()
        if not term:
            continue

        # Handle format: "term name ; GO:XXXXXXX" or just "GO:XXXXXXX"
        term = term.replace(";GO:", "; GO:")

        if "; GO:" in term:
            # Extract GO ID after the semicolon
            go_part = term.split("; GO:")[-1].strip()
            if go_part:
                # Handle case where there might be a comma-separated list
                for go_id in go_part.split(","):
                    go_id = go_id.strip()
                    if go_id:
                        normalized = ontology.normalize_go_id("GO:" + go_id.replace("GO:", ""))
                        if normalized:
                            slim_terms.add(normalized)
        elif term.upper().startswith("GO:"):
            normalized = ontology.normalize_go_id(term)
            if normalized:
                slim_terms.add(normalized)

    return slim_terms
