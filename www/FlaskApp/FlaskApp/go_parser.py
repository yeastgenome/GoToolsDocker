"""
GO Data Parser - Parse OBO ontology and GAF annotation files.

Provides data structures for GO Term Finder and GO Slim Mapper analysis.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Set, List, Tuple, Optional
import os
import pickle


# Namespace to aspect code mapping
NAMESPACE_TO_ASPECT = {
    "biological_process": "P",
    "molecular_function": "F",
    "cellular_component": "C",
}

ASPECT_TO_NAMESPACE = {v: k for k, v in NAMESPACE_TO_ASPECT.items()}

ASPECT_NAMES = {
    "P": "Biological Process",
    "F": "Molecular Function",
    "C": "Cellular Component",
}


@dataclass
class GOTerm:
    """Represents a GO term from the ontology."""
    go_id: str
    name: str = ""
    namespace: str = ""
    aspect: str = ""  # P, F, or C
    is_obsolete: bool = False
    alt_ids: Set[str] = field(default_factory=set)
    parents: Set[Tuple[str, str]] = field(default_factory=set)  # (parent_id, relation)


@dataclass
class GeneAnnotation:
    """Represents a gene annotation from a GAF file."""
    db: str
    db_object_id: str
    db_object_symbol: str
    qualifier: str
    go_id: str
    reference: str
    evidence_code: str
    aspect: str
    db_object_name: str
    synonyms: List[str]


class GOOntology:
    """
    Parsed GO ontology with term hierarchy and ancestor relationships.
    """

    def __init__(self):
        self.terms: Dict[str, GOTerm] = {}
        self.alt_id_to_primary: Dict[str, str] = {}
        self._parents_map: Dict[str, Set[str]] = {}  # go_id -> set of parent go_ids
        self._ancestors_cache: Dict[str, Set[str]] = {}

    def load_obo(self, obo_path: str, cache_path: Optional[str] = None) -> None:
        """
        Load GO ontology from an OBO file.

        Args:
            obo_path: Path to the OBO file
            cache_path: Optional pickle cache path
        """
        # Try loading from cache first
        if cache_path and os.path.exists(cache_path):
            obo_mtime = os.path.getmtime(obo_path)
            cache_mtime = os.path.getmtime(cache_path)
            if cache_mtime > obo_mtime:
                try:
                    with open(cache_path, "rb") as f:
                        cached = pickle.load(f)
                        self.terms = cached["terms"]
                        self.alt_id_to_primary = cached["alt_id_to_primary"]
                        self._parents_map = cached["parents_map"]
                        return
                except Exception:
                    pass  # Fall through to parse

        self._parse_obo(obo_path)
        self._build_parents_map()

        # Save cache
        if cache_path:
            try:
                with open(cache_path, "wb") as f:
                    pickle.dump({
                        "terms": self.terms,
                        "alt_id_to_primary": self.alt_id_to_primary,
                        "parents_map": self._parents_map,
                    }, f, protocol=pickle.HIGHEST_PROTOCOL)
            except Exception:
                pass  # Ignore cache write failures

    def _parse_obo(self, obo_path: str) -> None:
        """Parse an OBO file and populate terms."""
        current_term: Optional[GOTerm] = None
        in_term = False

        with open(obo_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")

                if not line:
                    continue

                if line == "[Term]":
                    if in_term and current_term:
                        self._commit_term(current_term)
                    current_term = None
                    in_term = True
                    continue

                if line.startswith("["):
                    if in_term and current_term:
                        self._commit_term(current_term)
                    in_term = False
                    current_term = None
                    continue

                if not in_term:
                    continue

                if line.startswith("id: "):
                    go_id = line[4:].strip()
                    current_term = GOTerm(go_id=go_id)
                elif current_term is None:
                    continue
                elif line.startswith("name: "):
                    current_term.name = line[6:].strip()
                elif line.startswith("namespace: "):
                    ns = line[11:].strip()
                    current_term.namespace = ns
                    current_term.aspect = NAMESPACE_TO_ASPECT.get(ns, "")
                elif line.startswith("is_obsolete: "):
                    current_term.is_obsolete = line[13:].strip().lower() == "true"
                elif line.startswith("alt_id: "):
                    current_term.alt_ids.add(line[8:].strip())
                elif line.startswith("is_a: "):
                    parent_id = line[6:].split(" ! ")[0].strip()
                    current_term.parents.add((parent_id, "is_a"))
                elif line.startswith("relationship: "):
                    rest = line[14:]
                    parts = rest.split()
                    if len(parts) >= 2 and parts[0] == "part_of":
                        current_term.parents.add((parts[1], "part_of"))

        if in_term and current_term:
            self._commit_term(current_term)

    def _commit_term(self, term: GOTerm) -> None:
        """Add a term to the ontology."""
        if not term.go_id:
            return
        self.terms[term.go_id] = term
        for alt_id in term.alt_ids:
            self.alt_id_to_primary[alt_id] = term.go_id

    def _build_parents_map(self) -> None:
        """Build the parent relationship map."""
        for go_id, term in self.terms.items():
            self._parents_map[go_id] = {
                parent_id for parent_id, _ in term.parents
            }

    def normalize_go_id(self, go_id: str) -> Optional[str]:
        """Normalize a GO ID and map alt_ids to primary."""
        if not go_id:
            return None
        go_id = go_id.strip().upper()
        if go_id.startswith("GO:"):
            go_id = "GO:" + go_id[3:].zfill(7)
        elif go_id.isdigit():
            go_id = "GO:" + go_id.zfill(7)
        else:
            return None

        # Map alt_id to primary
        go_id = self.alt_id_to_primary.get(go_id, go_id)
        return go_id if go_id in self.terms else None

    def get_ancestors(self, go_id: str) -> Set[str]:
        """
        Get all ancestors of a GO term (transitive closure of parents).

        Returns set of ancestor GO IDs (not including the term itself).
        """
        if go_id in self._ancestors_cache:
            return self._ancestors_cache[go_id]

        ancestors = set()
        to_visit = list(self._parents_map.get(go_id, set()))

        while to_visit:
            parent = to_visit.pop()
            if parent not in ancestors and parent in self.terms:
                ancestors.add(parent)
                to_visit.extend(self._parents_map.get(parent, set()))

        self._ancestors_cache[go_id] = ancestors
        return ancestors

    def get_term(self, go_id: str) -> Optional[GOTerm]:
        """Get a GO term by ID."""
        normalized = self.normalize_go_id(go_id)
        if normalized:
            return self.terms.get(normalized)
        return None

    def format_goid(self, go_id: str) -> str:
        """Format a GO ID as GO:XXXXXXX."""
        if go_id.startswith("GO:"):
            return "GO:" + go_id[3:].zfill(7)
        return "GO:" + go_id.zfill(7)


class GeneAssociations:
    """
    Parsed gene association file (GAF format) with gene-to-GO mappings.
    """

    def __init__(self, ontology: GOOntology):
        self.ontology = ontology
        # Gene name -> set of db_object_ids (for deduplication)
        self.gene_to_ids: Dict[str, Set[str]] = defaultdict(set)
        # db_object_id -> gene symbol
        self.id_to_symbol: Dict[str, str] = {}
        # db_object_id -> gene name (description)
        self.id_to_name: Dict[str, str] = {}
        # db_object_id -> set of synonyms
        self.id_to_synonyms: Dict[str, Set[str]] = defaultdict(set)
        # db_object_id -> set of (go_id, evidence_code, aspect)
        self.id_to_annotations: Dict[str, Set[Tuple[str, str, str]]] = defaultdict(set)
        # Name mapping for lookup (case-insensitive)
        self._name_to_id: Dict[str, str] = {}

    def load_gaf(self, gaf_path: str) -> None:
        """
        Load gene associations from a GAF file.

        Args:
            gaf_path: Path to the GAF file
        """
        with open(gaf_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("!") or not line.strip():
                    continue

                cols = line.rstrip("\n").split("\t")
                if len(cols) < 15:
                    continue

                db = cols[0]
                db_object_id = cols[1]
                db_object_symbol = cols[2]
                qualifier = cols[3]
                go_id = cols[4]
                evidence_code = cols[6]
                aspect = cols[8]
                db_object_name = cols[9] if len(cols) > 9 else ""
                synonyms_str = cols[10] if len(cols) > 10 else ""

                # Skip NOT annotations
                if qualifier and "NOT" in qualifier.upper():
                    continue

                # Normalize GO ID
                normalized_go_id = self.ontology.normalize_go_id(go_id)
                if not normalized_go_id:
                    continue

                # Record gene info
                self.id_to_symbol[db_object_id] = db_object_symbol
                self.id_to_name[db_object_id] = db_object_name

                # Parse synonyms
                synonyms = [s.strip() for s in synonyms_str.split("|") if s.strip()]
                self.id_to_synonyms[db_object_id].update(synonyms)

                # Build name mappings (case-insensitive)
                self._name_to_id[db_object_id.upper()] = db_object_id
                self._name_to_id[db_object_symbol.upper()] = db_object_id
                for syn in synonyms:
                    self._name_to_id[syn.upper()] = db_object_id

                # Record annotation
                self.id_to_annotations[db_object_id].add((
                    normalized_go_id,
                    evidence_code,
                    aspect,
                ))

                # Gene name mapping
                self.gene_to_ids[db_object_symbol.upper()].add(db_object_id)
                self.gene_to_ids[db_object_id.upper()].add(db_object_id)
                for syn in synonyms:
                    self.gene_to_ids[syn.upper()].add(db_object_id)

    def lookup_gene(self, gene_name: str) -> Optional[str]:
        """
        Look up a gene by name, symbol, or ID.

        Returns db_object_id if found, None otherwise.
        """
        gene_upper = gene_name.strip().upper()
        # Remove SGD: prefix if present
        if gene_upper.startswith("SGD:"):
            gene_upper = gene_upper[4:]
        return self._name_to_id.get(gene_upper)

    def get_gene_go_annotations(
        self,
        db_object_id: str,
        aspect: Optional[str] = None,
        evidence_codes: Optional[Set[str]] = None,
    ) -> Set[str]:
        """
        Get GO annotations for a gene.

        Args:
            db_object_id: The gene ID
            aspect: Optional aspect filter (P, F, or C)
            evidence_codes: Optional set of evidence codes to include

        Returns:
            Set of GO IDs annotated to this gene
        """
        go_ids = set()
        for go_id, ev_code, ann_aspect in self.id_to_annotations.get(db_object_id, set()):
            if aspect and ann_aspect != aspect:
                continue
            if evidence_codes and ev_code not in evidence_codes:
                continue
            go_ids.add(go_id)
        return go_ids

    def get_gene_go_annotations_with_ancestors(
        self,
        db_object_id: str,
        aspect: Optional[str] = None,
        evidence_codes: Optional[Set[str]] = None,
    ) -> Set[str]:
        """
        Get GO annotations for a gene including all ancestor terms.

        Args:
            db_object_id: The gene ID
            aspect: Optional aspect filter (P, F, or C)
            evidence_codes: Optional set of evidence codes to include

        Returns:
            Set of GO IDs (direct + inherited ancestors)
        """
        direct_go_ids = self.get_gene_go_annotations(
            db_object_id, aspect, evidence_codes
        )

        all_go_ids = set(direct_go_ids)
        for go_id in direct_go_ids:
            ancestors = self.ontology.get_ancestors(go_id)
            # Filter ancestors by aspect if needed
            if aspect:
                ancestors = {
                    a for a in ancestors
                    if self.ontology.terms.get(a, GOTerm("")).aspect == aspect
                }
            all_go_ids.update(ancestors)

        return all_go_ids

    def get_all_annotated_genes(
        self,
        aspect: Optional[str] = None,
        evidence_codes: Optional[Set[str]] = None,
    ) -> Set[str]:
        """
        Get all genes that have GO annotations.

        Args:
            aspect: Optional aspect filter
            evidence_codes: Optional evidence code filter

        Returns:
            Set of db_object_ids
        """
        genes = set()
        for db_object_id, annotations in self.id_to_annotations.items():
            for go_id, ev_code, ann_aspect in annotations:
                if aspect and ann_aspect != aspect:
                    continue
                if evidence_codes and ev_code not in evidence_codes:
                    continue
                genes.add(db_object_id)
                break
        return genes

    def get_gene_symbol(self, db_object_id: str) -> str:
        """Get gene symbol for a db_object_id."""
        return self.id_to_symbol.get(db_object_id, db_object_id)

    def get_gene_name(self, db_object_id: str) -> str:
        """Get gene name (description) for a db_object_id."""
        return self.id_to_name.get(db_object_id, "")


def load_slim_terms(slim_path: str, ontology: GOOntology) -> Set[str]:
    """
    Load GO Slim terms from a slim file (list format or OBO format).

    Args:
        slim_path: Path to the slim file
        ontology: GOOntology instance for normalization

    Returns:
        Set of GO IDs in the slim
    """
    slim_terms = set()

    # Check if it's an OBO file
    with open(slim_path, "r", encoding="utf-8") as f:
        head = f.read(256)
        is_obo = "[Term]" in head or head.startswith("format-version:")

    if is_obo:
        # Parse as OBO
        slim_ontology = GOOntology()
        slim_ontology._parse_obo(slim_path)
        for go_id, term in slim_ontology.terms.items():
            if not term.is_obsolete:
                normalized = ontology.normalize_go_id(go_id)
                if normalized:
                    slim_terms.add(normalized)
    else:
        # Parse as list file (format: "term name ; GO:XXXXXXX")
        import re
        with open(slim_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.split("#")[0].strip()
                if not line:
                    continue

                # Extract GO IDs
                matches = re.findall(r"GO:\s*\d{1,7}|\d{1,7}", line, flags=re.IGNORECASE)
                for match in matches:
                    if match.upper().startswith("GO:"):
                        num = re.sub(r"(?i)GO:\s*", "", match)
                    else:
                        num = match
                    go_id = "GO:" + num.zfill(7)
                    normalized = ontology.normalize_go_id(go_id)
                    if normalized:
                        slim_terms.add(normalized)
                        break  # Take first GO ID per line

    return slim_terms
