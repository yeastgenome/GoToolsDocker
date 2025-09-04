#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python rewrite of the classic Perl `map2slim` script.

Feature parity (practical subset):
- Read a GO slim OBO and a full GO ontology OBO (one or many files, or a cached pickle).
- Map an association file (GAF by default; GFF-ish with --gff) to slim terms:
    * For each input annotation term T, emit one output line for each "most pertinent"
      slim ancestor (nearest slim ancestor by shortest path; if multiple and one
      subsumes another, keep the more specific only).
- Optional count mode (-c/--count) prints, per slim term, the count of distinct
  products mapped directly (leaf) and via any descendant (all). With -t/--tab,
  indent according to the slim DAG's BFS depth (best-effort; DAGs aren't trees).
- Optional mapping dump (--outmap) for *all* terms in the full ontology, listing
  direct slim(s) and all slim ancestors. Optional --shownames adds names.
- Optional aspect filter (-a/--aspect F|P|C) applied to association rows.
- Optional ontology cache (--cache path) using pickle for faster reloads.
- Optional input mapping (--inmap) to preload memoized term→slim mappings.
- GFF mode (--gff) interprets col 3 as a GO term name (or SO:XXXX) and col 9 as product.
- Verbose mode (-v/--verbose) for progress logs.

Not (yet) implemented:
- Bucket terms (-b/--bucket). The original script can synthesize "OTHER-..." bucket
  nodes inside the slim; this implementation accepts the flag but **ignores** it,
  emitting a warning. If you need this, call it out and we can add it.
- Database mode (--dbname/--host) and --err redirection.
- Evidence code filtering (--evcode/--ev) is parsed but not applied.

Usage examples:
  Map an association file to slim:
    python map2slim.py goslim_generic.obo gene_ontology.obo gene_association.sgd \
        --out slimmed_gene_association.sgd --aspect F

  Count mode (tab-indented):
    python map2slim.py goslim_generic.obo gene_ontology.obo gene_association.sgd \
        --count --tab --aspect P --out counts.tsv

  Dump term→slim mappings for all terms:
    python map2slim.py goslim_generic.obo gene_ontology.obo --outmap map.txt --shownames

"""
from __future__ import annotations
import argparse
import sys
import os
import pickle
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Dict, Set, List, Tuple, Optional, Iterable

# ----------------------------
# Minimal, robust OBO parsing
# ----------------------------

@dataclass
class GOTerm:
    id: str
    name: str = ""
    namespace: str = ""
    is_obsolete: bool = False
    alt_ids: Set[str] = field(default_factory=set)
    parents: Set[Tuple[str, str]] = field(default_factory=set)  # (parent_id, relation)
    typedef: bool = False


def parse_obo(paths: List[str], verbose: bool = False) -> Tuple[Dict[str, GOTerm], Dict[str, str]]:
    """Parse one or more OBO files, collecting GO terms and basic relationships.
    Returns (terms, alt2primary).
    """
    terms: Dict[str, GOTerm] = {}
    alt2primary: Dict[str, str] = {}

    def commit(cur: Optional[GOTerm]):
        if cur and cur.id and not cur.typedef:
            terms[cur.id] = cur
            for a in cur.alt_ids:
                alt2primary[a] = cur.id

    for path in paths:
        if verbose:
            print(f"[obo] parsing {path}", file=sys.stderr)
        in_term = False
        cur: Optional[GOTerm] = None
        with open(path, "r", encoding="utf-8") as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                if not line:
                    continue
                if line == "[Term]":
                    if in_term:
                        commit(cur)
                    cur = None
                    in_term = True
                    continue
                if line.startswith("["):
                    if in_term:
                        commit(cur)
                        in_term = False
                    cur = None
                    # [Typedef] or others
                    continue
                if not in_term:
                    continue
                if line.startswith("id: "):
                    _id = line[4:].strip()
                    cur = GOTerm(_id)
                elif cur is None:
                    continue
                elif line.startswith("name: "):
                    cur.name = line[6:].strip()
                elif line.startswith("namespace: "):
                    cur.namespace = line[11:].strip()
                elif line.startswith("is_obsolete: "):
                    cur.is_obsolete = (line[13:].strip().lower() == "true")
                elif line.startswith("alt_id: "):
                    cur.alt_ids.add(line[8:].strip())
                elif line.startswith("is_a: "):
                    parent_id = line[6:].split(" ! ")[0].strip()
                    cur.parents.add((parent_id, "is_a"))
                elif line.startswith("relationship: "):
                    # relationship: part_of GO:0008150 ! biological_process
                    rest = line[len("relationship: "):]
                    parts = rest.split()
                    if len(parts) >= 2 and parts[0] == "part_of":
                        cur.parents.add((parts[1], "part_of"))
        if in_term:
            commit(cur)
    return terms, alt2primary


# ----------------------------
# Slim graph utilities
# ----------------------------

ASPECT_NS = {"F": "molecular_function", "P": "biological_process", "C": "cellular_component"}


def load_slim_any(path: str, alt2primary: Dict[str, str], terms: Dict[str, GOTerm], verbose: bool = False) -> Set[str]:
    """Load a slim set from either an OBO *or* a plain list of GO IDs.
    - If file appears to be OBO, parse and take all term IDs from it.
    - Otherwise, read lines and extract GO IDs, mapping alt IDs to primary.
    Obsolete or unknown terms are skipped.
    """
    slim: Set[str] = set()
    looks_obo = False
    try:
        with open(path, "r", encoding="utf-8") as fh:
            head = fh.read(256)
            looks_obo = "[Term]" in head or head.startswith("format-version:")
    except Exception:
        pass

    if looks_obo:
        if verbose:
            print(f"[slim] Detected OBO: {path}", file=sys.stderr)
        slim_terms, _ = parse_obo([path], verbose=False)
        for tid, t in slim_terms.items():
            if t and not t.is_obsolete and tid in terms:
                slim.add(tid)
        return slim

    # treat as list
    if verbose:
        print(f"[slim] Detected list file: {path}", file=sys.stderr)
    import re
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            m = re.search(r"(GO:\d{7}|\d{1,7})", line, flags=re.IGNORECASE)
            if not m:
                continue
            gid = m.group(1)
            if gid.isdigit():
                gid = "GO:" + gid.zfill(7)
            else:
                gid = "GO:" + gid.split("GO:")[1].zfill(7)
            gid = gid.upper()
            if gid in alt2primary:
                gid = alt2primary[gid]
            t = terms.get(gid)
            if not t or t.is_obsolete:
                continue
            slim.add(gid)
    return slim



def normalize_go_id(s: str) -> Optional[str]:
    s = (s or "").strip()
    if not s:
        return None
    import re
    m = re.search(r"(GO:\d{7})", s, flags=re.IGNORECASE)
    if m:
        go = m.group(1).upper()
        return "GO:" + go.split("GO:")[1].zfill(7)
    if s.isdigit() and len(s) <= 7:
        return "GO:" + s.zfill(7)
    return None


def build_parent_map(terms: Dict[str, GOTerm], allowed: Set[str]) -> Dict[str, Set[str]]:
    parents: Dict[str, Set[str]] = defaultdict(set)
    for tid, t in terms.items():
        for pid, rel in t.parents:
            if rel in allowed:
                parents[tid].add(pid)
    return parents


def build_child_map(parents: Dict[str, Set[str]]) -> Dict[str, Set[str]]:
    children: Dict[str, Set[str]] = defaultdict(set)
    for c, ps in parents.items():
        for p in ps:
            children[p].add(c)
    return children


def nearest_slim_ancestors(
    start: str,
    slim: Set[str],
    parents: Dict[str, Set[str]],
    namespace_of: Dict[str, str],
    aspect_ns: str,
) -> Tuple[Set[str], Set[str]]:
    """Return (direct_slim_set, all_slim_ancestors_set) for `start`.
    - direct_slim_set: nearest slim ancestors at minimal depth (if `start` in slim, {start}).
    - all_slim_ancestors_set: all slim ancestors reachable (any depth).
    """
    if start in slim:
        return {start}, {start}

    visited = {start}
    q = deque([(start, 0)])
    best_depth: Optional[int] = None
    direct: Set[str] = set()
    all_anc: Set[str] = set()

    while q:
        node, d = q.popleft()
        for p in parents.get(node, ()):  # traverse upward
            if p in visited:
                continue
            visited.add(p)
            if p in slim:
                all_anc.add(p)
                if best_depth is None or d + 1 < best_depth:
                    best_depth = d + 1
                    direct = {p}
                elif d + 1 == best_depth:
                    direct.add(p)
            else:
                # keep traversing regardless of namespace, but record only matches in aspect
                q.append((p, d + 1))
                if namespace_of.get(p) in (aspect_ns, None):
                    # If parent is in aspect, but not slim, it's still an ancestor that might lead to slim
                    pass
        if best_depth is not None and q and q[0][1] > best_depth:
            break

    # Remove redundant direct ancestors: if A and B are both in `direct` and A is ancestor of B,
    # keep the more specific (B). We'll test ancestor relation via upward closure from B.
    if len(direct) > 1:
        # compute upward closures once
        up: Dict[str, Set[str]] = {}
        for d0 in direct:
            seen = {d0}
            dq = deque([d0])
            anc = set()
            while dq:
                x = dq.popleft()
                for p in parents.get(x, ()):  # climbing further up
                    if p not in seen:
                        seen.add(p)
                        anc.add(p)
                        dq.append(p)
            up[d0] = anc
        keep = set(direct)
        for a in list(direct):
            for b in list(direct):
                if a == b:
                    continue
                if a in up.get(b, set()):  # a is ancestor of b => drop a
                    if a in keep:
                        keep.remove(a)
        direct = keep

    # Compute all slim ancestors (any depth)
    # We already captured during BFS when encountering slim nodes at any depth
    # but if no direct found, we still want all reachable slim ancestors.
    if not all_anc:
        # Explore fully upward to collect any slim ancestors
        visited = {start}
        dq = deque([start])
        while dq:
            x = dq.popleft()
            for p in parents.get(x, ("")):
                if p in visited:
                    continue
                visited.add(p)
                if p in slim:
                    all_anc.add(p)
                dq.append(p)

    return direct, all_anc


# ----------------------------
# Association file I/O
# ----------------------------

def is_comment(line: str) -> bool:
    return line.startswith("!") or line.strip() == ""


def parse_assoc_line(line: str) -> Optional[List[str]]:
    if is_comment(line):
        return None
    cols = line.rstrip("\n").split("\t")
    return cols if cols else None


# ----------------------------
# Core application
# ----------------------------

def load_inmap(path: str) -> Dict[str, Tuple[List[str], List[str]]]:
    memo: Dict[str, Tuple[List[str], List[str]]] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.lstrip().startswith("#"):
                continue
            # Format: "ACC => a b c // x y z"
            if "=>" in line and "//" in line:
                left, right = line.split("=>", 1)
                direct_s, all_s = right.split("//", 1)
                acc = left.strip().split()[0]
                direct = [t for t in direct_s.strip().split() if t]
                all_ = [t for t in all_s.strip().split() if t]
                memo[acc] = (direct, all_)
    return memo


def dump_outmap(
    out_path: str,
    terms: Dict[str, GOTerm],
    slim: Set[str],
    parents: Dict[str, Set[str]],
    shownames: bool,
    verbose: bool,
):
    namespace_of = {tid: t.namespace for tid, t in terms.items()}
    aspect_ns_set = set(ASPECT_NS.values())
    with open(out_path, "w", encoding="utf-8") as out:
        all_terms = sorted(terms.keys())
        for acc in all_terms:
            # restrict to proper GO namespaces
            if namespace_of.get(acc) not in aspect_ns_set:
                continue
            direct, all_anc = nearest_slim_ancestors(acc, slim, parents, namespace_of, namespace_of.get(acc, ""))
            if shownames:
                def fmt(ids: Iterable[str]) -> str:
                    return " ".join(f"{i} \"{terms.get(i).name if terms.get(i) else '?'}\"" for i in ids)
                line = f"{acc} => {fmt(sorted(direct))} // {fmt(sorted(all_anc))}\n"
            else:
                line = f"{acc} => {' '.join(sorted(direct))} // {' '.join(sorted(all_anc))}\n"
            out.write(line)
    if verbose:
        print(f"[outmap] wrote {out_path}", file=sys.stderr)


def process_assocs(
    assoc_path: str,
    out_path: Optional[str],
    terms: Dict[str, GOTerm],
    slim: Set[str],
    parents: Dict[str, Set[str]],
    aspect: Optional[str],
    gff: bool,
    count_mode: bool,
    tab: bool,
    verbose: bool,
):
    namespace_of = {tid: t.namespace for tid, t in terms.items()}

    # Count maps: slim_acc -> set(products)
    leafh: Dict[str, Set[str]] = defaultdict(set)
    allh: Dict[str, Set[str]] = defaultdict(set)

    seen_rows: Set[Tuple[str, ...]] = set()
    out_fh = open(out_path, "w", encoding="utf-8") if out_path else sys.stdout

    # read
    if assoc_path and assoc_path.endswith((".gz", ".Z")):
        import gzip
        opener = gzip.open  # type: ignore
        mode = "rt"
    else:
        opener = open  # type: ignore
        mode = "r"

    total_in = 0
    with opener(assoc_path, mode, encoding="utf-8") as fh:  # type: ignore[arg-type]
        for raw in fh:
            if is_comment(raw):
                continue
            total_in += 1
            cols = parse_assoc_line(raw)
            if not cols:
                continue

            # Extract needed fields
            if gff:
                # GFF-like: col2=source, col3=type (can be name or SO:XXXX), col9=attributes ~ product id
                type_val = cols[2] if len(cols) > 2 else ""
                prod = cols[8] if len(cols) > 8 else ""
                is_not = ""  # GFF path ignores NOT
                # Try as ID or name
                acc = None
                if type_val.startswith("SO:"):
                    acc = type_val
                else:
                    # look up by name
                    # (linear scan; could be indexed by name -> id if needed)
                    for tid, t in terms.items():
                        if t.name == type_val:
                            acc = tid
                            break
                    if acc is None:
                        if verbose:
                            print(f"[warn] GFF type not found in ontology: {type_val}", file=sys.stderr)
                        continue
            else:
                # GAF: 5th col is GO_ID, 2nd col is db_object_id, 9th col is aspect, 4th col may be 'NOT'
                is_not = cols[3] if len(cols) > 3 else ""
                acc = cols[4] if len(cols) > 4 else ""
                prod = cols[1] if len(cols) > 1 else ""

            # skip NOT
            if is_not and str(is_not).strip().lower() == "not":
                continue

            # aspect filter
            if aspect and not gff:
                gaf_aspect = cols[8].strip().upper() if len(cols) > 8 else ""
                if gaf_aspect != aspect:
                    continue

            acc_norm = normalize_go_id(acc)
            if not acc_norm:
                if verbose:
                    print(f"[warn] cannot normalize term: {acc}", file=sys.stderr)
                continue
            term = terms.get(acc_norm)
            if not term or term.is_obsolete:
                continue

            # Find slim ancestors
            aspect_ns = ASPECT_NS.get(aspect, namespace_of.get(acc_norm, "")) if aspect else namespace_of.get(acc_norm, "")
            direct, all_anc = nearest_slim_ancestors(acc_norm, slim, parents, namespace_of, aspect_ns)
            if not direct and not all_anc:
                # No path to slim; drop
                continue

            # Update counts
            for s in direct:
                leafh[s].add(prod)
            for s in all_anc:
                allh[s].add(prod)

            if not count_mode:
                # emit mapped rows (one per direct slim)
                if gff:
                    for s in sorted(direct):
                        new_cols = list(cols)
                        new_cols[2] = s
                        key = tuple(new_cols)
                        if key in seen_rows:
                            continue
                        seen_rows.add(key)
                        print("\t".join(new_cols), file=out_fh)
                else:
                    for s in sorted(direct):
                        new_cols = list(cols)
                        new_cols[4] = s
                        key = tuple(new_cols)
                        if key in seen_rows:
                            continue
                        seen_rows.add(key)
                        print("\t".join(new_cols), file=out_fh)

    if not count_mode:
        if out_fh is not sys.stdout:
            out_fh.close()
        if verbose:
            print(f"[done] processed {total_in} assoc lines", file=sys.stderr)
        return

    # COUNT MODE OUTPUT
    # We need to iterate slim DAG and print: acc name (full name) (slim name?) count_leaf count_all obsolete? type?
    # We'll approximate the Perl's output, using BFS depth for indentation when --tab.
    children = build_child_map(parents)

    # Determine slim nodes' depths for indentation (best-effort on DAG):
    # roots = slim terms with no slim parents
    slim_parents = {s: {p for p in parents.get(s, set()) if p in slim} for s in slim}
    roots = [s for s in slim if len(slim_parents.get(s, set())) == 0]

    depth: Dict[str, int] = {r: 0 for r in roots}
    dq = deque(roots)
    while dq:
        x = dq.popleft()
        for c in children.get(x, set()):
            if c in slim:
                nd = depth[x] + 1
                if c not in depth or nd < depth[c]:
                    depth[c] = nd
                    dq.append(c)

    # stable order: by depth, then id
    def sort_key(s: str) -> Tuple[int, str]:
        return (depth.get(s, 0), s)

    out = out_fh
    for s in sorted(slim, key=sort_key):
        t = terms.get(s)
        if not t:
            continue
        count_leaf = len(leafh.get(s, set()))
        count_all = len(allh.get(s, set()))
        indent = " " * (depth.get(s, 0) + 1) if tab else ""
        full_name = t.name
        # In Perl, the 3rd column shows the name in the full ontology (t2->name), often same.
        full_onto_name = full_name
        is_obs = "OBSOLETE" if t.is_obsolete else ""
        onto = t.namespace or ""
        print(f"{indent}{s} {full_name} ({full_onto_name})\t{count_leaf}\t{count_all}\t{is_obs}\t{onto}", file=out)
    if out is not sys.stdout:
        out.close()


# ----------------------------
# CLI
# ----------------------------

def main():
    p = argparse.ArgumentParser(
        description="Map a gene association file to a GO slim (Python rewrite of Perl map2slim)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Positional (like Perl): slimfile, ontfiles..., [assocfile]
    p.add_argument("slimfile", help="GO slim OBO file")
    p.add_argument("ontfiles", nargs="*", help="Full GO ontology OBO file(s), last arg may be assoc file if not using --outmap")

    # Options
    p.add_argument("--dbname", "-d")
    p.add_argument("--host", "-H")
    p.add_argument("--out", "-o", help="Output file (mapped assocs or counts). Defaults to stdout", default=None)
    p.add_argument("--err", "-e")
    p.add_argument("--force", action="store_true")
    p.add_argument("--ontdir", help="Directory containing ontology files (*.obo, *ontology)")
    p.add_argument("--outmap", help="Write term→slim mappings for all terms and exit")
    p.add_argument("--shownames", action="store_true", help="Include names in --outmap output")
    p.add_argument("--cache", help="Pickle cache path for ontology graph (read/write)")
    p.add_argument("--inmap", help="Preload memoized mappings from a previous --outmap file")
    p.add_argument("--gff", action="store_true", help="Treat assoc as GFF-like (col3=type, col9=product)")
    p.add_argument("--aspect", "-a", choices=["F", "P", "C"], help="Restrict to Aspect (GAF col 9)")
    p.add_argument("--count", "-c", action="store_true", help="Count mode instead of mapping rows")
    p.add_argument("--tab", "-t", action="store_true", help="Indent count-mode output by slim depth")
    p.add_argument("--bucket", "-b", help="(Not implemented) Bucket slim file output")
    p.add_argument("--evcode", "--ev", action="append", help="Evidence code filter (parsed, not applied)")
    p.add_argument("--verbose", "-v", action="store_true")

    args = p.parse_args()

    if args.bucket:
        print("[warn] --bucket is accepted but not implemented in this Python rewrite.", file=sys.stderr)

    # Discover files
    ontfiles: List[str] = list(args.ontfiles)
    assocfile: Optional[str] = None
    if args.ontdir:
        # override ontfiles by directory glob like Perl
        import glob
        ontfiles = sorted(glob.glob(os.path.join(args.ontdir, "*.obo")))
        if not ontfiles:
            ontfiles = sorted(glob.glob(os.path.join(args.ontdir, "*ontology")))
    else:
        # If --outmap is set, there is no assoc file; otherwise, last arg is assoc
        if not args.outmap and ontfiles:
            assocfile = ontfiles[-1]
            ontfiles = ontfiles[:-1]

    if not ontfiles:
        p.error("No ontology files provided. Supply one or more OBO files after slimfile, or use --ontdir.")

    # Load ontology, possibly from cache
    if args.cache and os.path.exists(args.cache):
        if args.verbose:
            print(f"[cache] loading ontology from {args.cache}", file=sys.stderr)
        with open(args.cache, "rb") as fh:
            terms, alt2primary = pickle.load(fh)
    else:
        terms, alt2primary = parse_obo([args.slimfile] + ontfiles, verbose=args.verbose)
        if args.cache:
            if args.verbose:
                print(f"[cache] writing ontology to {args.cache}", file=sys.stderr)
            with open(args.cache, "wb") as fh:
                pickle.dump((terms, alt2primary), fh, protocol=pickle.HIGHEST_PROTOCOL)

    # Identify slim set from either an OBO or a plain list of GO IDs
    slim: Set[str] = load_slim_any(args.slimfile, alt2primary, terms, verbose=args.verbose)
    if not slim:
        print(f"[warn] Slim set is empty after loading {args.slimfile}. Check file type and contents.", file=sys.stderr)

    # Build relation graph (is_a + part_of like Perl)
    allowed = {"is_a", "part_of"}
    parents = build_parent_map(terms, allowed)

    # Preload inmap memo if provided (applies only to --outmap shortcutting). We still recompute here for correctness.
    if args.inmap and args.verbose:
        _ = load_inmap(args.inmap)
        print(f"[inmap] loaded memo from {args.inmap} (not strictly required)", file=sys.stderr)

    # outmap mode
    if args.outmap:
        dump_outmap(args.outmap, terms, slim, parents, args.shownames, args.verbose)
        return

    if not assocfile:
        p.error("Association file is required unless --outmap is used.")

    # Map / Count
    process_assocs(
        assocfile, args.out, terms, slim, parents, args.aspect, args.gff, args.count, args.tab, args.verbose
    )


if __name__ == "__main__":
    main()
