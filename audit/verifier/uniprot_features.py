"""Load and index UniProt features for the deterministic verifier.

The on-disk cache (`audit/output/uniprot_features.json`, written by
`audit/fetch_uniprot_features.py`) maps each
accession to {length, features: [{type, start, end, description, evidence,
ligand?}]}. At verifier load time we build a per-accession position index:

    feat_idx[acc][pos] = [feature, feature, ...]

Range features (Domain, Region, Transmembrane, Topological domain, Signal
peptide, Chain, etc.) expand to occupy every residue in their range. Single-
residue features (Active site, Modified residue, Glycosylation) only land at
their start position. Disulfide bonds occupy both endpoints (each cysteine).

Evidence-code filter:
    Default: keep all annotations except purely-automatic computational
    inference (ECO:0000256, sequence-similarity-automatic). Manual curation,
    experimental, model-based, and human-curator inference all kept.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

# Default cache location: <repo>/audit/output/uniprot_features.json.
CACHE_PATH = (Path(__file__).resolve().parents[1]
              / "output" / "uniprot_features.json")

# Evidence codes we drop by default. Everything else is kept.
SKIP_EVIDENCE = frozenset({
    "ECO:0000256",  # purely automatic sequence-similarity assertion
})

# Range-style features cover [start..end] inclusive. Point-style features
# only land at start (start == end usually). The split is mostly for
# documentation; we expand any feature whose start != end.
RANGE_TYPES = frozenset({
    "Domain", "Region", "Repeat", "Zinc finger", "Transmembrane",
    "Intramembrane", "Topological domain", "Signal peptide", "Propeptide",
    "Transit peptide", "Chain", "Peptide", "Coiled coil", "DNA binding",
    "Motif",
})
POINT_TYPES = frozenset({
    "Active site", "Binding site", "Site", "Modified residue",
    "Glycosylation", "Disulfide bond", "Cross-link", "Lipidation",
})


def _passes_evidence(feature, allow_no_evidence=True):
    ev = feature.get("evidence", [])
    if not ev:
        return allow_no_evidence
    return any(c not in SKIP_EVIDENCE for c in ev)


def load_features(path=None):
    """Load the cache. Returns dict {acc: {length, features: [...]}}."""
    p = Path(path) if path else CACHE_PATH
    if not p.exists():
        raise FileNotFoundError(
            f"UniProt feature cache not found at {p}. "
            "Run fetch_uniprot_features.py first."
        )
    return json.loads(p.read_text())


def build_position_index(entry, evidence_filter=True):
    """Build a {position: [feature_slim]} index for one accession's features.

    For range features, every position in [start..end] gets a reference. For
    Disulfide bond, both endpoints (the two cysteines) get a reference with a
    `partner` field pointing at the other residue.
    """
    idx = defaultdict(list)
    for f in entry.get("features", []):
        if evidence_filter and not _passes_evidence(f):
            continue
        start = f.get("start")
        end = f.get("end")
        if start is None or end is None:
            continue
        ftype = f["type"]
        if ftype == "Disulfide bond" and start != end:
            # Two-cysteine bond — both residues annotated, with partner.
            for self_pos, partner in ((start, end), (end, start)):
                idx[self_pos].append({
                    "type": ftype,
                    "description": f.get("description", ""),
                    "partner": partner,
                    "evidence": f.get("evidence", []),
                })
            continue
        slim = {
            "type": ftype,
            "description": f.get("description", ""),
            "start": start,
            "end": end,
        }
        if "ligand" in f:
            slim["ligand"] = f["ligand"]
        if "evidence" in f:
            slim["evidence"] = f["evidence"]
        if start == end:
            idx[start].append(slim)
        else:
            for pos in range(start, end + 1):
                idx[pos].append(slim)
    return dict(idx)


class FeatureStore:
    """Wrapper around the on-disk cache that builds position indices lazily.

    Usage:
        store = FeatureStore.load_default()
        feats_at_pos = store.at(acc, pos)
        # → list of {type, description, ...} dicts; empty if no annotations
    """

    def __init__(self, raw):
        self._raw = raw
        self._idx_cache: dict[str, dict[int, list]] = {}

    @classmethod
    def load_default(cls):
        return cls(load_features())

    def has(self, acc):
        return acc in self._raw

    def length(self, acc):
        e = self._raw.get(acc)
        return e["length"] if e else None

    def features(self, acc):
        e = self._raw.get(acc)
        return e["features"] if e else []

    def index(self, acc):
        if acc not in self._idx_cache:
            entry = self._raw.get(acc)
            if not entry:
                return {}
            self._idx_cache[acc] = build_position_index(entry)
        return self._idx_cache[acc]

    def at(self, acc, pos):
        """Return the list of feature dicts annotated at `pos` (1-based)."""
        if not pos:
            return []
        return self.index(acc).get(pos, [])

    def has_type_at(self, acc, pos, ftype):
        return any(f["type"] == ftype for f in self.at(acc, pos))


if __name__ == "__main__":
    # Smoke test once the cache exists.
    store = FeatureStore.load_default()
    print(f"Loaded {sum(1 for _ in store._raw):,} accessions")
    # Spot check: BRAF V600 should be a phospho site (Phosphothreonine; by
    # autocatalysis), and lie within the kinase domain.
    feats = store.at("P15056", 600)
    print(f"\nBRAF position 600 features: {len(feats)}")
    for f in feats:
        print(f"  {f['type']}: {f.get('description')}  range=[{f.get('start')},"
              f"{f.get('end')}]")
