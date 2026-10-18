"""UniProt-feature-aware check functions for the deterministic verifier.

These checks use the existing trigger regexes from `triggers.py` to match the
LLM's claim language, but the *truth* source is the per-residue UniProt
feature annotation, not just the AA-identity heuristic.

Each `_check_uniprot_*` function returns a list of {claim_type, claimed,
ref, status} dicts. Wired into `verify_qa()` and `verify_qa_dna()` after the
existing checks.

CONSERVATIVE-PRECISION DESIGN:
We emit `verified` when annotation directly supports a claim, but we only
emit `contradicted` when there's a *direct subtype mismatch at the position*
— i.e., position has an annotation of the relevant type but its description
disagrees with the claim. We do NOT emit `contradicted` based on absence of
annotation, because UniProt feature coverage is uneven and "no annotation
here, but annotation exists elsewhere" produces too many false positives.

This means the feature checks are net coverage-additive (lots of new verified
claims) without polluting the contradiction signal that drives row filtering.

Additional disqualifiers (don't fire at all):
  - "near", "adjacent to", "close to", "around", "next to" within 30 chars
    before the trigger — those describe spatial proximity, not direct claims
  - Hedged / hypothetical / negated context (existing helpers from triggers.py)
  - When ref_matches_genome=False (gated by caller, not here) — gene/transcript
    resolution is suspect so feature lookup is unreliable

Domain claim handling has extra tolerance:
  - Substring/superstring match between claimed name and annotation description
  - A small synonym table for common protein-domain abbreviations
  - Ordinals ("second domain") matched against numeric suffixes ("...2")
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import triggers as T


# Inlined from mutation_verifier_dual to avoid a circular import — these are
# tiny enough to duplicate, and isolating feature_checks lets it be imported
# by both verifiers without creating an import cycle.

_SENT_SPLIT = re.compile(r"(?<=[.!?;])\s+|\n+")


def _sentence_span(answer, idx):
    starts = [0]
    for m in _SENT_SPLIT.finditer(answer):
        starts.append(m.end())
    s = 0
    for st in starts:
        if st <= idx:
            s = st
        else:
            break
    e = len(answer)
    for m in _SENT_SPLIT.finditer(answer, idx):
        e = m.start()
        break
    return s, e


def _is_hypothetical(answer, idx):
    s, e = _sentence_span(answer, idx)
    return bool(T.HYPOTHETICAL_REGEX.search(answer[s:e]))


def _has_negation_before(answer, start, window=30):
    ctx = answer[max(0, start - window):start]
    return bool(T.NEGATION_REGEX.search(ctx))


# ─────────────────────────────────────── claim regexes (feature-specific)

# Most of the claim language for active site / phospho / glyco / disulfide /
# transmem / topology / region is already in `triggers.py` (SPECIAL_TRIGGERS,
# REGION_TRIGGERS). We add a few extras here for kinds the existing tables
# undercover.

# "Domain X" / "in the X domain" — the description word(s) are extracted from
# whatever follows. We capture the word immediately preceding "domain".
DOMAIN_NAMED_RE = re.compile(
    r"\b(?:the\s+)?([A-Za-z][A-Za-z0-9-]{1,30}(?:\s+[A-Za-z][A-Za-z0-9-]{1,30})?)\s+"
    r"domain\b",
    re.IGNORECASE,
)
DOMAIN_GENERIC_RE = re.compile(r"\bdomain\b", re.IGNORECASE)

# Active site / catalytic claims — same as REGION_TRIGGERS["active_site"].
ACTIVE_SITE_RE = T.REGION_TRIGGER_REGEX["active_site"]

# Binding site claims.
BINDING_SITE_RE = T.REGION_TRIGGER_REGEX["binding_site"]

# Phospho / glyco / disulfide language reuse SPECIAL_TRIGGER_COMPILED.
PHOSPHO_RE = T.SPECIAL_TRIGGER_COMPILED["phospho_site"]["regex"]
GLYCO_RE = T.SPECIAL_TRIGGER_COMPILED["glycosylation"]["regex"]
DISULFID_RE = T.SPECIAL_TRIGGER_COMPILED["disulfide"]["regex"]

# Transmembrane and topology reuse REGION_TRIGGERS.
TRANSMEM_RE = T.REGION_TRIGGER_REGEX["transmembrane"]
EXTRACELLULAR_RE = T.REGION_TRIGGER_REGEX["extracellular"]
INTRACELLULAR_RE = T.REGION_TRIGGER_REGEX["intracellular"]

# Signal-peptide claims.
SIGNAL_PEPTIDE_RE = re.compile(
    r"\bsignal\s+(?:peptide|sequence)\b", re.IGNORECASE,
)


# Don't fire feature checks when the trigger is preceded by spatial-proximity
# language — "near the binding site" is not a claim that THIS residue is the
# binding site.
PROXIMITY_DISQUALIFIER_RE = re.compile(
    r"\b(?:near|adjacent\s+to|close\s+to|next\s+to|around|nearby|"
    r"close\s+to\s+the|in\s+the\s+vicinity\s+of|proximal\s+to|"
    r"flanking|surrounding)\s*$",
    re.IGNORECASE,
)

# Domain-name synonyms — the LLM's wording vs UniProt's wording for common
# domain families. Keys are LLM-style fragments; values are the substrings
# we'd expect to see in the UniProt description for an equivalent annotation.
DOMAIN_SYNONYMS = {
    "igsf": ("ig-like", "immunoglobulin"),
    "ig":   ("ig-like", "immunoglobulin"),
    "immunoglobulin": ("ig-like",),
    "brct": ("brct",),
    "kinase": ("kinase",),
    "catalytic": ("kinase", "catalytic"),
    "atp-binding": ("atp",),
    "dna-binding": ("dna",),
    "rna-binding": ("rna",),
    "egf": ("egf",),
    "ank":  ("ankyrin",),
    "ankyrin": ("ankyrin",),
    "sh2":  ("sh2",),
    "sh3":  ("sh3",),
    "ph":   ("ph",),
    "pdz":  ("pdz",),
    "ww":   ("ww",),
    "fha":  ("fha",),
    "rrm":  ("rrm",),
    "kelch": ("kelch",),
    "wd40": ("wd",),
    "wd":   ("wd",),
    "tm":   ("transmembrane",),
}

# Ordinal words → numeric suffixes commonly used in UniProt descriptions
# (e.g. "second TM domain" should match a TRANSMEM with "...2" or be the
# 2nd-occurrence of TRANSMEM in the protein).
ORDINAL_TO_NUM = {
    "first": 1, "1st": 1,
    "second": 2, "2nd": 2,
    "third": 3, "3rd": 3,
    "fourth": 4, "4th": 4,
    "fifth": 5, "5th": 5,
    "sixth": 6, "6th": 6,
    "seventh": 7, "7th": 7,
    "eighth": 8, "8th": 8,
    "ninth": 9, "9th": 9,
    "tenth": 10, "10th": 10,
}


def _has_feature_type(store, acc, ftype):
    """Does this protein have ANY feature of this type?"""
    for f in store.features(acc):
        if f["type"] == ftype:
            return True
    return False


def _has_subtype_at(features_at_pos, ftype, desc_keyword):
    """Is there a feature of `ftype` at this position whose description
    contains `desc_keyword` (case-insensitive)?"""
    kw = desc_keyword.lower()
    for f in features_at_pos:
        if f["type"] != ftype:
            continue
        if kw in (f.get("description") or "").lower():
            return True
    return False


def _disqualified_by_proximity(answer, start, window=30):
    ctx = answer[max(0, start - window):start]
    return bool(PROXIMITY_DISQUALIFIER_RE.search(ctx))


def _domain_name_matches(claimed, annotation_desc):
    """Conservative match between LLM's claimed domain name and a UniProt
    DOMAIN feature's description. Returns True if any of:
      - claimed is a substring of annotation (case-insensitive)
      - annotation is a substring of claimed
      - DOMAIN_SYNONYMS maps claimed to a token contained in annotation
      - claimed is an ordinal that matches the annotation's numeric suffix
    """
    c = (claimed or "").lower().strip()
    a = (annotation_desc or "").lower().strip()
    if not c or not a:
        return False
    if c in a or a in c:
        return True
    for syn_token in DOMAIN_SYNONYMS.get(c, ()):
        if syn_token in a:
            return True
    # Ordinals: "second" → look for trailing digit "2"
    if c in ORDINAL_TO_NUM:
        n = str(ORDINAL_TO_NUM[c])
        if a.rstrip().endswith(n):
            return True
    return False


# ─────────────────────────────────────── per-feature check functions

def _trigger_ok(answer, m):
    """Common precondition: hedge / negation / proximity disqualifier checks."""
    if _is_hypothetical(answer, m.start()):
        return False
    if _has_negation_before(answer, m.start()):
        return False
    if _disqualified_by_proximity(answer, m.start()):
        return False
    return True


def _check_uniprot_domain(answer, facts, store):
    """Verify only when position is in some Domain. Direct subtype check
    using `_domain_name_matches` (with synonym + ordinal handling)."""
    out = []
    acc = facts.get("uniprot_id")
    pos = facts.get("position")
    if not (acc and pos and store.has(acc)):
        return out

    pos_features = store.at(acc, pos)
    domains_here = [f for f in pos_features if f["type"] == "Domain"]
    if not domains_here:
        return out  # conservative: no Domain annotation here → don't fire

    seen_claims = set()
    for m in DOMAIN_NAMED_RE.finditer(answer):
        if not _trigger_ok(answer, m):
            continue
        claimed_name = m.group(1).strip().lower()
        if claimed_name in {"the", "this", "that", "such", "of", "a", "an",
                            "another", "next", "previous", "same", "different",
                            "mutation", "variant", "amino", "residue",
                            "protein", "extra", "intra", "single", "first",
                            "wild-type", "wildtype", "mutant", "type"}:
            continue
        if claimed_name in seen_claims:
            continue
        seen_claims.add(claimed_name)

        matched = any(
            _domain_name_matches(claimed_name, d.get("description"))
            for d in domains_here
        )
        if matched:
            out.append({"claim_type": "uniprot_domain",
                        "claimed": claimed_name,
                        "ref": [d.get("description", "") for d in domains_here],
                        "status": "verified"})
        # else: claim names a different domain than annotation says — but
        # synonym tables are incomplete, so we DON'T mark contradicted.
        # Skip silently to keep contradiction-precision high.
    return out


def _check_uniprot_active_site(answer, facts, store):
    """Verify only — never contradict (annotation may be incomplete)."""
    out = []
    acc = facts.get("uniprot_id")
    pos = facts.get("position")
    if not (acc and pos and store.has(acc)):
        return out
    m = ACTIVE_SITE_RE.search(answer)
    if not m or not _trigger_ok(answer, m):
        return out
    if any(f["type"] == "Active site" for f in store.at(acc, pos)):
        out.append({"claim_type": "uniprot_active_site",
                    "claimed": True, "ref": True, "status": "verified"})
    return out


def _check_uniprot_binding_site(answer, facts, store):
    out = []
    acc = facts.get("uniprot_id")
    pos = facts.get("position")
    if not (acc and pos and store.has(acc)):
        return out
    m = BINDING_SITE_RE.search(answer)
    if not m or not _trigger_ok(answer, m):
        return out
    if any(f["type"] == "Binding site" for f in store.at(acc, pos)):
        out.append({"claim_type": "uniprot_binding_site",
                    "claimed": True, "ref": True, "status": "verified"})
    return out


def _check_uniprot_phospho(answer, facts, store):
    """Phospho claim. Verified when annotation has Phospho-X here. Contradicted
    when annotation has Modified residue here but description is NOT phospho
    (e.g. claimed phospho but annotation says N-acetyl)."""
    out = []
    acc = facts.get("uniprot_id")
    pos = facts.get("position")
    if not (acc and pos and store.has(acc)):
        return out
    m = PHOSPHO_RE.search(answer)
    if not m or not _trigger_ok(answer, m):
        return out

    pos_features = store.at(acc, pos)
    mods_here = [f for f in pos_features if f["type"] == "Modified residue"]
    if not mods_here:
        return out  # no annotation here → not_applicable
    is_phospho = any("phospho" in (f.get("description") or "").lower()
                     for f in mods_here)
    status = "verified" if is_phospho else "contradicted"
    out.append({"claim_type": "uniprot_phospho_site",
                "claimed": "phospho",
                "ref": [f.get("description", "") for f in mods_here],
                "status": status})
    return out


def _check_uniprot_glyco(answer, facts, store):
    out = []
    acc = facts.get("uniprot_id")
    pos = facts.get("position")
    if not (acc and pos and store.has(acc)):
        return out
    m = GLYCO_RE.search(answer)
    if not m or not _trigger_ok(answer, m):
        return out
    if any(f["type"] == "Glycosylation" for f in store.at(acc, pos)):
        out.append({"claim_type": "uniprot_glycosylation",
                    "claimed": True, "ref": True, "status": "verified"})
    return out


def _check_uniprot_disulfide(answer, facts, store):
    """Disulfide claim. Verified when residue is annotated as part of a
    disulfide bond (with optional partner position from the index).
    Contradicted only when the WT residue isn't a Cys (impossible to form
    a disulfide); otherwise silent."""
    out = []
    acc = facts.get("uniprot_id")
    pos = facts.get("position")
    wt = facts.get("wt_aa")
    if not (acc and pos and store.has(acc)):
        return out
    m = DISULFID_RE.search(answer)
    if not m or not _trigger_ok(answer, m):
        return out
    pos_features = store.at(acc, pos)
    has_disulfide_here = any(f["type"] == "Disulfide bond"
                             for f in pos_features)
    if has_disulfide_here:
        out.append({"claim_type": "uniprot_disulfide",
                    "claimed": True, "ref": True, "status": "verified"})
    elif wt and wt != "C":
        # Can't form a disulfide on a non-Cys residue
        out.append({"claim_type": "uniprot_disulfide",
                    "claimed": True, "ref": f"wt={wt}_not_C",
                    "status": "contradicted"})
    return out


def _check_uniprot_transmem(answer, facts, store):
    out = []
    acc = facts.get("uniprot_id")
    pos = facts.get("position")
    if not (acc and pos and store.has(acc)):
        return out
    m = TRANSMEM_RE.search(answer)
    if not m or not _trigger_ok(answer, m):
        return out
    if any(f["type"] == "Transmembrane" for f in store.at(acc, pos)):
        out.append({"claim_type": "uniprot_transmembrane",
                    "claimed": True, "ref": True, "status": "verified"})
    return out


def _check_uniprot_topology(answer, facts, store):
    """LLM claims extracellular / intracellular / cytoplasmic at this residue.
    Verified iff matching topological domain annotation here. Contradicted
    iff annotation here says the OPPOSITE topology."""
    out = []
    acc = facts.get("uniprot_id")
    pos = facts.get("position")
    if not (acc and pos and store.has(acc)):
        return out
    pos_features = store.at(acc, pos)
    topo_here = [f for f in pos_features if f["type"] == "Topological domain"]
    if not topo_here:
        return out
    desc_lower = " ".join((t.get("description") or "").lower() for t in topo_here)

    for claim_label, regex, opp_words in (
        ("extracellular", EXTRACELLULAR_RE,
         ("cytoplasmic", "intracellular", "cytosolic")),
        ("intracellular", INTRACELLULAR_RE,
         ("extracellular", "lumenal")),
    ):
        m = regex.search(answer)
        if not m or not _trigger_ok(answer, m):
            continue
        same = (claim_label in desc_lower
                or (claim_label == "intracellular"
                    and any(w in desc_lower
                            for w in ("cytoplasmic", "intracellular",
                                      "cytosolic"))))
        opp = any(w in desc_lower for w in opp_words)
        if same:
            out.append({"claim_type": "uniprot_topology",
                        "claimed": claim_label,
                        "ref": [t.get("description", "") for t in topo_here],
                        "status": "verified"})
        elif opp:
            out.append({"claim_type": "uniprot_topology",
                        "claimed": claim_label,
                        "ref": [t.get("description", "") for t in topo_here],
                        "status": "contradicted"})
    return out


def _check_uniprot_signal_peptide(answer, facts, store):
    out = []
    acc = facts.get("uniprot_id")
    pos = facts.get("position")
    if not (acc and pos and store.has(acc)):
        return out
    m = SIGNAL_PEPTIDE_RE.search(answer)
    if not m or not _trigger_ok(answer, m):
        return out
    if any(f["type"] == "Signal peptide" for f in store.at(acc, pos)):
        out.append({"claim_type": "uniprot_signal_peptide",
                    "claimed": True, "ref": True, "status": "verified"})
    return out


# ─────────────────────────────────────── orchestration

def run_uniprot_feature_checks(answer, facts, store):
    """Run all 8 UniProt feature checks. Returns combined list of check dicts."""
    if not answer or not facts or not store:
        return []
    checks = []
    checks += _check_uniprot_domain(answer, facts, store)
    checks += _check_uniprot_active_site(answer, facts, store)
    checks += _check_uniprot_binding_site(answer, facts, store)
    checks += _check_uniprot_phospho(answer, facts, store)
    checks += _check_uniprot_glyco(answer, facts, store)
    checks += _check_uniprot_disulfide(answer, facts, store)
    checks += _check_uniprot_transmem(answer, facts, store)
    checks += _check_uniprot_topology(answer, facts, store)
    checks += _check_uniprot_signal_peptide(answer, facts, store)
    return checks


if __name__ == "__main__":
    """Smoke test: run all feature checks against a fabricated answer for a
    known well-annotated mutation."""
    from uniprot_features import FeatureStore
    store = FeatureStore.load_default()
    print(f"Loaded {sum(1 for _ in store._raw):,} accessions")

    # BRAF V600E — kinase domain, autophosphorylation site
    facts = {"uniprot_id": "P15056", "position": 600, "wt_aa": "V",
             "mut_aa": "E"}

    test_answers = [
        "The V600E mutation is located in the kinase domain.",
        "The V600E mutation disrupts a phosphorylation site.",
        "The V600E mutation creates a new active site.",
        "V600E is in the BRCT domain of BRAF.",
        "V600E creates a new disulfide bond.",
        "The mutation is in a transmembrane segment.",
        "The mutation is in a glycosylation site.",
    ]
    for a in test_answers:
        print(f"\n  answer: {a}")
        checks = run_uniprot_feature_checks(a, facts, store)
        for c in checks:
            sym = "✓" if c["status"] == "verified" else "✗" if c["status"] == "contradicted" else "?"
            print(f"    [{sym}] {c['claim_type']}: claimed={c['claimed']!r} "
                  f"ref={c['ref']!r}")
