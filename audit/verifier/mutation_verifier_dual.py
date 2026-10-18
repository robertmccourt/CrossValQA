"""Protein-side deterministic verifier for MutQA (both answers per row).

For every Q&A row, derive a fact dict from `(ref_sequence, mutated_sequence,
HGVS)` plus UniProt feature annotation, then regex-extract claims from
`qagen_answer_final` (phase-1 answer) and `llm2_answer_final` (phase-2 answer)
and emit per-claim {verified, contradicted} verdicts.

Input:  gold rows JSONL from `audit/build_gold_rows_v2.py`.
Output: <out_dir>/mutation_verification.jsonl (one row per input row; carries
        `row_id`, `qagen_label`, `llm2_label`) and
        <out_dir>/mutation_verification_summary.md (by bucket, split, claim type).

Example (from the repo root):
    python audit/verifier/mutation_verifier_dual.py \
        --input audit/output/_gold_rows_v2.jsonl \
        --out_dir audit/output/protein_v3 \
        --features audit/output/uniprot_features.json
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import triggers as T
from aa_properties import (
    AA_3TO1,
    AA_FULL_TO_1,
    aa_change,
    aa_properties,
)
from topic_bucket import bucket_topic
from uniprot_features import CACHE_PATH as FEATURES_CACHE, FeatureStore, load_features
from feature_checks import run_uniprot_feature_checks


# ----------------------------------------------------------- HGVS parser ---

_HGVS_PROT_SHORT = re.compile(
    r"(?<![A-Za-z0-9])p\.([A-Z])(\d+)([A-Z*])(?![A-Za-z0-9])"
)
_HGVS_PROT_LONG = re.compile(
    r"(?<![A-Za-z0-9])p\.(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|"
    r"Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter)(\d+)"
    r"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|"
    r"Trp|Tyr|Val|Ter|\*|fs|del|dup|ins)",
    re.IGNORECASE,
)
_HGVS_DEL = re.compile(r"p\.(?:[A-Z][a-z]{2}|[A-Z])(\d+)(?:_(?:[A-Z][a-z]{2}|[A-Z])\d+)?del", re.IGNORECASE)
_HGVS_INS = re.compile(r"p\.(?:[A-Z][a-z]{2}|[A-Z])(\d+)_[A-Z]\w*\d+ins", re.IGNORECASE)
_HGVS_DUP = re.compile(r"p\.(?:[A-Z][a-z]{2}|[A-Z])(\d+)(?:_[A-Z]\w*\d+)?dup", re.IGNORECASE)
_HGVS_FS = re.compile(r"p\.([A-Z][a-z]{2}|[A-Z])(\d+)([A-Z][a-z]{2}|[A-Z])?fs", re.IGNORECASE)
_HGVS_NONSENSE = re.compile(r"p\.([A-Z][a-z]{2}|[A-Z])(\d+)(?:Ter|\*)", re.IGNORECASE)


def _norm_aa(token: str) -> str:
    """1-letter uppercase. Accepts 1-letter, 3-letter, or 'Ter'/'*'."""
    if not token:
        return ""
    t = token.strip()
    if t == "*" or t.lower() == "ter":
        return "*"
    if len(t) == 1:
        return t.upper()
    return AA_3TO1.get(t.title(), "")


def parse_hgvs(s: str) -> dict | None:
    """Parse an HGVS-ish mutation descriptor into {wt, pos, mut, kind}.
    Returns None on failure. `kind` ∈ {substitution, nonsense, frameshift,
    deletion, insertion, duplication, unknown}.
    """
    if not s:
        return None
    s = s.strip()

    # nonsense first (catches p.Trp45*, p.W45Ter, p.Trp45Ter)
    m = _HGVS_NONSENSE.search(s)
    if m:
        return {"wt": _norm_aa(m.group(1)), "pos": int(m.group(2)),
                "mut": "*", "kind": "nonsense"}

    # frameshift
    m = _HGVS_FS.search(s)
    if m:
        return {"wt": _norm_aa(m.group(1)), "pos": int(m.group(2)),
                "mut": _norm_aa(m.group(3)) if m.group(3) else "",
                "kind": "frameshift"}

    # duplication
    m = _HGVS_DUP.search(s)
    if m:
        return {"wt": "", "pos": int(m.group(1)), "mut": "", "kind": "duplication"}

    # insertion
    m = _HGVS_INS.search(s)
    if m:
        return {"wt": "", "pos": int(m.group(1)), "mut": "", "kind": "insertion"}

    # deletion
    m = _HGVS_DEL.search(s)
    if m:
        return {"wt": "", "pos": int(m.group(1)), "mut": "", "kind": "deletion"}

    # short substitution: p.R253G
    m = _HGVS_PROT_SHORT.search(s)
    if m:
        wt = _norm_aa(m.group(1))
        mut = _norm_aa(m.group(3))
        if wt and mut and wt != "*":
            kind = "nonsense" if mut == "*" else "substitution"
            return {"wt": wt, "pos": int(m.group(2)), "mut": mut, "kind": kind}

    # long substitution: p.Arg253Gly
    m = _HGVS_PROT_LONG.search(s)
    if m:
        wt = _norm_aa(m.group(1))
        mut_tok = m.group(3)
        if mut_tok.lower() in {"fs", "del", "dup", "ins"}:
            return {"wt": wt, "pos": int(m.group(2)), "mut": "",
                    "kind": mut_tok.lower() if mut_tok.lower() != "fs" else "frameshift"}
        mut = _norm_aa(mut_tok)
        if wt and mut:
            kind = "nonsense" if mut == "*" else "substitution"
            return {"wt": wt, "pos": int(m.group(2)), "mut": mut, "kind": kind}

    return None


# ---------------------------------------------------------- Facts builder ---


def mutation_facts(ref_seq: str, mut_seq: str, hgvs: str,
                   mutation_type: str | None = None,
                   uniprot_id: str | None = None) -> dict | None:
    parsed = parse_hgvs(hgvs)
    if parsed is None:
        return None

    facts: dict[str, Any] = {
        "wt_aa": parsed["wt"],
        "mut_aa": parsed["mut"],
        "position": parsed["pos"],
        "kind": parsed["kind"],
        "mutation_type_label": (mutation_type or "").strip().lower(),
        "protein_length": len(ref_seq) if ref_seq else None,
        "ref_aa_at_pos": None,
        "sanity_ok": False,
        "uniprot_id": uniprot_id,
    }

    # Sequence sanity for substitutions / nonsense / frameshift only.
    if ref_seq and parsed["pos"] and parsed["pos"] <= len(ref_seq):
        facts["ref_aa_at_pos"] = ref_seq[parsed["pos"] - 1]
        if parsed["wt"] and parsed["wt"] != "*":
            facts["sanity_ok"] = (facts["ref_aa_at_pos"] == parsed["wt"])
        else:
            facts["sanity_ok"] = True

    # Position class.
    if facts["protein_length"]:
        L = facts["protein_length"]
        p = parsed["pos"]
        third = max(1, L // 3)
        if p <= third:
            facts["region_class"] = "n_terminal"
        elif p > L - third:
            facts["region_class"] = "c_terminal"
        else:
            facts["region_class"] = "middle"
        facts["dist_from_n_term"] = p
        facts["dist_from_c_term"] = L - p + 1
    else:
        facts["region_class"] = None

    # AA-change derived facts (only for true substitutions).
    if parsed["kind"] == "substitution" and parsed["wt"] and parsed["mut"]:
        facts["change"] = aa_change(parsed["wt"], parsed["mut"])
    else:
        facts["change"] = {"applicable": False}

    # Glycosylation sequon check (N-X-S/T, X != P).
    if (parsed["wt"] == "N" and ref_seq and parsed["pos"] + 1 < len(ref_seq)):
        x = ref_seq[parsed["pos"]]
        s_or_t = ref_seq[parsed["pos"] + 1] if parsed["pos"] + 1 < len(ref_seq) else ""
        facts["wt_in_glyco_sequon"] = (x != "P" and s_or_t in "ST")
    else:
        facts["wt_in_glyco_sequon"] = False

    return facts


# -------------------------------------------------- Sentence / context ---

_SENT_SPLIT = re.compile(r"(?<=[.!?;])\s+|\n+")


def _sentence_span(answer: str, idx: int) -> tuple[int, int]:
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


def _is_hypothetical(answer: str, idx: int) -> bool:
    s, e = _sentence_span(answer, idx)
    return bool(T.HYPOTHETICAL_REGEX.search(answer[s:e]))


def _has_negation_before(answer: str, start: int, window: int = 30) -> bool:
    ctx = answer[max(0, start - window): start]
    return bool(T.NEGATION_REGEX.search(ctx))


def _has_pre_disqualifier(answer: str, start: int, family: str, window: int = 30) -> bool:
    pat = T.DISQUALIFIER_PRE.get(family)
    if pat is None:
        return False
    ctx = answer[max(0, start - window): start]
    return bool(pat.search(ctx))


def _has_post_disqualifier(answer: str, end: int, family: str, window: int = 30) -> bool:
    pat = T.DISQUALIFIER_POST.get(family)
    if pat is None:
        return False
    ctx = answer[end: end + window]
    return bool(pat.search(ctx))


def _dedup_checks(checks: list[dict]) -> list[dict]:
    seen = set()
    out = []
    for c in checks:
        key = (c["claim_type"], str(c.get("claimed")))
        if key in seen:
            continue
        seen.add(key)
        out.append(c)
    return out


# --------------------------------------------------- Per-family verifiers ---


def _check_muttype(answer: str, facts: dict) -> list[dict]:
    """Match mutation-kind triggers and emit a single verdict.

    Prior bug: iterated kinds in dict order and emitted a verdict at the
    first match, so an answer that legitimately mentions BOTH the correct
    kind and an adjacent one (e.g. a frameshift answer that also mentions
    "premature stop codon") would be flagged contradicted whenever the
    "wrong" kind happened to come first in the dict.

    Fix: prefer the reference kind. If any non-hedged, non-negated match
    of the reference kind appears in the answer, emit `verified`. Only
    when no reference-kind match exists do we fall back to the first
    other kind that fires, marking that contradiction.
    """
    out = []
    ref_kind = facts["kind"]

    def _fires(kind: str):
        pat = T.MUTTYPE_TRIGGER_REGEX.get(kind)
        if pat is None:
            return None
        regex = pat["regex"] if isinstance(pat, dict) else pat
        for m in regex.finditer(answer):
            if _is_hypothetical(answer, m.start()):
                continue
            if _has_negation_before(answer, m.start()):
                continue
            return m
        return None

    # Step 1: did the answer name the reference kind anywhere?
    if ref_kind and _fires(ref_kind):
        out.append({"claim_type": "mutation_kind", "claimed": ref_kind,
                    "ref": ref_kind, "status": "verified"})
        return out

    # Step 2: otherwise, fall back to the first other kind that fires.
    for kind in T.MUTTYPE_TRIGGER_REGEX:
        if kind == ref_kind:
            continue
        if _fires(kind):
            out.append({"claim_type": "mutation_kind", "claimed": kind,
                        "ref": ref_kind, "status": "contradicted"})
            break
    return out


def _check_hgvs_inline(answer: str, facts: dict) -> list[dict]:
    """Look for inline 'R253G' / 'p.Arg253Gly' style claims and check identity."""
    out = []
    ref_wt, ref_mut, ref_pos = facts["wt_aa"], facts["mut_aa"], facts["position"]
    if not ref_pos:
        return out

    for m in T.HGVS_SHORT_PATTERN.finditer(answer):
        wt_c, pos_c, mut_c = m.group(1), int(m.group(2)), m.group(3)
        if pos_c != ref_pos:
            continue
        if _is_hypothetical(answer, m.start()):
            continue
        if ref_wt and wt_c != ref_wt:
            out.append({"claim_type": "wt_residue_at_pos",
                        "claimed": wt_c, "ref": ref_wt, "status": "contradicted"})
        elif ref_wt:
            out.append({"claim_type": "wt_residue_at_pos",
                        "claimed": wt_c, "ref": ref_wt, "status": "verified"})
        if ref_mut and ref_mut != "*" and mut_c != ref_mut:
            out.append({"claim_type": "mut_residue_at_pos",
                        "claimed": mut_c, "ref": ref_mut, "status": "contradicted"})
        elif ref_mut and ref_mut != "*":
            out.append({"claim_type": "mut_residue_at_pos",
                        "claimed": mut_c, "ref": ref_mut, "status": "verified"})

    for m in T.HGVS_LONG_PATTERN.finditer(answer):
        wt_c = _norm_aa(m.group(1))
        pos_c = int(m.group(2))
        mut_c = _norm_aa(m.group(3))
        if pos_c != ref_pos:
            continue
        if _is_hypothetical(answer, m.start()):
            continue
        if ref_wt and wt_c and wt_c != ref_wt:
            out.append({"claim_type": "wt_residue_at_pos",
                        "claimed": wt_c, "ref": ref_wt, "status": "contradicted"})
        if ref_mut and mut_c and ref_mut != "*" and mut_c != ref_mut:
            out.append({"claim_type": "mut_residue_at_pos",
                        "claimed": mut_c, "ref": ref_mut, "status": "contradicted"})
    return out


def _check_position(answer: str, facts: dict) -> list[dict]:
    out = []
    ref = facts["position"]
    if not ref:
        return out
    for m in T.POSITION_NUM_PATTERN.finditer(answer):
        try:
            n = int(m.group(1))
        except ValueError:
            continue
        if _is_hypothetical(answer, m.start()):
            continue
        status = "verified" if n == ref else "contradicted"
        out.append({"claim_type": "position", "claimed": n, "ref": ref, "status": status})
    return out


def _check_protein_length(answer: str, facts: dict) -> list[dict]:
    out = []
    ref = facts.get("protein_length")
    if not ref:
        return out
    for m in T.PROTEIN_LENGTH_PATTERN.finditer(answer):
        n = int(m.group(1))
        if _is_hypothetical(answer, m.start()):
            continue
        # Tolerate ±1 (some sources count signal peptide differently).
        status = "verified" if abs(n - ref) <= 1 else "contradicted"
        out.append({"claim_type": "protein_length", "claimed": n,
                    "ref": ref, "status": status})
    return out


def _check_charge(answer: str, facts: dict) -> list[dict]:
    out = []
    ch = facts["change"]
    if not ch.get("applicable"):
        return out
    ref = ch["charge_change"]
    for regex, claim_type, expected in T.CHARGE_TRIGGER_COMPILED:
        for m in regex.finditer(answer):
            if _is_hypothetical(answer, m.start()):
                continue
            if _has_pre_disqualifier(answer, m.start(), "charge_change"):
                continue
            status = "verified" if expected == ref else "contradicted"
            out.append({"claim_type": "charge_change", "claimed": expected,
                        "ref": ref, "status": status})
            break  # one fire per pattern is enough
    return out


def _check_polarity(answer: str, facts: dict) -> list[dict]:
    out = []
    ch = facts["change"]
    if not ch.get("applicable"):
        return out
    ref = ch["polarity_change"]
    for regex, claim_type, expected in T.POLARITY_TRIGGER_COMPILED:
        for m in regex.finditer(answer):
            if _is_hypothetical(answer, m.start()):
                continue
            if _has_pre_disqualifier(answer, m.start(), "polarity_change"):
                continue
            status = "verified" if expected == ref else "contradicted"
            out.append({"claim_type": "polarity_change", "claimed": expected,
                        "ref": ref, "status": status})
            break
    return out


def _check_hydrophobicity(answer: str, facts: dict) -> list[dict]:
    out = []
    ch = facts["change"]
    if not ch.get("applicable"):
        return out
    ref = ch["hydrophobicity_change"]
    for regex, claim_type, expected in T.HYDROPHOB_TRIGGER_COMPILED:
        for m in regex.finditer(answer):
            if _is_hypothetical(answer, m.start()):
                continue
            if _has_pre_disqualifier(answer, m.start(), "hydrophobicity_change"):
                continue
            status = "verified" if expected == ref else "contradicted"
            out.append({"claim_type": "hydrophobicity_change", "claimed": expected,
                        "ref": ref, "status": status})
            break
    return out


def _check_size(answer: str, facts: dict) -> list[dict]:
    out = []
    ch = facts["change"]
    if not ch.get("applicable"):
        return out
    ref = ch["size_change"]
    for regex, claim_type, expected in T.SIZE_TRIGGER_COMPILED:
        for m in regex.finditer(answer):
            if _is_hypothetical(answer, m.start()):
                continue
            if _has_pre_disqualifier(answer, m.start(), "size_change"):
                continue
            status = "verified" if expected == ref else "contradicted"
            out.append({"claim_type": "size_change", "claimed": expected,
                        "ref": ref, "status": status})
            break
    return out


def _check_special(answer: str, facts: dict) -> list[dict]:
    """Special-residue claims: only valid if WT or mut AA is in the group's
    `valid_when_aa_in` set. If the trigger fires AND the AA matches, mark
    `verified`; if the trigger fires AND neither WT nor mut is in the set,
    `contradicted`."""
    out = []
    wt = facts["wt_aa"]
    mut = facts["mut_aa"]
    aa_pair = {a for a in (wt, mut) if a and a != "*"}
    if not aa_pair:
        return out
    for group, spec in T.SPECIAL_TRIGGER_COMPILED.items():
        regex = spec["regex"]
        valid_set = spec["valid_when_aa_in"]
        m = regex.search(answer)
        if not m:
            continue
        if _is_hypothetical(answer, m.start()):
            continue
        if _has_pre_disqualifier(answer, m.start(), "special_residue"):
            continue
        if _has_post_disqualifier(answer, m.end(), "special_residue"):
            continue
        applicable = bool(aa_pair & valid_set)
        # Glycosylation needs the sequon.
        if group == "glycosylation":
            applicable = applicable and facts.get("wt_in_glyco_sequon", False)
        status = "verified" if applicable else "contradicted"
        out.append({"claim_type": f"special_{group}",
                    "claimed": True, "ref": applicable, "status": status})
    return out


def _check_region(answer: str, facts: dict) -> list[dict]:
    """N-/C-terminal claims are checked against the position-class fact;
    other regions (active site, TM, etc.) require external annotation we
    don't have, so they fire as `not_applicable` (no row emitted)."""
    out = []
    rc = facts.get("region_class")
    if not rc:
        return out
    for region, regex in T.REGION_TRIGGER_REGEX.items():
        if region not in {"n_terminal", "c_terminal"}:
            continue
        m = regex.search(answer)
        if not m:
            continue
        if _is_hypothetical(answer, m.start()):
            continue
        if _has_pre_disqualifier(answer, m.start(), "region"):
            continue
        if _has_post_disqualifier(answer, m.end(), "region"):
            continue
        status = "verified" if region == rc else "contradicted"
        out.append({"claim_type": "region", "claimed": region,
                    "ref": rc, "status": status})
    return out


def verify_qa(answer: str, facts: dict, feat_store=None) -> list[dict]:
    if not answer or not facts:
        return []
    checks: list[dict] = []
    checks += _check_muttype(answer, facts)
    checks += _check_hgvs_inline(answer, facts)
    checks += _check_position(answer, facts)
    checks += _check_protein_length(answer, facts)
    checks += _check_charge(answer, facts)
    checks += _check_polarity(answer, facts)
    checks += _check_hydrophobicity(answer, facts)
    checks += _check_size(answer, facts)
    checks += _check_special(answer, facts)
    checks += _check_region(answer, facts)
    # UniProt feature checks — gated on having an accession + position +
    # passing the wt-vs-ref-seq sanity check (so we don't fire on rows
    # where the protein sequence doesn't match the claimed wt residue).
    if (feat_store is not None
            and facts.get("uniprot_id")
            and facts.get("position")
            and facts.get("sanity_ok", True)):
        checks += run_uniprot_feature_checks(answer, facts, feat_store)
    return _dedup_checks(checks)


def aggregate_label(checks: list[dict]) -> str:
    if not checks:
        return "not_applicable"
    if any(c["status"] == "contradicted" for c in checks):
        return "contradicted"
    return "verified"


# ------------------------------------------------------------------ Main ---


def iterate_jsonl(path: Path):
    """Yield rows from a `.jsonl` (line-delimited) or `.json` (single list) file."""
    if path.suffix == ".json":
        # Streamed list parse via ijson if available; else load whole file.
        try:
            import ijson  # type: ignore
            with path.open("rb") as f:
                for row in ijson.items(f, "item"):
                    yield row
                return
        except ImportError:
            pass
        data = json.load(path.open())
        for row in data:
            yield row
        return
    with path.open() as f:
        for line in f:
            line = line.strip()
            if line:
                yield json.loads(line)


def _resolve_answer(row: dict, field: str) -> str:
    return row.get(field) or row.get(field.replace("_final", "")) or ""


ANSWER_SOURCES = ("qagen", "llm2")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--input", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--qagen_field", default="qagen_answer_final",
                   help="Field for the first answer (LLM1 / phase1).")
    p.add_argument("--llm2_field", default="llm2_answer_final",
                   help="Field for the second answer (LLM2 / phase2).")
    p.add_argument("--mutation_field", default="HGVSp_Short",
                   help="Field carrying the HGVS-protein descriptor; falls "
                        "back to `mutation` if missing.")
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--features", type=Path, default=FEATURES_CACHE,
                   help="UniProt feature cache from fetch_uniprot_features.py "
                        "(default: audit/output/uniprot_features.json). If "
                        "missing, the uniprot_* checks are skipped.")
    p.add_argument("--no_features", action="store_true",
                   help="Skip the UniProt-feature checks (uniprot_* claim "
                        "families) even if the cache exists.")
    args = p.parse_args()

    src = Path(args.input)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    per_qa_path = out_dir / "mutation_verification.jsonl"

    try:
        if args.no_features:
            raise FileNotFoundError("--no_features given.")
        feat_store = FeatureStore(load_features(args.features))
        print(f"Loaded UniProt features for {sum(1 for _ in feat_store._raw):,} accessions")
    except FileNotFoundError as e:
        print(f"WARNING: {e}\nProceeding without UniProt feature checks.")
        feat_store = None

    bucket_totals = {
        s: defaultdict(
            lambda: Counter({"checkable_qa": 0, "total_qa": 0,
                             "verified_claims": 0, "contradicted_claims": 0,
                             "total_claims": 0})
        )
        for s in ANSWER_SOURCES
    }
    by_claim_type = {s: Counter() for s in ANSWER_SOURCES}
    by_claim_type_verified = {s: Counter() for s in ANSWER_SOURCES}

    parse_fail = 0
    sanity_fail = 0
    total_qa = 0

    source_field = {"qagen": args.qagen_field, "llm2": args.llm2_field}

    with per_qa_path.open("w") as fout:
        for i, row in enumerate(iterate_jsonl(src)):
            if args.limit and i >= args.limit:
                break
            total_qa += 1
            protein = (row.get("protein") or row.get("Hugo_Symbol")
                       or row.get("gene"))
            mutation = row.get(args.mutation_field) or row.get("mutation") or ""
            split = row.get("split") or "unknown"
            qtype = row.get("question_type")
            bucket = bucket_topic(qtype) if qtype else (
                row.get("mutation_type", "unknown") or "unknown").lower()

            ref_seq = row.get("ref_sequence") or ""
            mut_seq = row.get("mutated_sequence") or ""
            mtype = row.get("mutation_type")
            uniprot_id = row.get("uniprot_primary")

            facts = mutation_facts(ref_seq, mut_seq, mutation, mtype, uniprot_id)
            if facts is None:
                parse_fail += 1
                fout.write(json.dumps({
                    "row_id": row.get("_row_id", i), "protein": protein,
                    "mutation": mutation, "parse_failed": True,
                }) + "\n")
                for s in ANSWER_SOURCES:
                    bucket_totals[s][(bucket, split)]["total_qa"] += 1
                continue
            if not facts["sanity_ok"]:
                sanity_fail += 1

            results = {}
            for s in ANSWER_SOURCES:
                ans = _resolve_answer(row, source_field[s])
                ck = verify_qa(ans, facts, feat_store)
                results[s] = {"checks": ck, "label": aggregate_label(ck)}

            out_row = {
                "row_id": row.get("_row_id", i),
                "protein": protein, "mutation": mutation,
                "question_type": qtype, "mutation_type": mtype,
                "bucket": bucket, "split": split,
                "facts": {k: v for k, v in facts.items() if k != "change"},
                "sanity_ok": facts["sanity_ok"],
            }
            for s in ANSWER_SOURCES:
                out_row[f"{s}_checks"] = results[s]["checks"]
                out_row[f"{s}_label"] = results[s]["label"]
            fout.write(json.dumps(out_row, ensure_ascii=False) + "\n")

            for s in ANSWER_SOURCES:
                checks = results[s]["checks"]
                btot = bucket_totals[s][(bucket, split)]
                btot["total_qa"] += 1
                if checks:
                    btot["checkable_qa"] += 1
                    btot["total_claims"] += len(checks)
                    for c in checks:
                        by_claim_type[s][c["claim_type"]] += 1
                        if c["status"] == "verified":
                            btot["verified_claims"] += 1
                            by_claim_type_verified[s][c["claim_type"]] += 1
                        else:
                            btot["contradicted_claims"] += 1

            if (i + 1) % 5000 == 0:
                print(f"  processed {i+1} rows", flush=True)

    # Write summary --------------------------------------------------------
    summary_path = out_dir / "mutation_verification_summary.md"
    lines: list[str] = []
    lines.append("# MutQA Deterministic Verification — Summary (dual-answer)\n")
    lines.append(
        f"Source: `{src}`  \nRows processed: **{total_qa:,}**  "
        f"(HGVS parse failures: {parse_fail}; ref-seq sanity failures: {sanity_fail})\n"
        f"Answer sources: `qagen`={args.qagen_field}, `llm2`={args.llm2_field}\n"
    )

    def _row(vals: Counter):
        total = vals["total_qa"]
        checkable = vals["checkable_qa"]
        verified = vals["verified_claims"]
        contradicted = vals["contradicted_claims"]
        claims = vals["total_claims"]
        cov = checkable / total if total else 0.0
        ver = verified / claims if claims else 0.0
        con = contradicted / claims if claims else 0.0
        return (f"{total:,}", f"{checkable:,}", f"{cov:.3%}",
                f"{claims:,}", f"{ver:.3%}", f"{con:.3%}")

    for s in ANSWER_SOURCES:
        lines.append(f"\n---\n\n# Source: `{s}` ({source_field[s]})\n")
        overall = Counter()
        for v in bucket_totals[s].values():
            overall.update(v)
        lines.append("## Overall\n")
        lines.append("| total Q&A | checkable | coverage | claims | verified% | contradicted% |")
        lines.append("|---|---|---|---|---|---|")
        lines.append("| " + " | ".join(_row(overall)) + " |\n")

        by_bucket = defaultdict(Counter)
        for (bucket, split), v in bucket_totals[s].items():
            by_bucket[bucket].update(v)
        lines.append("## By bucket\n")
        lines.append("| bucket | total Q&A | checkable | coverage | claims | verified% | contradicted% |")
        lines.append("|---|---|---|---|---|---|---|")
        for bucket in sorted(by_bucket.keys()):
            lines.append(f"| {bucket} | " + " | ".join(_row(by_bucket[bucket])) + " |")
        lines.append("")

        by_split = defaultdict(Counter)
        for (bucket, split), v in bucket_totals[s].items():
            by_split[split].update(v)
        lines.append("## By split\n")
        lines.append("| split | total Q&A | checkable | coverage | claims | verified% | contradicted% |")
        lines.append("|---|---|---|---|---|---|---|")
        for split in sorted(by_split.keys()):
            lines.append(f"| {split} | " + " | ".join(_row(by_split[split])) + " |")
        lines.append("")

        lines.append("## Per claim type\n")
        lines.append("| claim_type | total fired | verified | verified% |")
        lines.append("|---|---|---|---|")
        for ct in sorted(by_claim_type[s].keys()):
            total = by_claim_type[s][ct]
            ver = by_claim_type_verified[s][ct]
            rate = ver / total if total else 0.0
            lines.append(f"| {ct} | {total:,} | {ver:,} | {rate:.3%} |")

    summary_path.write_text("\n".join(lines))
    print(f"\nWrote:\n  {per_qa_path}\n  {summary_path}")


if __name__ == "__main__":
    main()
