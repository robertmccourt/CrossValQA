#!/usr/bin/env python3
"""Prepare the CrossValQA v2 release file (step 1 of 4): projection + ref_sequence.

Reads `m2tqa_gene_corrected_v2.jsonl`, projects every record to the lean
release schema (V1_LEAN_FIELDS), normalizes the 10 straggler question_type
variants that the earlier normalizer missed, and populates `ref_sequence`
via (a) a uniprot->sequence cache built from the v1 release, (b) a local
fetch cache, and (c) a UniProt REST fetch for accessions still missing.
`gene` is taken from `corrected_gene`; `qagen_answer` is carried over from
the v1 release by (pmid, mutation, question).

All other enrichment fields (pdb/alphafold, mutation_position, wt_aa/mut_aa,
has_structure, structure_evaluable, mutated_sequence) are left as null/empty
placeholders; enrich_release_v2.py fills them.

Input:  uniprot/output/m2tqa_gene_corrected_v2.jsonl  (423,415 rows)
        release/input/crossvalqa_release.jsonl        (v1 enriched release, 416K rows;
                                                       built outside this repo)
        UniProt REST (rest.uniprot.org) unless --skip-fetch
Output: release/output/crossvalqa_release_v2.jsonl
        release/output/uniprot_ref_sequences.json     (fetch cache)
        release/output/uniprot_unfetched.txt

Example:
    python release/prepare_release_v2.py --limit 1000 --skip-fetch   # pilot
    python release/prepare_release_v2.py
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.request
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = REPO_ROOT / "release/input"
OUTPUT_DIR = REPO_ROOT / "release/output"
V1_RELEASE = INPUT_DIR / "crossvalqa_release.jsonl"
V2_SRC = REPO_ROOT / "uniprot/output/m2tqa_gene_corrected_v2.jsonl"
DST = OUTPUT_DIR / "crossvalqa_release_v2.jsonl"
SEQ_CACHE = OUTPUT_DIR / "uniprot_ref_sequences.json"
UNFETCHED_REPORT = OUTPUT_DIR / "uniprot_unfetched.txt"

QTYPE_STRAGGLER_MAP = {
    "functional": "function",
    "functional effect": "function",
    "functional_effect": "function",
    "clinical features": "clinical",
    "clinical_features": "clinical",
    "frequency/prevalence": "frequency",
    "disease association": "disease",
    "disease_association": "disease",
    "pathogenicity|disease": "pathogenicity",
    "frequency|clinical": "frequency",
}
CANONICAL_QTYPES = {
    "pathogenicity", "other", "function", "mechanism",
    "disease", "frequency", "clinical", "inheritance",
}

# Lean v1 release schema. Fields not in this set are stripped from v2.
V1_LEAN_FIELDS = [
    "pmid", "gene", "mutation", "question", "question_type",
    "qagen_answer", "phase1_answer", "phase2_answer",
    "phase2_confidence", "phase2_answerable",
    "grounding", "grounding_method", "uniprot_primary",
    "cited_sentences", "phase1_evidence_ids", "phase2_evidence_ids",
    "ref_sequence", "mutated_sequence",
    "split",
    "pdb_ids", "pdb_urls", "alphafold_id", "alphafold_url",
    "mutation_position", "mutation_locations",
    "wt_aa", "mut_aa",
    "has_structure", "structure_evaluable",
]


def load_v1_caches():
    """Build (uniprot->sequence) and ((pmid,mutation,question)->qagen_answer) caches from the v1 release."""
    seq_cache: dict[str, str] = {}
    qagen_cache: dict[tuple, str] = {}
    with open(V1_RELEASE) as f:
        for line in f:
            r = json.loads(line)
            up = r.get("uniprot_primary")
            seq = r.get("ref_sequence")
            if up and seq and up not in seq_cache:
                seq_cache[up] = seq
            qa = r.get("qagen_answer")
            if qa is not None:
                key = (r.get("pmid"), r.get("mutation"), r.get("question"))
                qagen_cache[key] = qa
    return seq_cache, qagen_cache


def collect_needed_uniprots(seq_cache: dict[str, str]) -> list[str]:
    needed: set[str] = set()
    with open(V2_SRC) as f:
        for line in f:
            r = json.loads(line)
            up = r.get("uniprot_primary")
            if up and up not in seq_cache:
                needed.add(up)
    return sorted(needed)


def fetch_uniprot_batch(accessions: list[str]) -> dict[str, str]:
    """Fetch FASTA for a batch of UniProt accessions. Returns dict acc->seq."""
    if not accessions:
        return {}
    url = (
        "https://rest.uniprot.org/uniprotkb/accessions"
        f"?accessions={','.join(accessions)}&format=fasta&size=500"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "CrossValQA-release-prep/1.0"})
    result: dict[str, str] = {}
    with urllib.request.urlopen(req, timeout=60) as resp:
        text = resp.read().decode("utf-8", errors="replace")
    current_acc = None
    current_seq: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if current_acc is not None and current_seq:
                result[current_acc] = "".join(current_seq)
            # Parse header e.g. ">sp|Q8BGQ1|VPS33B_MOUSE ..." → Q8BGQ1
            parts = line[1:].split("|")
            current_acc = parts[1] if len(parts) >= 2 else None
            current_seq = []
        else:
            if current_acc is not None:
                current_seq.append(line.strip())
    if current_acc is not None and current_seq:
        result[current_acc] = "".join(current_seq)
    return result


def fetch_missing_sequences(needed: list[str], batch_size: int = 200, sleep_s: float = 0.5) -> dict[str, str]:
    fetched: dict[str, str] = {}
    total = len(needed)
    print(f"Fetching {total:,} uniprot sequences (batch={batch_size})")
    for i in range(0, total, batch_size):
        batch = needed[i:i + batch_size]
        try:
            got = fetch_uniprot_batch(batch)
        except Exception as exc:
            print(f"  batch {i}-{i + len(batch)}: ERROR {exc}", file=sys.stderr)
            got = {}
        fetched.update(got)
        print(f"  {i + len(batch):>6}/{total:>6}  (+{len(got):>3} in batch, {len(fetched):>5} cum)")
        if i + batch_size < total:
            time.sleep(sleep_s)
    return fetched


def normalize_qtype(qt: str | None) -> str | None:
    if qt is None:
        return None
    if qt in CANONICAL_QTYPES:
        return qt
    mapped = QTYPE_STRAGGLER_MAP.get(qt)
    if mapped:
        return mapped
    return qt  # leave exotic values alone; reported in summary


def project_to_v1_schema(
    v2_rec: dict,
    seq_cache: dict[str, str],
    qagen_cache: dict[tuple, str],
) -> tuple[dict, bool]:
    """Return (lean v1-schema record, qtype_was_renormalized)."""
    up = v2_rec.get("uniprot_primary")
    original_qt = v2_rec.get("question_type")
    norm_qt = normalize_qtype(original_qt)
    renormalized = (norm_qt != original_qt) and (norm_qt is not None)
    key = (v2_rec.get("pmid"), v2_rec.get("mutation"), v2_rec.get("question"))
    out = {
        "pmid": v2_rec.get("pmid"),
        # v1 release `gene` is the corrected/recovered symbol; v2 keeps the raw
        # symbol in `gene` and the corrected one in `corrected_gene`.
        "gene": v2_rec.get("corrected_gene") or v2_rec.get("gene"),
        "mutation": v2_rec.get("mutation"),
        "question": v2_rec.get("question"),
        "question_type": norm_qt,
        "qagen_answer": qagen_cache.get(key),
        "phase1_answer": v2_rec.get("phase1_answer"),
        "phase2_answer": v2_rec.get("phase2_answer"),
        "phase2_confidence": v2_rec.get("phase2_confidence"),
        "phase2_answerable": v2_rec.get("phase2_answerable"),
        "grounding": v2_rec.get("grounding"),
        "grounding_method": v2_rec.get("grounding_method"),
        "uniprot_primary": up,
        "cited_sentences": v2_rec.get("cited_sentences"),
        "phase1_evidence_ids": v2_rec.get("phase1_evidence_ids"),
        "phase2_evidence_ids": v2_rec.get("phase2_evidence_ids"),
        "ref_sequence": seq_cache.get(up) if up else None,
        "mutated_sequence": None,   # step 2
        "split": None,              # create_splits.py writes this
        "pdb_ids": [],              # step 2
        "pdb_urls": [],             # step 2
        "alphafold_id": None,       # step 2
        "alphafold_url": None,      # step 2
        "mutation_position": None,  # step 2
        "mutation_locations": None, # step 2
        "wt_aa": None,              # step 2
        "mut_aa": None,             # step 2
        "has_structure": None,      # step 2
        "structure_evaluable": None,  # step 2
    }
    return out, renormalized


def main(argv: list[str] | None = None) -> int:
    global V1_RELEASE, V2_SRC, DST, SEQ_CACHE, UNFETCHED_REPORT
    ap = argparse.ArgumentParser(description="Project v2 to the release schema and fill ref_sequence.")
    ap.add_argument("--limit", type=int, default=None, help="Pilot: process only the first N records")
    ap.add_argument("--skip-fetch", action="store_true", help="Skip UniProt REST fetch; use only the v1 cache")
    ap.add_argument("--use-seq-cache", action="store_true", help="Load prior fetch cache if it exists (default on)")
    ap.add_argument("--v1-release", type=Path, default=V1_RELEASE, help=f"v1 enriched release (default: {V1_RELEASE.relative_to(REPO_ROOT)})")
    ap.add_argument("--v2", type=Path, default=V2_SRC, help=f"v2 dataset (default: {V2_SRC.relative_to(REPO_ROOT)})")
    ap.add_argument("--output", type=Path, default=DST, help=f"output JSONL (default: {DST.relative_to(REPO_ROOT)})")
    ap.add_argument("--seq-cache", type=Path, default=SEQ_CACHE, help=f"UniProt sequence fetch cache (default: {SEQ_CACHE.relative_to(REPO_ROOT)})")
    ap.add_argument("--unfetched-report", type=Path, default=UNFETCHED_REPORT, help=f"list of accessions not fetched (default: {UNFETCHED_REPORT.relative_to(REPO_ROOT)})")
    args = ap.parse_args(argv)
    V1_RELEASE, V2_SRC, DST = args.v1_release, args.v2, args.output
    SEQ_CACHE, UNFETCHED_REPORT = args.seq_cache, args.unfetched_report

    print(f"[1/5] Loading v1 caches from {V1_RELEASE}")
    seq_cache, qagen_cache = load_v1_caches()
    print(f"  seq cache: {len(seq_cache):,} uniprots")
    print(f"  qagen cache: {len(qagen_cache):,} (pmid,mutation,question) keys")

    if SEQ_CACHE.exists():
        prior = json.loads(SEQ_CACHE.read_text())
        added = 0
        for acc, seq in prior.items():
            if acc not in seq_cache:
                seq_cache[acc] = seq
                added += 1
        print(f"  loaded {added:,} additional sequences from {SEQ_CACHE}")

    print(f"[2/5] Scanning v2 for missing uniprots")
    needed = collect_needed_uniprots(seq_cache)
    print(f"  {len(needed):,} uniprots need fetching")

    if needed and not args.skip_fetch:
        fetched = fetch_missing_sequences(needed)
        # merge into cache
        for acc, seq in fetched.items():
            seq_cache[acc] = seq
        # persist the newly-fetched portion so reruns are cheap
        SEQ_CACHE.parent.mkdir(parents=True, exist_ok=True)
        prior = {}
        if SEQ_CACHE.exists():
            prior = json.loads(SEQ_CACHE.read_text())
        prior.update(fetched)
        SEQ_CACHE.write_text(json.dumps(prior))
        print(f"  cached {len(prior):,} fetched sequences to {SEQ_CACHE}")
        still_missing = [acc for acc in needed if acc not in seq_cache]
        UNFETCHED_REPORT.parent.mkdir(parents=True, exist_ok=True)
        UNFETCHED_REPORT.write_text("\n".join(still_missing))
        print(f"  {len(still_missing):,} uniprots still unfetched -> {UNFETCHED_REPORT}")
    elif args.skip_fetch:
        print("  --skip-fetch: using cache only")

    print(f"[3/5] Streaming v2 -> v1-schema JSONL")
    DST.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    renormalized_count = 0
    exotic_qtypes: dict[str, int] = {}
    grounding_counts: dict[str, int] = {}
    have_ref = 0
    with open(V2_SRC) as fin, open(DST, "w") as fout:
        for line in fin:
            if args.limit is not None and total >= args.limit:
                break
            r = json.loads(line)
            out, renormalized = project_to_v1_schema(r, seq_cache, qagen_cache)
            if renormalized:
                renormalized_count += 1
            qt = out["question_type"]
            if qt not in CANONICAL_QTYPES and qt is not None:
                exotic_qtypes[qt] = exotic_qtypes.get(qt, 0) + 1
            grounding_counts[out["grounding"]] = grounding_counts.get(out["grounding"], 0) + 1
            if out["ref_sequence"] is not None:
                have_ref += 1
            fout.write(json.dumps(out) + "\n")
            total += 1

    print(f"[4/5] Summary")
    print(f"  wrote {total:,} records to {DST}")
    print(f"  qtype re-normalized (stragglers caught): {renormalized_count:,}")
    if exotic_qtypes:
        print(f"  !! {sum(exotic_qtypes.values()):,} records still have non-canonical qtype:")
        for k, v in sorted(exotic_qtypes.items(), key=lambda x: -x[1]):
            print(f"       {v:>6}  {k!r}")
    else:
        print(f"  all qtypes in canonical 8")
    print(f"  grounding distribution: {grounding_counts}")
    print(f"  ref_sequence populated: {have_ref:,} / {total:,} ({100*have_ref/max(total,1):.1f}%)")

    print(f"[5/5] Verify: head -1 | json.tool")
    with open(DST) as f:
        first = json.loads(f.readline())
    for k in V1_LEAN_FIELDS:
        if k not in first:
            print(f"  MISSING field in output: {k}", file=sys.stderr)
            return 1
    print(f"  all {len(V1_LEAN_FIELDS)} v1-schema fields present on first record")
    return 0


if __name__ == "__main__":
    sys.exit(main())
