#!/usr/bin/env python3
"""Build the v2 working set: drop v2-verifier contradictions and unanswerable rows.

Two filters over `m2tqa_gene_corrected_v2.jsonl` (423,415 rows):
  1. Drop a row when the deterministic protein verifier (first, v2 run) labels
     `qagen_label == "contradicted"` OR `llm2_label == "contradicted"`.
     The verifier's `row_id` is the 0-based line index in the v2 file;
     rows with `parse_failed` are kept.                        -> 416,101 rows
  2. Drop a row when `phase2_answerable is False` (the Phase 2 model
     explicitly declined to answer).                          -> 411,060 rows

Kept lines are written byte-for-byte unchanged, in input order. This
reproduces `m2tqa_gene_corrected_v2_verified_answerable.jsonl`, the input to
build_v3_filtered.py and to the DNA verifier.

Input:  uniprot/output/m2tqa_gene_corrected_v2.jsonl
        + protein-verifier labels (v2 run), e.g. mutation_verification.jsonl
Output: release/output/m2tqa_gene_corrected_v2_verified_answerable.jsonl

Example:
    python release/build_v2_verified_answerable.py \\
        --protein-labels audit/output/protein_v2/mutation_verification.jsonl
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
V2_FULL = REPO_ROOT / "uniprot/output/m2tqa_gene_corrected_v2.jsonl"
OUTPUT = REPO_ROOT / "release/output/m2tqa_gene_corrected_v2_verified_answerable.jsonl"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--v2", type=Path, default=V2_FULL, help=f"v2 dataset (default: {V2_FULL.relative_to(REPO_ROOT)})")
    ap.add_argument("--protein-labels", type=Path, required=True,
                    help="protein-verifier output for the v2 run (JSONL with row_id, qagen_label, llm2_label)")
    ap.add_argument("--output", type=Path, default=OUTPUT, help=f"output JSONL (default: {OUTPUT.relative_to(REPO_ROOT)})")
    args = ap.parse_args()

    print(f"Reading protein verifier labels: {args.protein_labels}")
    drop_rows = set()
    n_seen = 0
    with args.protein_labels.open() as f:
        for line in f:
            r = json.loads(line)
            n_seen += 1
            if r.get("parse_failed"):
                continue
            if (r.get("qagen_label") == "contradicted"
                    or r.get("llm2_label") == "contradicted"):
                drop_rows.add(int(r["row_id"]))
    print(f"  verifier rows: {n_seen:,}")
    print(f"  contradicted rows: {len(drop_rows):,}")

    n_in = n_contra = n_unanswerable = n_kept = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.v2.open() as fin, args.output.open("w") as fout:
        for i, line in enumerate(fin):
            n_in += 1
            if i in drop_rows:
                n_contra += 1
                continue
            if json.loads(line).get("phase2_answerable") is False:
                n_unanswerable += 1
                continue
            fout.write(line)
            n_kept += 1

    print(f"\nv2 rows in:                 {n_in:,}")
    print(f"  dropped (contradicted):   {n_contra:,}")
    print(f"  after filter 1:           {n_in - n_contra:,}")
    print(f"  dropped (unanswerable):   {n_unanswerable:,}")
    print(f"  kept:                     {n_kept:,}")
    print(f"\nWrote: {args.output}")


if __name__ == "__main__":
    main()
