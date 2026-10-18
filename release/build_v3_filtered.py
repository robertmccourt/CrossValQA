"""Apply the combined deterministic-verifier filter, producing the v3 working set.

Drop a row of the v2 working set (v2_verified_answerable) when EITHER:
  - the protein verifier (v3 run) labels qagen_label==contradicted OR
    llm2_label==contradicted, OR
  - the DNA verifier (v3 run) labels qagen_label==contradicted OR
    llm2_label==contradicted.
Rows with parse_failed are never dropped.

Row alignment: the protein verifier ran on rows aligned to the full v2 file
(m2tqa_gene_corrected_v2.jsonl), NOT to v2_verified_answerable, so its row_id
is mapped back to a (pmid, mutation, question) key via v2 row order. The DNA
verifier ran on v2_verified_answerable, so its row_id indexes that file
directly. Both drop sets are unioned and applied by key.

Input:  uniprot/output/m2tqa_gene_corrected_v2.jsonl            (423,415 rows)
        release/output/m2tqa_gene_corrected_v2_verified_answerable.jsonl (411,060)
        protein-verifier labels (mutation_verification.jsonl)
        DNA-verifier labels (mutation_verification_dna.jsonl)
Output: release/output/m2tqa_v3_dna_verified.jsonl               (401,406 rows)

Example:
    python release/build_v3_filtered.py \\
        --protein-labels audit/output/protein_v3/mutation_verification.jsonl \\
        --dna-labels audit/output/dna_v3/mutation_verification_dna.jsonl
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
V2_FULL = REPO_ROOT / "uniprot/output/m2tqa_gene_corrected_v2.jsonl"
V2_WORKING = REPO_ROOT / "release/output/m2tqa_gene_corrected_v2_verified_answerable.jsonl"
OUTPUT = REPO_ROOT / "release/output/m2tqa_v3_dna_verified.jsonl"


def row_key(row):
    """Stable per-QA key. Each (pmid, mutation, question) is unique in v2."""
    return (str(row.get("pmid", row.get("pubmed_id", ""))),
            row.get("mutation", ""),
            row.get("question", ""))


def parse_args():
    ap = argparse.ArgumentParser(description="Drop rows the protein or DNA verifier labels contradicted (v2 working set -> v3).")
    ap.add_argument("--v2", type=Path, default=V2_FULL,
                    help=f"full v2 file the protein verifier is row-aligned to (default: {V2_FULL.relative_to(REPO_ROOT)})")
    ap.add_argument("--working", type=Path, default=V2_WORKING,
                    help=f"v2 working set (default: {V2_WORKING.relative_to(REPO_ROOT)})")
    ap.add_argument("--protein-labels", type=Path, required=True,
                    help="protein-verifier output (v3 run), row-aligned to --v2")
    ap.add_argument("--dna-labels", type=Path, required=True,
                    help="DNA-verifier output (v3 run), row-aligned to --working")
    ap.add_argument("--output", type=Path, default=OUTPUT,
                    help=f"output JSONL (default: {OUTPUT.relative_to(REPO_ROOT)})")
    return ap.parse_args()


def main():
    args = parse_args()
    V2_FULL, V2_WORKING, OUTPUT = args.v2, args.working, args.output
    PROTEIN_VERIFIER, DNA_VERIFIER = args.protein_labels, args.dna_labels
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    # Step 1: build the protein-verifier drop set, keyed by (pmid, mutation, question).
    # The protein verifier output uses row_id from _gold_rows_v2.jsonl, which
    # was row-aligned to v2.jsonl. We map back via v2.jsonl row order.
    print("Loading v2.jsonl row index...")
    v2_keys = []
    with V2_FULL.open() as f:
        for line in f:
            v2_keys.append(row_key(json.loads(line)))
    print(f"  v2 rows: {len(v2_keys):,}")

    print("Reading protein verifier output...")
    drop_protein = set()
    n_protein_seen = 0
    with PROTEIN_VERIFIER.open() as f:
        for line in f:
            r = json.loads(line)
            n_protein_seen += 1
            if r.get("parse_failed"):
                continue
            if (r.get("qagen_label") == "contradicted"
                    or r.get("llm2_label") == "contradicted"):
                rid = r.get("row_id")
                if rid is not None:
                    rid = int(rid)
                    if rid < len(v2_keys):
                        drop_protein.add(v2_keys[rid])
    print(f"  protein verifier rows: {n_protein_seen:,}")
    print(f"  protein-contradicted keys: {len(drop_protein):,}")

    # Step 2: DNA verifier output uses row_id from v2_verified_answerable, so
    # we can read its row_ids directly into a working-row index set.
    print("Loading v2_verified_answerable row index...")
    working_keys = []
    working_rows = []
    with V2_WORKING.open() as f:
        for line in f:
            r = json.loads(line)
            working_keys.append(row_key(r))
            working_rows.append(r)
    print(f"  v2_verified_answerable rows: {len(working_keys):,}")

    print("Reading DNA verifier output...")
    drop_dna_keys = set()
    n_dna_seen = 0
    n_dna_contra = 0
    with DNA_VERIFIER.open() as f:
        for line in f:
            r = json.loads(line)
            n_dna_seen += 1
            if r.get("parse_failed"):
                continue
            if (r.get("qagen_label") == "contradicted"
                    or r.get("llm2_label") == "contradicted"):
                n_dna_contra += 1
                rid = r.get("row_id")
                if rid is not None:
                    rid = int(rid)
                    if rid < len(working_keys):
                        drop_dna_keys.add(working_keys[rid])
    print(f"  DNA verifier rows: {n_dna_seen:,}")
    print(f"  DNA-contradicted rows: {n_dna_contra:,}")
    print(f"  DNA-contradicted keys: {len(drop_dna_keys):,}")

    # Step 3: union the two drop sets and apply to v2_verified_answerable.
    drop_union = drop_protein | drop_dna_keys
    drop_only_protein = drop_protein - drop_dna_keys
    drop_only_dna = drop_dna_keys - drop_protein
    drop_both = drop_protein & drop_dna_keys
    print(f"\n  drop only protein:  {len(drop_only_protein):,}")
    print(f"  drop only DNA:      {len(drop_only_dna):,}")
    print(f"  drop both:          {len(drop_both):,}")
    print(f"  drop union:         {len(drop_union):,}")

    # Step 4: write the filtered file.
    n_in = 0
    n_kept = 0
    n_dropped = 0
    with OUTPUT.open("w") as out:
        for r in working_rows:
            n_in += 1
            if row_key(r) in drop_union:
                n_dropped += 1
                continue
            out.write(json.dumps(r, ensure_ascii=False) + "\n")
            n_kept += 1

    print(f"\nv2_verified_answerable rows in:  {n_in:,}")
    print(f"  kept (v3): {n_kept:,}")
    print(f"  dropped:   {n_dropped:,}  ({100*n_dropped/n_in:.2f}%)")
    print(f"\nWrote: {OUTPUT}")


if __name__ == "__main__":
    main()
