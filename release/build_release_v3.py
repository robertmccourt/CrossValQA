"""Build the 39-column public-release file from the v3 filtered working set.

Joins, by (pmid, mutation, question):
  - the v3 working set (m2tqa_v3_dna_verified.jsonl, 401,406 rows), which
    supplies the phase / cross-validation / traceability columns, and
  - the enriched v2 release (crossvalqa_release_v2.jsonl, 32 columns), which
    supplies sequences, parsed mutation fields, isoform and 3D-structure
    columns.
Working-set values win when present and non-null; otherwise the enriched
value is used. Rows with no enriched match keep null structural fields.

Adds the derived `notation_type` column (protein / dna / other).
Drops `qagen_answer` from the enriched file (it is identical to phase2_answer).
Projects to the 39 columns in the 8 ordered groups of COLUMN_ORDER.

Input:  release/output/m2tqa_v3_dna_verified.jsonl
        release/output/crossvalqa_release_v2.jsonl
Output: release/output/crossvalqa_release_v3.jsonl   (401,406 rows)

Example:
    python release/build_release_v3.py
    python release/build_release_v3.py --working W.jsonl --enriched E.jsonl --output O.jsonl
"""
from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "release/output"
WORKING = OUTPUT_DIR / "m2tqa_v3_dna_verified.jsonl"
ENRICHED = OUTPUT_DIR / "crossvalqa_release_v2.jsonl"
OUTPUT = OUTPUT_DIR / "crossvalqa_release_v3.jsonl"

# 39 columns in 8 ordered groups.
COLUMN_ORDER = [
    # Group 1 IDENTITY (6)
    "pmid", "gene", "mutation", "notation_type", "question", "question_type",
    # Group 2 PHASE 1 (3)
    "phase1_answer", "phase1_rationale", "phase1_evidence_ids",
    # Group 3 PHASE 2 (3)
    "phase2_answer", "phase2_rationale", "phase2_evidence_ids",
    # Group 4 CROSS-VALIDATION (6)
    "grounding", "grounding_method",
    "r1_supports_a2", "r1_a2_reasoning",
    "r2_supports_a1", "r2_a1_reasoning",
    # Group 5 PROTEIN CONTEXT (10)
    "uniprot_primary", "ref_sequence", "mutated_sequence",
    "wt_aa", "mut_aa", "mutation_position", "mutation_locations",
    "matched_isoform", "matched_isoform_status", "matched_isoform_reason",
    # Group 6 3D STRUCTURE (6)
    "has_structure", "structure_evaluable",
    "alphafold_id", "alphafold_url",
    "pdb_ids", "pdb_urls",
    # Group 7 SOURCE TRACEABILITY (3)
    "cited_sentences", "context_sentence_ids", "mutation_sentence_ids",
    # Group 8 AUXILIARY (2)
    "clinvar_variant_name", "split",
]
assert len(COLUMN_ORDER) == 39, f"expected 39 cols, got {len(COLUMN_ORDER)}"

DNA_PREFIX = re.compile(r"^[cgmn]\.", re.IGNORECASE)
BARE_NUC = re.compile(r"^\d+[ACGT]>[ACGT]$|^[ACGT]\d+[ACGT]$")
PROT_PREFIX = re.compile(r"^p\.")


def notation_type_of(mutation):
    if not mutation:
        return "other"
    if PROT_PREFIX.match(mutation):
        return "protein"
    if DNA_PREFIX.match(mutation) or BARE_NUC.match(mutation):
        return "dna"
    return "other"


def row_key(row):
    return (str(row.get("pmid", "")),
            row.get("mutation", ""),
            row.get("question", ""))


def parse_args():
    ap = argparse.ArgumentParser(description="Project the v3 working set + enriched v2 release to the 39-column release.")
    ap.add_argument("--working", type=Path, default=WORKING,
                    help=f"v3 working set (default: {WORKING.relative_to(REPO_ROOT)})")
    ap.add_argument("--enriched", type=Path, default=ENRICHED,
                    help=f"enriched v2 release (default: {ENRICHED.relative_to(REPO_ROOT)})")
    ap.add_argument("--output", type=Path, default=OUTPUT,
                    help=f"output JSONL (default: {OUTPUT.relative_to(REPO_ROOT)})")
    return ap.parse_args()


def main():
    args = parse_args()
    WORKING, ENRICHED, OUTPUT = args.working, args.enriched, args.output
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    print(f"Reading enriched release: {ENRICHED}")
    enriched_idx = {}
    with ENRICHED.open() as f:
        for line in f:
            r = json.loads(line)
            enriched_idx[row_key(r)] = r
    print(f"  {len(enriched_idx):,} enriched rows indexed")

    print(f"\nReading working v3: {WORKING}")
    n_in = 0
    n_matched = 0
    n_unmatched = 0
    notation_counts = Counter()
    missing_field_counts = Counter()
    with WORKING.open() as wf, OUTPUT.open("w") as out:
        for line in wf:
            n_in += 1
            r = json.loads(line)
            key = row_key(r)
            enriched = enriched_idx.get(key)
            if enriched is None:
                n_unmatched += 1
                # Fall back to working-only fields; structural ones will be null
                enriched = {}
            else:
                n_matched += 1

            # Build the merged row preferring working fields for the
            # phase/cross-val columns and enriched for the structural ones.
            merged = {}
            mutation = r.get("mutation", "")
            ntype = notation_type_of(mutation)
            notation_counts[ntype] += 1

            for col in COLUMN_ORDER:
                if col == "notation_type":
                    merged[col] = ntype
                    continue
                if col == "gene":
                    # enriched.gene == corrected_gene (checked on a sample).
                    # Prefer enriched; fall back to working corrected_gene.
                    merged[col] = (enriched.get("gene")
                                   or r.get("corrected_gene")
                                   or r.get("gene"))
                    continue
                # Working-side fields (phase rationales, cross-val reasoning,
                # source-traceability, etc.) usually only in working file.
                # Structural fields (sequences, 3D, isoform) only in enriched.
                # Try working first, fall back to enriched for everything,
                # since the field name is unambiguous.
                if col in r and r[col] is not None:
                    merged[col] = r[col]
                elif col in enriched:
                    merged[col] = enriched[col]
                else:
                    merged[col] = None
                    missing_field_counts[col] += 1

            out.write(json.dumps(merged, ensure_ascii=False) + "\n")

    print(f"\nRows written:    {n_in:,}")
    print(f"  matched in enriched:   {n_matched:,}")
    print(f"  unmatched (no merge):  {n_unmatched:,}")
    print(f"\nNotation type distribution:")
    for nt, c in notation_counts.most_common():
        print(f"  {nt}: {c:,} ({100*c/n_in:.1f}%)")
    if missing_field_counts:
        print(f"\nFields missing in some rows (working OR enriched):")
        for col, c in missing_field_counts.most_common():
            print(f"  {col}: {c:,} ({100*c/n_in:.1f}%)")
    print(f"\nWrote: {OUTPUT}")
    print(f"  size: {OUTPUT.stat().st_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
