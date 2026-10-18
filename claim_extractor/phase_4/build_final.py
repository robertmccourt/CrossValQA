#!/usr/bin/env python3
"""
Step 0: Build m2tqa_final.jsonl from the merged Phase 3 output

Joins the merged Phase 3 output (phase_3/merge_recovery.py) with a
(pmid, mutation) -> gene map and writes the Phase 4 starting file that
extract_cited_sentences.py reads.

Transform (reverse-engineered from the original run and verified record-for-record
against the original m2tqa_final.jsonl, 473,243 / 473,243 records identical):
  - every Phase 3 record is kept, in input order, with all fields unchanged
    (no filtering: cross-grounded, grounded AND ungrounded records are kept;
    release filtering happens much later)
  - one field is appended:
      gene: gene_map.get(f"{pmid}|{mutation}", "")   (empty string when unmapped)

Input:
- phase3_validated_merged.jsonl (473,243 records; output of phase_3/merge_recovery.py.
  In the original run this merged file was saved as phase3_validated.jsonl.)
- pmid_mutation_gene_map.json ({"<pmid>|<mutation>": "<GENE>"}, 153,630 keys)

Output:
- m2tqa_final.jsonl (same records with 'gene' field added)

NOTE: pmid_mutation_gene_map.json is a REQUIRED INPUT that is distributed
separately; the code that originally produced it did not survive and is not part
of this repository. Supply it with --gene-map. Format: one JSON object mapping
"<pmid>|<mutation>" (PMID string, '|', the record's mutation string exactly as in
the Phase 3 output) to a gene symbol string, 153,630 keys.
Its 'gene' values are a noisy prior (an early text-mining heuristic) and are NOT
the final gene assignment: downstream steps add recovered_gene
(run_gene_recovery.py), anchor_gene_* (phase_4/anchoring/tiered_anchor_v3.py) and
corrected_gene (uniprot/scripts/correct_gene_assignments.py), which override it.
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

PHASE3_MERGED = config.CLAIM_EXTRACTOR_DIR / "phase_3" / "output" / "phase3_validated_merged.jsonl"
GENE_MAP_PATH = config.CLAIM_EXTRACTOR_DIR / "phase_4" / "output" / "pmid_mutation_gene_map.json"
OUTPUT_PATH = config.CLAIM_EXTRACTOR_DIR / "phase_4" / "output" / "m2tqa_final.jsonl"


def main():
    parser = argparse.ArgumentParser(description="Build m2tqa_final.jsonl from merged Phase 3 output")
    parser.add_argument("--input", default=str(PHASE3_MERGED),
                        help="Merged Phase 3 output (phase_3/merge_recovery.py)")
    parser.add_argument("--gene-map", default=str(GENE_MAP_PATH),
                        help='JSON map {"<pmid>|<mutation>": "<GENE>"}')
    parser.add_argument("--output", default=str(OUTPUT_PATH))
    args = parser.parse_args()

    print("=" * 70)
    print("STEP 0: BUILD m2tqa_final.jsonl")
    print("=" * 70)

    if not Path(args.gene_map).exists():
        sys.exit(f"ERROR: gene map not found: {args.gene_map}\n"
                 "This file is not produced by any script in this repository; supply it with --gene-map.")

    print(f"\nLoading gene map from {args.gene_map}...")
    with open(args.gene_map) as f:
        gene_map = json.load(f)
    print(f"Loaded {len(gene_map):,} (pmid, mutation) -> gene entries")

    print(f"Processing records from {args.input}...")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with_gene = 0
    grounding_counts = Counter()

    with open(args.input) as f_in, open(args.output, "w") as f_out:
        for line in f_in:
            record = json.loads(line)
            key = f"{record['pmid']}|{record['mutation']}"
            record["gene"] = gene_map.get(key, "")
            if record["gene"]:
                with_gene += 1
            grounding_counts[record.get("grounding")] += 1
            f_out.write(json.dumps(record) + "\n")
            count += 1

    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"Records written: {count:,}")
    if count:
        print(f"With gene: {with_gene:,} ({100 * with_gene / count:.1f}%)")
    for grounding, c in grounding_counts.most_common():
        print(f"  {grounding}: {c:,}")
    print(f"\nOutput written to: {args.output}")


if __name__ == "__main__":
    main()
