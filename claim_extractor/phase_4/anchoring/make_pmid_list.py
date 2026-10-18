#!/usr/bin/env python3
"""
Anchoring step 1: write the list of unique PMIDs in the Phase 4 output.

Input:
- phase_4/output/m2tqa_gene_recovered.jsonl (run_gene_recovery.py output)

Output:
- anchoring/output/pmid_list.txt (one PMID per line, unique, sorted as strings)

This list drives both the PubTator Central fetch (fetch_pubtator_bioc.py) and
the article batching (prepare_ner_input.py).
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import config  # noqa: E402

QA_FILE = config.PHASE4_GENE_RECOVERED_FILE
ANCHORING_OUTPUT_DIR = config.PHASE4_DIR / "anchoring" / "output"
PMID_FILE = ANCHORING_OUTPUT_DIR / "pmid_list.txt"


def main():
    parser = argparse.ArgumentParser(description="Write unique PMIDs of the Phase 4 output")
    parser.add_argument("--input", type=Path, default=QA_FILE,
                        help="m2tqa_gene_recovered.jsonl (run_gene_recovery.py output)")
    parser.add_argument("--output", type=Path, default=PMID_FILE)
    args = parser.parse_args()

    pmids = set()
    n = 0
    with open(args.input) as f:
        for line in f:
            qa = json.loads(line)
            n += 1
            pmid = str(qa.get("pmid", ""))
            if pmid:
                pmids.add(pmid)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        for pmid in sorted(pmids):
            f.write(pmid + "\n")

    print(f"Records: {n:,}")
    print(f"Unique PMIDs: {len(pmids):,}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
