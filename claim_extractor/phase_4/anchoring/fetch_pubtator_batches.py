#!/usr/bin/env python3
"""
Anchoring step 4 (Tier 1 NER input): fetch gene and variant annotations from
PubTator Central API (network).

Reads batch_*.txt files, extracts PMIDs, fetches annotations from PubTator API,
and saves results in PubTator format (.out files).

Input:
- anchoring/output/ner_input/batch_NNNN.txt (prepare_ner_input.py; only the PMIDs are used)

Output:
- anchoring/output/ner_output/batch_NNNN.out, PubTator text format:
    PMID|t|Title
    PMID|a|Abstract
    PMID<TAB>Start<TAB>End<TAB>Mention<TAB>Type<TAB>ConceptID
  Annotations are PubTator3's (title + abstract offsets), not full-text GNormPlus.
  This is what produced the Tier 1 input of the released dataset; see README.md.

Resumable: batches whose .out file already exists and is non-empty are skipped.
"""

import argparse
import re
import sys
import time
import requests
from pathlib import Path
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import config  # noqa: E402

# Paths
ANCHORING_OUTPUT_DIR = config.PHASE4_DIR / "anchoring" / "output"
INPUT_DIR = ANCHORING_OUTPUT_DIR / "ner_input"
OUTPUT_DIR = ANCHORING_OUTPUT_DIR / "ner_output"

# PubTator API endpoint
PUBTATOR_API = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/pubtator"

# Rate limiting
SLEEP_BETWEEN_REQUESTS = 0.5

# Batch size for API requests (PubTator recommends max 100 per request)
API_BATCH_SIZE = 100


def extract_pmids_from_batch(batch_file: Path) -> list:
    """Extract all unique PMIDs from a batch file."""
    pmids = set()
    with open(batch_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            # PubTator format: PMID|t|Title or PMID|a|Abstract
            match = re.match(r'^(\d+)\|[ta]\|', line)
            if match:
                pmids.add(match.group(1))
    return sorted(pmids)


def fetch_pubtator_annotations(pmids: list) -> str:
    """
    Fetch annotations from PubTator Central API.

    Args:
        pmids: List of PMIDs to fetch

    Returns:
        PubTator format text with annotations
    """
    if not pmids:
        return ""

    try:
        # POST request with JSON body
        response = requests.post(
            PUBTATOR_API,
            json={"pmids": pmids},
            timeout=120
        )
        response.raise_for_status()
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"  Error fetching PMIDs: {e}")
        return ""


def process_batch(batch_file: Path, output_file: Path):
    """Process a single batch file and save annotations."""
    # Extract PMIDs
    pmids = extract_pmids_from_batch(batch_file)

    if not pmids:
        print(f"  No PMIDs found in {batch_file.name}")
        return 0

    print(f"  Found {len(pmids)} PMIDs in {batch_file.name}")

    # Fetch annotations in smaller batches
    all_annotations = []

    for i in range(0, len(pmids), API_BATCH_SIZE):
        batch_pmids = pmids[i:i + API_BATCH_SIZE]

        annotations = fetch_pubtator_annotations(batch_pmids)
        if annotations:
            all_annotations.append(annotations)

        # Rate limiting
        if i + API_BATCH_SIZE < len(pmids):
            time.sleep(SLEEP_BETWEEN_REQUESTS)

    # Write output
    if all_annotations:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("\n".join(all_annotations))
        return len(pmids)

    return 0


def main():
    parser = argparse.ArgumentParser(description="Fetch PubTator3 annotations (PubTator format) per NER input batch")
    parser.add_argument("--input-dir", type=Path, default=INPUT_DIR, help="Directory of batch_*.txt files")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Directory for batch_*.out files")
    args = parser.parse_args()
    input_dir, output_dir = args.input_dir, args.output_dir

    print("=" * 60)
    print("PubTator Central API Annotation Fetcher")
    print("=" * 60)
    print()

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all batch files
    batch_files = sorted(input_dir.glob("batch_*.txt"))

    if not batch_files:
        print(f"No batch files found in {input_dir}")
        return

    print(f"Found {len(batch_files)} batch files")
    print()

    total_pmids = 0
    successful_batches = 0

    for batch_file in tqdm(batch_files, desc="Processing batches"):
        output_file = output_dir / f"{batch_file.stem}.out"

        # Skip if already processed
        if output_file.exists() and output_file.stat().st_size > 0:
            print(f"  Skipping {batch_file.name} (already processed)")
            successful_batches += 1
            continue

        pmid_count = process_batch(batch_file, output_file)

        if pmid_count > 0:
            total_pmids += pmid_count
            successful_batches += 1

        # Rate limiting between batches
        time.sleep(SLEEP_BETWEEN_REQUESTS)

    print()
    print("=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"Successful batches: {successful_batches}/{len(batch_files)}")
    print(f"Total PMIDs processed: {total_pmids}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
