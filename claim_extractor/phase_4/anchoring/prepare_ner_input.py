#!/usr/bin/env python3
"""
Anchoring step 3: convert articles to PubTator format for GNormPlus processing.

GNormPlus accepts PubTator format:
PMID|t|Title text
PMID|a|Abstract/body text

For full text, we concatenate into the abstract field.

Input:
- the source article corpus, {pmid: {"article": <full text>, ...}}
  (config.SOURCE_ARTICLES_FILE; the same file Phase 0 starts from)
- anchoring/output/pmid_list.txt (make_pmid_list.py)

Output:
- anchoring/output/ner_input/batch_NNNN.txt (BATCH_SIZE articles each) + manifest.txt

These batch files are the input of the Tier 1 NER step: either GNormPlus run on the
full text (see README.md) or fetch_pubtator_batches.py, which reads only the PMIDs
from them and downloads PubTator3 annotations instead.
"""

import argparse
import json
import os
import sys
from pathlib import Path
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import config  # noqa: E402

ANCHORING_OUTPUT_DIR = config.PHASE4_DIR / "anchoring" / "output"
ARTICLES_FILE = config.SOURCE_ARTICLES_FILE
OUTPUT_DIR = ANCHORING_OUTPUT_DIR / "ner_input"
PMID_FILE = ANCHORING_OUTPUT_DIR / "pmid_list.txt"

# GNormPlus performs better with smaller batches
BATCH_SIZE = 1000  # Articles per file


def load_target_pmids(pmid_file):
    """Load PMIDs we need to process."""
    with open(pmid_file) as f:
        return set(line.strip() for line in f if line.strip())


def clean_text(text: str) -> str:
    """Clean text for PubTator format."""
    if not text:
        return ""
    # Replace problematic characters
    text = text.replace('\n', ' ').replace('\r', ' ')
    text = text.replace('\t', ' ')
    text = text.replace('|', ' ')  # Pipe is delimiter
    # Normalize whitespace
    text = ' '.join(text.split())
    return text


def extract_title_and_body(full_text: str) -> tuple:
    """Extract title and body from full text."""
    if not full_text:
        return "", ""

    lines = full_text.split('\n')

    # First non-empty line is typically the title
    title = ""
    body_start = 0
    for i, line in enumerate(lines):
        line = line.strip()
        if line:
            title = line
            body_start = i + 1
            break

    # Rest is body
    body = ' '.join(lines[body_start:])

    return clean_text(title), clean_text(body)


def main():
    parser = argparse.ArgumentParser(description="Write articles as PubTator-format batch files for NER")
    parser.add_argument("--articles", type=Path, default=ARTICLES_FILE,
                        help='Source article corpus JSON {pmid: {"article": ...}}')
    parser.add_argument("--pmid-file", type=Path, default=PMID_FILE, help="pmid_list.txt (make_pmid_list.py)")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()

    print("Loading target PMIDs...")
    target_pmids = load_target_pmids(args.pmid_file)
    print(f"Target PMIDs: {len(target_pmids):,}")

    print("Loading articles...")
    with open(args.articles) as f:
        all_articles = json.load(f)
    print(f"Total articles in file: {len(all_articles):,}")

    # Filter to target PMIDs
    articles = {}
    for pmid, data in all_articles.items():
        if pmid in target_pmids:
            if isinstance(data, dict) and data.get('article'):
                articles[pmid] = data['article']

    print(f"Articles with text for target PMIDs: {len(articles):,}")

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Write in batches
    pmids = list(articles.keys())
    batch_num = 0
    total_written = 0

    for i in range(0, len(pmids), BATCH_SIZE):
        batch_pmids = pmids[i:i+BATCH_SIZE]
        batch_file = args.output_dir / f"batch_{batch_num:04d}.txt"

        with open(batch_file, 'w', encoding='utf-8') as f:
            for pmid in batch_pmids:
                full_text = articles[pmid]
                title, body = extract_title_and_body(full_text)

                if title or body:
                    f.write(f"{pmid}|t|{title}\n")
                    f.write(f"{pmid}|a|{body}\n")
                    f.write("\n")  # Blank line separator
                    total_written += 1

        batch_num += 1

    print(f"\nCreated {batch_num} batch files in {args.output_dir}")
    print(f"Total articles written: {total_written:,}")

    # Also create a manifest
    manifest_file = args.output_dir / "manifest.txt"
    with open(manifest_file, 'w') as f:
        for i in range(batch_num):
            f.write(f"batch_{i:04d}.txt\n")

    print(f"Manifest: {manifest_file}")


if __name__ == "__main__":
    main()
