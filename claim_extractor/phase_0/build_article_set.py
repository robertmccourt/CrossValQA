#!/usr/bin/env python3
"""
Build the Phase 0 article set (articles_with_mutations.json).

Filters the output of batch_extractor.py down to the articles used by
verify_mutations.py and Phases 1-2. Filters are applied in this order:
  1. Drop articles with > MAX_MUTATIONS extracted mutations (regex outliers,
     e.g. large variant tables).
  2. Drop articles whose text is shorter than MIN_ARTICLE_CHARS (empty or
     abstract-only records).
  3. Drop articles with zero extracted mutations.

Phase 0 run order:
    python batch_extractor.py        # source collection -> *.mutations_extracted.json
    python build_article_set.py      # -> articles_with_mutations.json
    python verify_mutations.py       # -> verified_mutations.json

Usage:
    python build_article_set.py [--input PATH] [--output PATH] [--no-mutations-output PATH]
"""

import argparse
import json
import sys
from pathlib import Path

# Add parent directory for config
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import EXTRACTED_ARTICLES_FILE, ARTICLES_FILE

MAX_MUTATIONS = 200       # articles with more extracted mutations are removed
MIN_ARTICLE_CHARS = 5000  # articles with shorter text are removed


def main():
    parser = argparse.ArgumentParser(description="Filter extracted articles into the Phase 0 article set")
    parser.add_argument(
        "--input",
        default=str(EXTRACTED_ARTICLES_FILE),
        help="Input JSON from batch_extractor.py (PMID -> {article, proteins})"
    )
    parser.add_argument(
        "--output",
        default=str(ARTICLES_FILE),
        help="Output JSON file path (articles_with_mutations.json)"
    )
    parser.add_argument(
        "--no-mutations-output",
        default=None,
        help="Optional path to save articles removed for having zero mutations (not used downstream)"
    )
    args = parser.parse_args()

    print(f"Loading data from {args.input}...")
    with open(args.input, "r") as f:
        data = json.load(f)
    original_count = len(data)
    print(f"Original articles: {original_count:,}")

    # 1. Mutation-count outliers
    removed = 0
    for pmid in list(data.keys()):
        if len(data[pmid].get("proteins", [])) > MAX_MUTATIONS:
            del data[pmid]
            removed += 1
    print(f"Removed: {removed:,} articles with > {MAX_MUTATIONS} mutations")
    print(f"Remaining: {len(data):,} articles")

    # 2. Short / empty articles
    removed = 0
    for pmid in list(data.keys()):
        article = data[pmid].get("article") or ""
        if len(article) < MIN_ARTICLE_CHARS:
            del data[pmid]
            removed += 1
    print(f"Removed: {removed:,} articles with < {MIN_ARTICLE_CHARS:,} chars")
    print(f"Remaining: {len(data):,} articles")

    # 3. Articles with no extracted mutations
    with_mutations = {}
    without_mutations = {}
    for pmid, entry in data.items():
        if len(entry.get("proteins", [])) > 0:
            with_mutations[pmid] = entry
        else:
            without_mutations[pmid] = entry
    print(f"Articles with mutations: {len(with_mutations):,}")
    print(f"Articles without mutations (removed): {len(without_mutations):,}")

    print(f"\nSaving {len(with_mutations):,} articles to {args.output}...")
    with open(args.output, "w") as f:
        json.dump(with_mutations, f)

    if args.no_mutations_output:
        print(f"Saving {len(without_mutations):,} articles to {args.no_mutations_output}...")
        with open(args.no_mutations_output, "w") as f:
            json.dump(without_mutations, f)
    print("Saved.")


if __name__ == "__main__":
    main()
