#!/usr/bin/env python3
"""
Step 1: Extract Cited Sentences

Extracts cited sentence text from articles based on evidence IDs and joins them
with records. This pre-computation step allows gene recovery to run without
needing the full article corpus.

Sentence IDs are resolved with phase_0/prefilter.py::split_into_sentences, the
same splitter that assigned the IDs upstream.

Input:
- m2tqa_final.jsonl (473,243 records; output of build_final.py)
- articles_with_mutations.json (86,317 articles)

Output:
- m2tqa_with_cited.jsonl (same records with 'cited_sentences' field added)
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "phase_0"))
import config  # noqa: E402
from prefilter import split_into_sentences  # noqa: E402

FINAL_DATASET = config.CLAIM_EXTRACTOR_DIR / "phase_4" / "output" / "m2tqa_final.jsonl"
ARTICLE_PATH = config.PIPELINE_ROOT / "database" / "articles_with_mutations.json"
OUTPUT_PATH = config.CLAIM_EXTRACTOR_DIR / "phase_4" / "output" / "m2tqa_with_cited.jsonl"


def extract_cited_sentences(record, article_sentences):
    """Extract cited sentence text from sentence IDs."""
    cited = []
    ev_ids = set()

    for field in ['phase1_evidence_ids', 'phase2_evidence_ids', 'mutation_sentence_ids']:
        ids = record.get(field, [])
        if isinstance(ids, list):
            for item in ids:
                if isinstance(item, int):
                    ev_ids.add(item)
                elif isinstance(item, dict) and 'id' in item:
                    ev_ids.add(item['id'])

    for sid in sorted(ev_ids):
        if 0 <= sid < len(article_sentences):
            cited.append(article_sentences[sid])

    return cited


def main():
    parser = argparse.ArgumentParser(description="Attach cited sentence text to each record")
    parser.add_argument("--input", default=str(FINAL_DATASET), help="m2tqa_final.jsonl (build_final.py output)")
    parser.add_argument("--articles", default=str(ARTICLE_PATH), help="articles_with_mutations.json")
    parser.add_argument("--output", default=str(OUTPUT_PATH))
    args = parser.parse_args()

    print("=" * 70)
    print("STEP 1: EXTRACT CITED SENTENCES")
    print("=" * 70)

    print("\nLoading articles...")
    with open(args.articles) as f:
        articles = json.load(f)
    print(f"Loaded {len(articles)} articles")

    # Pre-split all articles into sentences
    print("Splitting articles into sentences...")
    article_sentences = {}
    for pmid, data in articles.items():
        text = data.get('article', '')
        article_sentences[pmid] = split_into_sentences(text) if text else []

    print(f"Processing records from {args.input}...")
    count = 0
    missing_pmid = 0
    total_cited = 0

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.input) as f_in, open(args.output, 'w') as f_out:
        for line in f_in:
            record = json.loads(line)
            pmid = str(record.get('pmid', ''))

            if pmid in article_sentences:
                cited = extract_cited_sentences(record, article_sentences[pmid])
                record['cited_sentences'] = cited
                total_cited += len(cited)
            else:
                record['cited_sentences'] = []
                missing_pmid += 1

            f_out.write(json.dumps(record) + '\n')
            count += 1

            if count % 50000 == 0:
                print(f"  Processed {count}...")

    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"Records processed: {count}")
    print(f"Missing PMIDs: {missing_pmid}")
    print(f"Total cited sentences extracted: {total_cited}")
    print(f"Avg cited sentences per record: {total_cited / count:.2f}")
    print(f"\nOutput written to: {args.output}")


if __name__ == "__main__":
    main()
