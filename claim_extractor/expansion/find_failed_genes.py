#!/usr/bin/env python3
"""
Identify genes that failed to produce grounded QA pairs from a pipeline run.

Compares the target gene list against Phase 3 output to find genes that either:
1. Had no articles with extractable mutations
2. Had mutations but failed verification/gene-match
3. Passed through phases but produced no grounded QA pairs

Outputs a targets JSONL file for re-fetching with more papers per gene
(this is how the batch2_round2 targets were built).

Input:  expansion/output/<batch articles>.json + <batch>/phase3/phase3_validated.jsonl
Output: expansion/output/<batch>_failed_genes.jsonl

Usage:
    # Find failed genes from Batch 2
    python find_failed_genes.py --batch batch2

    # Custom paths
    python find_failed_genes.py \
        --articles expansion/output/batch2_articles.json \
        --phase3 expansion/output/batch2/phase3/phase3_validated.jsonl \
        --output expansion/output/batch2_failed_genes.jsonl
"""

import json
import argparse
from pathlib import Path
from collections import defaultdict

EXPANSION_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = EXPANSION_DIR / "output"


def main():
    parser = argparse.ArgumentParser(description="Find genes with no grounded QA pairs")
    parser.add_argument("--batch", choices=["batch1", "batch2"],
                        help="Batch name (sets paths automatically)")
    parser.add_argument("--articles", type=str, default=None,
                        help="Articles JSON file (to get gene list)")
    parser.add_argument("--phase3", type=str, default=None,
                        help="Phase 3 validated JSONL file")
    parser.add_argument("--output", type=str, default=None,
                        help="Output targets JSONL for re-fetch")
    parser.add_argument("--include-ungrounded", action="store_true",
                        help="Also re-fetch genes that only produced ungrounded QAs")
    args = parser.parse_args()

    # Set paths from batch name
    if args.batch:
        if args.batch == "batch1":
            articles_path = OUTPUT_DIR / "expansion_articles.json"
        else:
            articles_path = OUTPUT_DIR / "batch2_articles.json"
        phase3_path = OUTPUT_DIR / args.batch / "phase3" / "phase3_validated.jsonl"
        output_path = OUTPUT_DIR / f"{args.batch}_failed_genes.jsonl"
    else:
        articles_path = Path(args.articles)
        phase3_path = Path(args.phase3)
        output_path = Path(args.output)

    # Load all target genes from articles
    print(f"Loading articles from {articles_path}...")
    with open(articles_path) as f:
        articles = json.load(f)

    all_genes = set()
    gene_pmids = defaultdict(set)
    for pmid, info in articles.items():
        gene = info.get("gene", "")
        if gene:
            all_genes.add(gene)
            gene_pmids[gene].add(pmid)

    print(f"Total genes in articles: {len(all_genes)}")

    # Load Phase 3 results and find genes with grounded QAs
    grounded_genes = set()
    ungrounded_only_genes = set()
    gene_stats = defaultdict(lambda: {"cross-grounded": 0, "grounded": 0, "ungrounded": 0})

    if phase3_path.exists():
        print(f"Loading Phase 3 results from {phase3_path}...")
        pmid_to_gene = {pmid: info.get("gene", "") for pmid, info in articles.items()}

        with open(phase3_path) as f:
            for line in f:
                if not line.strip():
                    continue
                record = json.loads(line)
                pmid = record.get("pmid", "")
                grounding = record.get("grounding", "ungrounded")
                gene = pmid_to_gene.get(pmid, "")
                if gene:
                    gene_stats[gene][grounding] += 1

        for gene, stats in gene_stats.items():
            if stats["cross-grounded"] > 0 or stats["grounded"] > 0:
                grounded_genes.add(gene)
            elif stats["ungrounded"] > 0:
                ungrounded_only_genes.add(gene)
    else:
        print(f"WARNING: Phase 3 output not found at {phase3_path}")
        print("Assuming all genes failed (Phase 3 hasn't run yet)")

    # Determine failed genes
    if args.include_ungrounded:
        failed_genes = all_genes - grounded_genes
    else:
        failed_genes = all_genes - grounded_genes - ungrounded_only_genes

    no_output_genes = all_genes - grounded_genes - ungrounded_only_genes

    # Summary
    print(f"\n{'='*60}")
    print(f"GENE OUTCOME SUMMARY")
    print(f"{'='*60}")
    print(f"Total genes in articles:         {len(all_genes):,}")
    print(f"Genes with grounded QAs:         {len(grounded_genes):,}")
    print(f"Genes with ungrounded QAs only:  {len(ungrounded_only_genes):,}")
    print(f"Genes with no pipeline output:   {len(no_output_genes):,}")
    print(f"Failed genes for re-fetch:       {len(failed_genes):,}")

    # Write failed genes as targets
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        for gene in sorted(failed_genes):
            f.write(json.dumps({"gene": gene}) + "\n")

    print(f"\nSaved to {output_path}")

    # Also print some examples
    if grounded_genes:
        print(f"\nSample grounded genes: {sorted(list(grounded_genes))[:10]}")
    if failed_genes:
        print(f"Sample failed genes:   {sorted(list(failed_genes))[:10]}")


if __name__ == "__main__":
    main()
