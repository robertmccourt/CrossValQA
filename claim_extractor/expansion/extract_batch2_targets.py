#!/usr/bin/env python3
"""
Extract Batch 2 target genes for pipeline expansion.

Batch 2 = all protein-coding genes NOT in current dataset AND NOT in Batch 1.
These genes have no ClinVar-cited PMIDs, so we search PMC directly.

Input:  uniprot/output/m2tqa_gene_corrected.jsonl, mappings/Homo_sapiens.gene_info,
        expansion/output/expansion_articles.json (Batch 1 articles)
Output: expansion/output/batch2_targets.jsonl ({"gene": ...} per line)

Usage:
    python extract_batch2_targets.py
"""

import argparse
import json
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import UNIPROT_GENE_CORRECTED_FILE, GENE_INFO_FILE

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
BATCH1_ARTICLES = OUTPUT_DIR / "expansion_articles.json"


def load_current_genes() -> set:
    """Load gene symbols from the current dataset."""
    genes = set()
    with open(UNIPROT_GENE_CORRECTED_FILE) as f:
        for line in f:
            r = json.loads(line)
            g = r.get('corrected_gene') or r.get('recovered_gene') or r.get('anchor_gene_name', '')
            if g:
                genes.add(g)
    return genes


def load_batch1_genes() -> set:
    """Load genes already fetched in Batch 1."""
    genes = set()
    if BATCH1_ARTICLES.exists():
        with open(BATCH1_ARTICLES) as f:
            data = json.load(f)
        for entry in data.values():
            if "gene" in entry:
                genes.add(entry["gene"])
    return genes


def load_protein_coding_genes() -> set:
    """Load protein-coding gene symbols from NCBI."""
    genes = set()
    with open(GENE_INFO_FILE) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 9 and parts[9] == 'protein-coding':
                symbol = parts[2]
                if symbol and ';' not in symbol and ' ' not in symbol and len(symbol) <= 15:
                    genes.add(symbol)
    return genes


def main():
    argparse.ArgumentParser(
        description="Batch 2 targets: protein-coding genes not in the dataset and not in Batch 1."
    ).parse_args()

    print("Loading current dataset genes...")
    current_genes = load_current_genes()
    print(f"  Current dataset: {len(current_genes):,} genes")

    print("Loading Batch 1 genes...")
    batch1_genes = load_batch1_genes()
    print(f"  Batch 1: {len(batch1_genes):,} genes")

    print("Loading protein-coding genes from NCBI...")
    coding_genes = load_protein_coding_genes()
    print(f"  Protein-coding: {len(coding_genes):,} genes")

    # Batch 2 = protein-coding minus current minus batch1
    batch2_genes = sorted(coding_genes - current_genes - batch1_genes)
    print(f"\n  Batch 2 targets: {len(batch2_genes):,} genes")

    # Write targets (same format as Batch 1 for fetcher compatibility)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file = OUTPUT_DIR / "batch2_targets.jsonl"
    with open(output_file, 'w') as f:
        for gene in batch2_genes:
            f.write(json.dumps({"gene": gene}) + "\n")

    print(f"  Output: {output_file}")
    print(f"  Records: {len(batch2_genes):,}")


if __name__ == "__main__":
    main()
