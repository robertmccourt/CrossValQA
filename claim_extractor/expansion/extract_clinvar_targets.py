#!/usr/bin/env python3
"""
Batch 1 targets: pick one ClinVar-cited PMID per protein-coding gene that is
missing from the current MutQA dataset.

Selects NCBI protein-coding genes not in the pre-expansion dataset, maps them
to PMIDs cited by ClinVar for SNVs / small indels, drops PMIDs already in the
Phase 0 article database, ranks genes by number of available PMIDs, and keeps
the smallest new PMID per gene.

Input:  uniprot/output/m2tqa_gene_corrected.jsonl (pre-expansion dataset)
        mappings/Homo_sapiens.gene_info (NCBI)
        expansion/data/variant_summary.txt.gz, expansion/data/var_citations.txt (ClinVar FTP)
        database/articles_with_mutations.json (existing Phase 0 articles)
Output: expansion/output/expansion_targets.jsonl, expansion/output/expansion_pmids.txt

Usage:
    python extract_clinvar_targets.py
    python extract_clinvar_targets.py --max-genes 8000 --variant-summary PATH --var-citations PATH
"""

import json
import gzip
import argparse
from collections import defaultdict
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import (
    UNIPROT_GENE_CORRECTED_FILE, GENE_INFO_FILE, ARTICLES_FILE,
)

CLINVAR_DIR = Path(__file__).resolve().parent / "data"
VARIANT_SUMMARY = CLINVAR_DIR / "variant_summary.txt.gz"
VAR_CITATIONS = CLINVAR_DIR / "var_citations.txt"
OUTPUT_DIR = Path(__file__).resolve().parent / "output"


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


def load_protein_coding_genes() -> set:
    """Load protein-coding gene symbols from NCBI."""
    genes = set()
    with open(GENE_INFO_FILE) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 9 and parts[9] == 'protein-coding':
                symbol = parts[2]
                # Filter: single clean gene symbol
                if symbol and ';' not in symbol and ' ' not in symbol and len(symbol) <= 15:
                    genes.add(symbol)
    return genes


def build_allele_to_gene_variant() -> dict:
    """Build AlleleID -> (gene, variant_name) mapping from ClinVar variant_summary."""
    mapping = {}
    with gzip.open(VARIANT_SUMMARY, 'rt') as f:
        header = f.readline().strip().split('\t')
        allele_idx = header.index('#AlleleID')
        gene_idx = header.index('GeneSymbol')
        name_idx = header.index('Name')
        type_idx = header.index('Type')

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) <= max(allele_idx, gene_idx, name_idx, type_idx):
                continue

            allele_id = parts[allele_idx]
            gene = parts[gene_idx]
            name = parts[name_idx]
            vtype = parts[type_idx]

            # Only single nucleotide variants and small indels
            if vtype not in ('single nucleotide variant', 'Deletion', 'Duplication',
                             'Insertion', 'Indel'):
                continue

            if gene and ';' not in gene and allele_id not in mapping:
                mapping[allele_id] = (gene, name)

    return mapping


def build_gene_pmid_map(allele_map: dict, target_genes: set) -> dict:
    """Build gene -> set of PMIDs from ClinVar citations."""
    gene_pmids = defaultdict(set)

    with open(VAR_CITATIONS) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6 or parts[4] != 'PubMed':
                continue

            allele_id = parts[0]
            pmid = parts[5]

            if allele_id in allele_map:
                gene, variant = allele_map[allele_id]
                if gene in target_genes:
                    gene_pmids[gene].add(pmid)

    return gene_pmids


def load_existing_pmids() -> set:
    """Load PMIDs already in our article database."""
    with open(ARTICLES_FILE) as f:
        data = json.load(f)
    return set(data.keys())


def main():
    global VARIANT_SUMMARY, VAR_CITATIONS
    parser = argparse.ArgumentParser(description="Extract ClinVar targets for expansion")
    parser.add_argument("--max-genes", type=int, default=8000, help="Max genes to target")
    parser.add_argument("--variant-summary", type=Path, default=VARIANT_SUMMARY,
                        help="ClinVar variant_summary.txt.gz (default: expansion/data/)")
    parser.add_argument("--var-citations", type=Path, default=VAR_CITATIONS,
                        help="ClinVar var_citations.txt (default: expansion/data/)")
    args = parser.parse_args()
    VARIANT_SUMMARY, VAR_CITATIONS = args.variant_summary, args.var_citations

    # Step 1: Identify missing genes
    print("Loading current genes...")
    current_genes = load_current_genes()
    print(f"  Current genes: {len(current_genes):,}")

    print("Loading protein-coding genes from NCBI...")
    coding_genes = load_protein_coding_genes()
    print(f"  Protein-coding genes: {len(coding_genes):,}")

    # Step 2: Build ClinVar mappings
    print("Building ClinVar allele -> gene/variant mapping...")
    allele_map = build_allele_to_gene_variant()
    print(f"  Allele mappings: {len(allele_map):,}")

    # Target genes: protein-coding, not in current dataset
    target_genes = coding_genes - current_genes
    print(f"  Target genes (missing + protein-coding): {len(target_genes):,}")

    print("Building gene -> PMID mapping from ClinVar citations...")
    gene_pmids = build_gene_pmid_map(allele_map, target_genes)
    print(f"  Genes with ClinVar-cited PMIDs: {len(gene_pmids):,}")

    # Step 3: Filter to NEW PMIDs only (not already processed)
    print("Loading existing PMIDs...")
    existing_pmids = load_existing_pmids()
    print(f"  Existing PMIDs: {len(existing_pmids):,}")

    genes_with_new_pmids = {}
    for gene, pmids in gene_pmids.items():
        new_pmids = pmids - existing_pmids
        if new_pmids:
            genes_with_new_pmids[gene] = new_pmids

    print(f"  Genes with at least 1 NEW PMID: {len(genes_with_new_pmids):,}")

    # Step 4: Rank genes by number of available PMIDs (more = better studied)
    ranked = sorted(genes_with_new_pmids.items(), key=lambda x: -len(x[1]))
    selected = ranked[:args.max_genes]
    print(f"\n  Selected top {len(selected):,} genes")

    # Step 5: Pick 1 PMID per gene (first available)
    targets = []
    selected_pmids = set()
    for gene, pmids in selected:
        pmid = sorted(pmids)[0]  # deterministic: pick smallest PMID
        targets.append({
            "gene": gene,
            "pmid": pmid,
            "available_pmids": len(pmids),
        })
        selected_pmids.add(pmid)

    print(f"  Unique PMIDs: {len(selected_pmids):,}")

    # Step 6: Save
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file = OUTPUT_DIR / "expansion_targets.jsonl"
    with open(output_file, 'w') as f:
        for t in targets:
            f.write(json.dumps(t) + "\n")
    print(f"\n  Output: {output_file}")
    print(f"  Records: {len(targets):,}")

    # Also save just the PMIDs for the fetcher
    pmid_file = OUTPUT_DIR / "expansion_pmids.txt"
    with open(pmid_file, 'w') as f:
        for pmid in sorted(selected_pmids):
            f.write(pmid + "\n")
    print(f"  PMID list: {pmid_file}")
    print(f"  Unique PMIDs to fetch: {len(selected_pmids):,}")

    # Stats
    pmid_counts = [t["available_pmids"] for t in targets]
    print(f"\n  PMIDs per gene: min={min(pmid_counts)}, median={sorted(pmid_counts)[len(pmid_counts)//2]}, max={max(pmid_counts)}")


if __name__ == "__main__":
    main()
