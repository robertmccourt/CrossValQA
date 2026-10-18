#!/usr/bin/env python3
"""
Generate statistics for the claim extractor dataset.
Reads the Phase 4 gene-recovery output (m2tqa_gene_recovered.jsonl, which carries
the recovered_gene / recovery_tier fields) and writes stats.md.

Usage:
    python generate_stats.py [--input PATH] [--output PATH]
"""

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import PHASE4_GENE_RECOVERED_FILE, CLAIM_EXTRACTOR_DIR

# Configuration (overridable via --input / --output)
DATASET_PATH = PHASE4_GENE_RECOVERED_FILE
OUTPUT_PATH = CLAIM_EXTRACTOR_DIR / "stats.md"


def classify_mutation(mutation: str) -> str:
    """Classify mutation into type based on nomenclature."""
    mutation = mutation.strip()

    # HGVS coding DNA (c.)
    if re.match(r'^c\.', mutation, re.IGNORECASE):
        return "Coding DNA (c.)"

    # HGVS protein (p.)
    if re.match(r'^p\.', mutation, re.IGNORECASE):
        return "Protein (p.)"

    # HGVS genomic (g.)
    if re.match(r'^g\.', mutation, re.IGNORECASE):
        return "Genomic (g.)"

    # HGVS mitochondrial (m.)
    if re.match(r'^m\.', mutation, re.IGNORECASE):
        return "Mitochondrial (m.)"

    # HGVS RNA (r.)
    if re.match(r'^r\.', mutation, re.IGNORECASE):
        return "RNA (r.)"

    # rsID (SNP database ID)
    if re.match(r'^rs\d+', mutation, re.IGNORECASE):
        return "rsID (dbSNP)"

    # Amino acid substitution (e.g., V600E, R248W)
    if re.match(r'^[A-Z]\d+[A-Z]$', mutation):
        return "AA Substitution (e.g., V600E)"

    # Amino acid with fs/del/ins (e.g., L858R, E746_A750del)
    if re.match(r'^[A-Z]\d+', mutation) and any(x in mutation.lower() for x in ['del', 'ins', 'dup', 'fs']):
        return "AA Indel/Frameshift"

    # Nucleotide substitution (e.g., A>G, 1234A>G)
    if re.search(r'\d*[ACGT]>[ACGT]', mutation, re.IGNORECASE):
        return "Nucleotide Substitution"

    # IVS notation (intronic)
    if re.match(r'^IVS', mutation, re.IGNORECASE):
        return "Intronic (IVS)"

    # Exon-level (e.g., exon 19 deletion)
    if 'exon' in mutation.lower():
        return "Exon-level"

    # Simple nucleotide notation (e.g., T89C, G12D without p. prefix)
    if re.match(r'^[ACGT]\d+[ACGT]$', mutation, re.IGNORECASE):
        return "Simple Nucleotide"

    # Deletion patterns
    if 'del' in mutation.lower():
        return "Deletion"

    # Insertion patterns
    if 'ins' in mutation.lower():
        return "Insertion"

    # Duplication patterns
    if 'dup' in mutation.lower():
        return "Duplication"

    # Complex/other
    return "Other"


def classify_variant_type(mutation: str) -> str:
    """Classify mutation by biological variant type (SNV, indel, etc.)."""
    mutation_lower = mutation.lower().strip()
    mutation_upper = mutation.upper().strip()

    # Frameshift (check first as it's specific)
    if 'fs' in mutation_lower or 'frameshift' in mutation_lower:
        return "Frameshift"

    # Nonsense/Stop-gain mutations
    if any(x in mutation for x in ['*', 'Ter', 'X']) and not mutation_lower.startswith('exon'):
        # Check for patterns like p.Arg100*, p.R100X, p.Arg100Ter
        if re.search(r'(p\.)?[A-Z][a-z]{0,2}\d+(\*|Ter|X)($|[^a-zA-Z])', mutation):
            return "Nonsense (Stop-gain)"
        if re.search(r'\d+[X\*]$', mutation_upper):
            return "Nonsense (Stop-gain)"

    # Splice site variants
    if any(x in mutation_lower for x in ['splice', 'ivs', '+1', '+2', '-1', '-2']):
        return "Splice Site"
    # HGVS splice notation: c.123+1G>A or c.123-2A>G
    if re.search(r'c\.\d+[+-]\d+', mutation_lower):
        return "Splice Site"

    # Complex indels (delins)
    if 'delins' in mutation_lower or 'indel' in mutation_lower:
        return "Complex Indel (delins)"

    # Deletions
    if 'del' in mutation_lower:
        # Check if it's a large deletion (exon-level or multi-codon)
        if 'exon' in mutation_lower:
            return "Large Deletion (exon-level)"
        if re.search(r'_[A-Z]?\d+', mutation):  # Range like E746_A750del
            return "In-frame Deletion"
        return "Deletion"

    # Insertions
    if 'ins' in mutation_lower:
        return "Insertion"

    # Duplications
    if 'dup' in mutation_lower:
        return "Duplication"

    # SNV patterns - single nucleotide changes
    # Pattern: c.123A>G, 123A>G, A>G
    if re.search(r'[ACGT]>[ACGT]', mutation_upper):
        return "SNV (Single Nucleotide Variant)"

    # Protein missense: p.Val600Glu, p.V600E, V600E
    # Single amino acid substitution (not del/ins/fs)
    if re.match(r'^(p\.)?[A-Z][a-z]{0,2}\d+[A-Z][a-z]{0,2}$', mutation):
        return "Missense (Amino Acid Substitution)"
    if re.match(r'^[A-Z]\d+[A-Z]$', mutation_upper):
        return "Missense (Amino Acid Substitution)"

    # Simple nucleotide notation that implies SNV: T89C, G12D
    if re.match(r'^[ACGT]\d+[ACGT]$', mutation_upper):
        return "SNV (Single Nucleotide Variant)"

    # rsID - these are typically SNVs but we can't be certain
    if re.match(r'^rs\d+', mutation_lower):
        return "dbSNP Variant (type unknown)"

    # Synonymous/silent mutations
    if '=' in mutation or 'syn' in mutation_lower:
        return "Synonymous (Silent)"

    # Check for missense-like patterns in HGVS protein notation
    if mutation_lower.startswith('p.'):
        # p.Ala123Gly pattern (three-letter codes)
        if re.match(r'^p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}$', mutation):
            return "Missense (Amino Acid Substitution)"

    return "Other/Unclassified"


def main():
    global DATASET_PATH, OUTPUT_PATH
    parser = argparse.ArgumentParser(description="Generate dataset statistics (stats.md)")
    parser.add_argument("--input", default=str(DATASET_PATH), help="Phase 4 gene-recovered JSONL")
    parser.add_argument("--output", default=str(OUTPUT_PATH), help="Output markdown path")
    args = parser.parse_args()
    DATASET_PATH = Path(args.input)
    OUTPUT_PATH = Path(args.output)

    print(f"Loading dataset from {DATASET_PATH}...")

    # Counters and collectors
    question_types = Counter()
    mutation_types = Counter()  # Nomenclature-based (c., p., etc.)
    variant_types = Counter()   # Biological type (SNV, indel, etc.)
    grounding_types = Counter()
    confidence_levels = Counter()
    recovery_tiers = Counter()

    genes = Counter()
    recovered_genes = Counter()
    final_genes = Counter()  # Best available gene

    pmids = set()
    unique_mutations = set()
    pmid_mutation_pairs = set()

    # Per-gene mutation counts
    gene_mutations = defaultdict(set)

    total_records = 0
    records_with_gene = 0
    records_with_recovered_gene = 0

    # Process dataset
    with open(DATASET_PATH) as f:
        for line in f:
            record = json.loads(line)
            total_records += 1

            # Question type
            qtype = record.get("question_type", "unknown")
            question_types[qtype] += 1

            # Mutation
            mutation = record.get("mutation", "")
            unique_mutations.add(mutation)
            mutation_types[classify_mutation(mutation)] += 1
            variant_types[classify_variant_type(mutation)] += 1

            # PMID
            pmid = record.get("pmid", "")
            pmids.add(pmid)
            pmid_mutation_pairs.add((pmid, mutation))

            # Grounding
            grounding = record.get("grounding", "unknown")
            grounding_types[grounding] += 1

            # Confidence
            confidence = record.get("phase2_confidence", "unknown")
            confidence_levels[confidence] += 1

            # Gene handling
            gene = record.get("gene")
            recovered_gene = record.get("recovered_gene")
            recovery_tier = record.get("recovery_tier")

            if gene:
                genes[gene] += 1
                records_with_gene += 1

            if recovered_gene:
                recovered_genes[recovered_gene] += 1
                records_with_recovered_gene += 1

            if recovery_tier:
                recovery_tiers[recovery_tier] += 1

            # Determine final/best gene
            final_gene = recovered_gene or gene
            if final_gene:
                final_genes[final_gene] += 1
                gene_mutations[final_gene].add(mutation)

    print(f"Processed {total_records:,} records")

    # Generate markdown report
    md = []
    md.append("# Claim Extractor Dataset Statistics\n")
    md.append(f"**Dataset:** `{DATASET_PATH.name}`\n")
    md.append(f"**Generated:** Auto-generated statistics\n")
    md.append("")

    # Overview
    md.append("## Overview\n")
    md.append("| Metric | Value |")
    md.append("|--------|-------|")
    md.append(f"| Total QA pairs | {total_records:,} |")
    md.append(f"| Unique mutations | {len(unique_mutations):,} |")
    md.append(f"| Unique PMIDs (articles) | {len(pmids):,} |")
    md.append(f"| Unique PMID-mutation pairs | {len(pmid_mutation_pairs):,} |")
    md.append(f"| Unique genes (final) | {len(final_genes):,} |")
    md.append(f"| Avg QA pairs per mutation | {total_records / len(unique_mutations):.2f} |")
    md.append(f"| Avg mutations per article | {len(pmid_mutation_pairs) / len(pmids):.2f} |")
    md.append("")

    # Question Type Distribution
    md.append("## Question Type Distribution\n")
    md.append("| Question Type | Count | Percentage |")
    md.append("|---------------|-------|------------|")
    for qtype, count in question_types.most_common():
        pct = count / total_records * 100
        md.append(f"| {qtype} | {count:,} | {pct:.1f}% |")
    md.append("")

    # Mutation Nomenclature Distribution
    md.append("## Mutation Nomenclature Distribution\n")
    md.append("| Nomenclature Type | Count | Percentage |")
    md.append("|-------------------|-------|------------|")
    for mtype, count in mutation_types.most_common():
        pct = count / total_records * 100
        md.append(f"| {mtype} | {count:,} | {pct:.1f}% |")
    md.append("")

    # Biological Variant Type Distribution
    md.append("## Biological Variant Type Distribution\n")
    md.append("| Variant Type | Count | Percentage |")
    md.append("|--------------|-------|------------|")
    for vtype, count in variant_types.most_common():
        pct = count / total_records * 100
        md.append(f"| {vtype} | {count:,} | {pct:.1f}% |")
    md.append("")

    # Grounding Distribution
    md.append("## Grounding Quality Distribution\n")
    md.append("| Grounding Level | Count | Percentage |")
    md.append("|-----------------|-------|------------|")
    for gtype, count in grounding_types.most_common():
        pct = count / total_records * 100
        md.append(f"| {gtype} | {count:,} | {pct:.1f}% |")
    md.append("")

    # Confidence Distribution
    md.append("## Phase 2 Confidence Distribution\n")
    md.append("| Confidence | Count | Percentage |")
    md.append("|------------|-------|------------|")
    for conf, count in confidence_levels.most_common():
        pct = count / total_records * 100
        md.append(f"| {conf} | {count:,} | {pct:.1f}% |")
    md.append("")

    # Gene Coverage
    md.append("## Gene Coverage\n")
    md.append("| Metric | Value |")
    md.append("|--------|-------|")
    md.append(f"| Records with original gene | {records_with_gene:,} ({records_with_gene/total_records*100:.1f}%) |")
    md.append(f"| Records with recovered gene | {records_with_recovered_gene:,} ({records_with_recovered_gene/total_records*100:.1f}%) |")
    md.append(f"| Total unique genes (original) | {len(genes):,} |")
    md.append(f"| Total unique genes (recovered) | {len(recovered_genes):,} |")
    md.append(f"| Total unique genes (final/best) | {len(final_genes):,} |")
    md.append("")

    # Gene Recovery Tiers
    if recovery_tiers:
        md.append("### Gene Recovery Tiers\n")
        md.append("| Recovery Tier | Count | Percentage |")
        md.append("|---------------|-------|------------|")
        for tier, count in recovery_tiers.most_common():
            pct = count / total_records * 100
            md.append(f"| {tier} | {count:,} | {pct:.1f}% |")
        md.append("")

    # Known non-gene terms that may appear in gene field
    NON_GENES = {"VUS", "SIFT", "ACMG", "CADD", "REVEL", "PolyPhen", "GERP", "PhyloP",
                 "gnomAD", "ExAC", "ClinVar", "HGMD", "dbSNP", "COSMIC", "NA", "N/A",
                 "Unknown", "unknown", "NONE", "None", "null"}

    # Top Genes (filtered)
    md.append("### Top 30 Genes by QA Pair Count\n")
    md.append("> **Note:** Some entries like 'VUS' (Variant of Uncertain Significance) and 'SIFT' (prediction tool) ")
    md.append("> are misclassified as genes in the original data. These are marked with * below.\n")
    md.append("")
    md.append("| Rank | Gene | QA Pairs | Unique Mutations |")
    md.append("|------|------|----------|------------------|")
    rank = 0
    for gene, count in final_genes.most_common(50):
        if rank >= 30:
            break
        n_mutations = len(gene_mutations[gene])
        marker = "*" if gene in NON_GENES else ""
        md.append(f"| {rank + 1} | {gene}{marker} | {count:,} | {n_mutations:,} |")
        rank += 1
    md.append("")

    # Count non-genes
    non_gene_count = sum(final_genes[g] for g in NON_GENES if g in final_genes)
    actual_genes = {g: c for g, c in final_genes.items() if g not in NON_GENES}
    md.append(f"**Data quality note:** {non_gene_count:,} records ({non_gene_count/total_records*100:.1f}%) have ")
    md.append(f"non-gene identifiers in the gene field. There are **{len(actual_genes):,}** verified gene symbols.\n")
    md.append("")

    # Sample mutations by nomenclature type
    md.append("## Sample Mutations by Nomenclature\n")

    # Collect unique samples for both classification types
    nomenclature_samples = defaultdict(set)
    variant_samples = defaultdict(set)
    with open(DATASET_PATH) as f:
        for line in f:
            record = json.loads(line)
            mutation = record.get("mutation", "")

            mtype = classify_mutation(mutation)
            if len(nomenclature_samples[mtype]) < 5:
                nomenclature_samples[mtype].add(mutation)

            vtype = classify_variant_type(mutation)
            if len(variant_samples[vtype]) < 5:
                variant_samples[vtype].add(mutation)

    for mtype in mutation_types.keys():
        samples = list(nomenclature_samples.get(mtype, []))
        if samples:
            md.append(f"**{mtype}:** `{'`, `'.join(samples[:5])}`\n")
    md.append("")

    # Sample mutations by biological variant type
    md.append("## Sample Mutations by Variant Type\n")
    for vtype in variant_types.keys():
        samples = list(variant_samples.get(vtype, []))
        if samples:
            md.append(f"**{vtype}:** `{'`, `'.join(samples[:5])}`\n")
    md.append("")

    # Summary statistics
    md.append("## Summary\n")
    md.append(f"This dataset contains **{total_records:,}** validated QA pairs covering "
              f"**{len(unique_mutations):,}** unique genetic mutations across "
              f"**{len(final_genes):,}** genes from **{len(pmids):,}** PubMed articles.\n")
    md.append("")

    cross_grounded = grounding_types.get("cross-grounded", 0)
    grounded = grounding_types.get("grounded", 0)
    md.append(f"- **{cross_grounded/total_records*100:.1f}%** of QA pairs are cross-grounded (highest quality)\n")
    md.append(f"- **{(cross_grounded + grounded)/total_records*100:.1f}%** of QA pairs pass validation\n")

    high_conf = confidence_levels.get("high", 0)
    md.append(f"- **{high_conf/total_records*100:.1f}%** have high confidence from Phase 2 LLM\n")

    # Write output
    output_text = "\n".join(md)
    with open(OUTPUT_PATH, "w") as f:
        f.write(output_text)

    print(f"\nStatistics written to {OUTPUT_PATH}")
    print("\n" + "=" * 60)
    print(output_text)


if __name__ == "__main__":
    main()
