#!/usr/bin/env python3
"""
Correct gene assignments across the MutQA dataset (formerly fix_gene_assignments.py).

Uses a 3-tier confidence strategy:
  GOLD: PMID→ClinVar resource mutation match (gene + mutation verified externally)
  HIGH: Both recovered_gene and anchor_gene agree, OR single gene in QA text
  LOW:  Best guess from available signals (resource gene set, valid NCBI symbol, fallback)

Reads:  m2tqa_grounded_final.jsonl (416K records; master_grounding_pipeline.py output)
Writes: m2tqa_gene_corrected.jsonl (same records with corrected_gene + gene_confidence fields)

Dependencies (must be supplied; no network calls):
  - --pmid-resource: PMID→ClinVar resource. The MIRIAD-derived PubMed corpus JSON
    {pmid: {"proteins": ["NM_014855.3(AP5Z1):c.1033C>T (p.Arg345Ter)", ...], "article": ...}},
    i.e. the upstream corpus BEFORE phase_0/batch_extractor.py overwrites "proteins" with
    regex-extracted variants (do NOT pass database/articles_with_mutations.json: its
    "proteins" field holds regex-extracted strings, not ClinVar HGVS, so it will not
    reproduce the original GOLD tier).
    Only the "proteins" lists are used. The original run used a June 2025 copy of
    pubmed_to_proteins_with_articles.json; the Dec 2025 *.backfilled.json version of the
    same corpus has the same format but has not been checked for identical GOLD output.
  - --gene-info: NCBI mappings/Homo_sapiens.gene_info
"""

import argparse
import json
import re
import sys
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "claim_extractor"))
import config  # noqa: E402

# ============================================================
# Paths
# ============================================================
INPUT_FILE = config.PIPELINE_ROOT / "uniprot" / "output" / "m2tqa_grounded_final.jsonl"
OUTPUT_FILE = config.PIPELINE_ROOT / "uniprot" / "output" / "m2tqa_gene_corrected.jsonl"
PMID_RESOURCE = config.PIPELINE_ROOT / "database" / "pubmed_to_proteins_with_articles.json"
GENE_INFO = config.PIPELINE_ROOT / "mappings" / "Homo_sapiens.gene_info"

# ============================================================
# Amino acid conversion
# ============================================================
AA_3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': '*',
}

def protein_3to1(p_str):
    def replace_aa(m):
        return AA_3TO1.get(m.group(0), m.group(0))
    return re.sub(r'[A-Z][a-z]{2}', replace_aa, p_str)


def mutation_forms(m):
    """Generate all plausible notation forms for a mutation string."""
    forms = {m}
    if 'del' in m:
        forms.add(re.sub(r'del[ACGT]+$', 'del', m))
    if 'dup' in m:
        forms.add(re.sub(r'dup[ACGT]+$', 'dup', m))
    if m.startswith('p.'):
        forms.add(protein_3to1(m))
    if not m.startswith(('c.', 'p.', 'g.', 'rs')):
        forms.add('c.' + m)
        forms.add('p.' + m)
    return forms


# ============================================================
# Load resources
# ============================================================

def load_gene_info():
    """Load NCBI gene_info → alias_map and valid_genes."""
    alias_map = {}
    with open(GENE_INFO) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 4:
                symbol = parts[2]
                alias_map[symbol] = symbol
                if parts[4] != '-':
                    for syn in parts[4].split('|'):
                        alias_map[syn.strip()] = symbol
    return alias_map


def load_pmid_resource():
    """Load PMID→ClinVar resource, build mutation index and gene sets."""
    with open(PMID_RESOURCE) as f:
        data = json.load(f)

    pmid_index = {}   # pmid -> {mutation_form -> (gene, clinvar_name)}
    pmid_genes = {}   # pmid -> set of genes

    for pmid, info in data.items():
        variants = {}
        genes = set()
        for pstr in info.get('proteins', []):
            gene_match = re.search(r'\((\w+)\)', pstr)
            if not gene_match:
                continue
            gene = gene_match.group(1)
            genes.add(gene)

            coding = re.search(r':(c\.[^\s\)]+)', pstr)
            protein = re.search(r'\((p\.[A-Za-z0-9*=_]+)\)\s*$', pstr)

            if coding:
                variants[coding.group(1)] = (gene, pstr)
            if protein:
                p3 = protein.group(1)
                p1 = protein_3to1(p3)
                variants[p3] = (gene, pstr)
                if p1 != p3:
                    variants[p1] = (gene, pstr)

        if variants:
            pmid_index[pmid] = variants
        pmid_genes[pmid] = genes

    return pmid_index, pmid_genes


# ============================================================
# QA text gene extraction
# ============================================================

# Tokens that match gene symbol patterns but aren't genes in this context
_QA_GENE_BLACKLIST = {
    'SIFT', 'NOT', 'YES', 'THE', 'FOR', 'AND', 'BUT', 'ALL', 'HAS', 'WAS',
    'ARE', 'CAN', 'HAD', 'HER', 'HIS', 'ITS', 'LET', 'MAY', 'NEW', 'OLD',
    'OUR', 'OWN', 'SAY', 'SHE', 'TOO', 'USE', 'DNA', 'RNA', 'PCR', 'NGS',
    'SNV', 'CNV', 'LOH', 'VUS', 'MAF', 'VAF', 'SNP', 'BMI', 'WHO', 'FDA',
    'USA', 'EUR', 'ASN', 'AFR', 'SAS', 'EAS', 'AMR', 'GnomAD',
}

_canonical_genes = None


def load_canonical_genes(alias_map):
    """Extract set of canonical gene symbols from alias_map."""
    global _canonical_genes
    if _canonical_genes is None:
        _canonical_genes = set(alias_map.values())
    return _canonical_genes


def extract_gene_from_qa(record, alias_map):
    """
    Extract gene symbol from QA text fields.

    Returns canonical gene symbol if exactly one valid gene is found,
    otherwise returns None.
    """
    canonical = load_canonical_genes(alias_map)

    text_parts = [
        record.get('question', ''),
        record.get('phase1_answer', ''),
        record.get('phase1_rationale', ''),
        record.get('phase2_answer', ''),
        record.get('phase2_rationale', ''),
    ]
    combined = ' '.join(text_parts)

    # Find uppercase tokens that look like gene symbols
    tokens = re.findall(r'\b([A-Z][A-Z0-9]{1,14})\b', combined)
    found = set()
    for t in tokens:
        if t in _QA_GENE_BLACKLIST:
            continue
        if t in canonical:
            found.add(t)
        elif t in alias_map:
            found.add(alias_map[t])

    if len(found) == 1:
        return found.pop()
    return None


# ============================================================
# Gene correction logic
# ============================================================

def correct_gene(record, alias_map, pmid_index, pmid_genes):
    """
    Determine the best gene assignment for a record.

    Returns (gene, confidence, method, clinvar_name).
    confidence: 'gold', 'high', or 'low'
    method: specific method used within the confidence tier
    """
    pmid = str(record.get('pmid', ''))
    mutation = record.get('mutation', '').strip()
    rg = record.get('recovered_gene', '')
    ag = record.get('anchor_gene_name', '')

    # === GOLD: PMID resource mutation match ===
    clinvar_name = None
    if pmid in pmid_index:
        forms = mutation_forms(mutation)
        for f in forms:
            if f in pmid_index[pmid]:
                gene, clinvar_name = pmid_index[pmid][f]
                return gene, 'gold', 'pmid_mutation_match', clinvar_name

    # === HIGH: Both agree + valid NCBI ===
    if rg and ag and rg == ag and rg in alias_map:
        return alias_map[rg], 'high', 'both_agree', None

    # === HIGH: Single gene mentioned in QA text ===
    qa_gene = extract_gene_from_qa(record, alias_map)
    if qa_gene:
        return qa_gene, 'high', 'qa_text', None

    # === LOW: Best guess from remaining signals ===

    # Try PMID resource gene set (pick whichever candidate appears)
    if rg and ag and rg != ag and pmid in pmid_genes:
        resource_genes = pmid_genes[pmid]
        rg_canonical = alias_map.get(rg, rg)
        ag_canonical = alias_map.get(ag, ag)
        rg_in = rg_canonical in resource_genes or rg in resource_genes
        ag_in = ag_canonical in resource_genes or ag in resource_genes

        if rg_in and not ag_in:
            best = rg_canonical if rg_canonical in resource_genes else rg
            return best, 'low', 'resource_rg', None
        elif ag_in and not rg_in:
            best = ag_canonical if ag_canonical in resource_genes else ag
            return best, 'low', 'resource_ag', None
        elif rg_in and ag_in:
            return rg_canonical, 'low', 'resource_both', None

    # recovered_gene is valid NCBI symbol
    if rg and rg in alias_map:
        return alias_map[rg], 'low', 'rg_valid', None

    # anchor_gene is valid NCBI symbol
    if ag and ag in alias_map:
        return alias_map[ag], 'low', 'ag_valid', None

    # Fallback
    return rg or ag or '', 'low', 'fallback', None


# ============================================================
# Main
# ============================================================

def main():
    global INPUT_FILE, OUTPUT_FILE, PMID_RESOURCE, GENE_INFO
    parser = argparse.ArgumentParser(description="Assign corrected_gene + gene_confidence (gold/high/low)")
    parser.add_argument("--input", type=Path, default=INPUT_FILE,
                        help="m2tqa_grounded_final.jsonl (master_grounding_pipeline.py output)")
    parser.add_argument("--output", type=Path, default=OUTPUT_FILE)
    parser.add_argument("--pmid-resource", type=Path, default=PMID_RESOURCE,
                        help="PMID->ClinVar corpus JSON with HGVS 'proteins' lists (see module docstring)")
    parser.add_argument("--gene-info", type=Path, default=GENE_INFO, help="NCBI Homo_sapiens.gene_info")
    args = parser.parse_args()
    INPUT_FILE, OUTPUT_FILE = args.input, args.output
    PMID_RESOURCE, GENE_INFO = args.pmid_resource, args.gene_info

    print("Loading NCBI gene info...")
    alias_map = load_gene_info()
    print(f"  {len(alias_map):,} symbols/aliases")

    print("Loading PMID→ClinVar resource...")
    pmid_index, pmid_genes = load_pmid_resource()
    print(f"  {len(pmid_index):,} PMIDs with mutation index")

    print(f"Processing {INPUT_FILE}...")
    confidence_counts = Counter()
    method_counts = Counter()
    corrections = 0
    total = 0

    with open(INPUT_FILE) as fin, open(OUTPUT_FILE, 'w') as fout:
        for line in fin:
            record = json.loads(line)
            total += 1

            old_gene = record.get('gene', '')
            best_gene, confidence, method, clinvar_name = correct_gene(
                record, alias_map, pmid_index, pmid_genes
            )

            record['corrected_gene'] = best_gene
            record['gene_confidence'] = confidence
            record['gene_correction_method'] = method
            if clinvar_name:
                record['clinvar_variant_name'] = clinvar_name

            if best_gene != old_gene:
                corrections += 1

            confidence_counts[confidence] += 1
            method_counts[method] += 1
            fout.write(json.dumps(record) + '\n')

            if total % 100000 == 0:
                print(f"  {total:,} processed...", flush=True)

    print(f"\nDone: {total:,} records")
    print(f"Corrections: {corrections:,} ({100*corrections/total:.1f}%)")

    print(f"\nConfidence summary:")
    for conf in ['gold', 'high', 'low']:
        c = confidence_counts[conf]
        print(f"  {conf.upper():6s}: {c:8,d} ({100*c/total:.1f}%)")

    print(f"\nMethod breakdown:")
    for method, c in sorted(method_counts.items(), key=lambda x: -x[1]):
        print(f"  {method:25s}: {c:8,d} ({100*c/total:.1f}%)")

    print(f"\nOutput: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
