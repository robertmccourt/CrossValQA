#!/usr/bin/env python3
"""
Normalize question_type to the 8 main categories used in the paper.

The Phase 1 LLM labels question_type freely, so the raw field has ~180 values
(compound labels such as "frequency|other", and long-tail labels such as
"zygosity" or "genotype"). This step maps every record to one of:

    function, pathogenicity, disease, other, frequency, clinical, mechanism, inheritance

Rules (applied in this order):
  1. a value that is already a main category is left unchanged;
  2. a compound value ("a|b", "a/b") maps to its first token that is a main category;
  3. any other value is looked up in LONG_TAIL_MAP (semantic-similarity mapping).
For every changed record the raw label is kept in a new field
`question_type_original`, appended at the end of the record. Values that match no
rule are left unchanged and reported.

Runs on the output of correct_gene_assignments.py (m2tqa_gene_corrected.jsonl),
before the expansion merge (claim_extractor/expansion/merge_expansion.py).
Expansion records are merged in later without this step; their remaining
non-main labels ("functional", "functional effect", ...) are normalized at
release time by release/prepare_release_v2.py.

Provenance: the original normalization was done ad hoc (no script survived).
LONG_TAIL_MAP and the compound rule were reconstructed from the archived
before/after files (m2tqa_gene_corrected_pre_qt_norm.jsonl ->
m2tqa_gene_corrected.jsonl); applying this script to the "before" file reproduces
the "after" file byte for byte (416,188 records, 3,892 changed).

Input:  uniprot/output/m2tqa_gene_corrected_pre_qt_norm.jsonl
        (= correct_gene_assignments.py output, renamed before normalizing)
Output: uniprot/output/m2tqa_gene_corrected.jsonl
"""

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "claim_extractor"))
import config  # noqa: E402

INPUT_FILE = config.UNIPROT_OUTPUT_DIR / "m2tqa_gene_corrected_pre_qt_norm.jsonl"
OUTPUT_FILE = config.UNIPROT_GENE_CORRECTED_FILE

MAIN_TYPES = ["function", "pathogenicity", "disease", "other",
              "frequency", "clinical", "mechanism", "inheritance"]

# Non-compound labels -> main category
LONG_TAIL_MAP = {
    "acmg_codes": "pathogenicity",
    "acmg_criteria": "pathogenicity",
    "association": "disease",
    "availability": "other",
    "biochemical": "function",
    "biological_markers": "other",
    "brain anomalies": "other",
    "brain_anomalies": "other",
    "carrier status": "inheritance",
    "carrier_status": "inheritance",
    "category": "other",
    "causality": "pathogenicity",
    "classification": "pathogenicity",
    "classification_and_association": "other",
    "clinical features": "clinical",
    "clinical_features": "clinical",
    "co-occurrence": "other",
    "compound_heterozygote_status": "other",
    "confirmation": "other",
    "conservation": "function",
    "criteria": "other",
    "demographic": "frequency",
    "demographics": "frequency",
    "detection": "clinical",
    "developmental_delay": "other",
    "diagnosis": "clinical",
    "disease association": "disease",
    "disease_and_clinical_features": "disease",
    "disease_association": "disease",
    "distribution": "other",
    "drug_response": "clinical",
    "epilepsy_classification": "other",
    "evidence": "other",
    "family_history": "inheritance",
    "frequency_and_prediction": "frequency",
    "functional": "function",
    "functional_effect": "function",
    "gene": "other",
    "gene_association": "disease",
    "gene_location": "function",
    "genetic_information": "other",
    "genetics": "other",
    "genotype": "function",
    "genotype-phenotype": "other",
    "geographic_origin": "frequency",
    "germline_testing": "other",
    "hearing_level": "other",
    "heterozygosity": "inheritance",
    "immunological": "other",
    "inclusion_criteria": "other",
    "interpretation": "pathogenicity",
    "known_mutation": "other",
    "location": "function",
    "location_of_effect": "other",
    "methodology": "other",
    "microcephaly": "other",
    "molecular_analysis": "function",
    "molecular_consequence": "function",
    "mutation_description": "other",
    "mutation_type": "function",
    "negative": "other",
    "negative_association": "other",
    "negative_disease_association": "disease",
    "negative_finding": "other",
    "negative_framing": "other",
    "novelty": "other",
    "occurrence": "other",
    "onset": "clinical",
    "origin": "frequency",
    "outcome": "clinical",
    "pathology": "other",
    "pathway": "mechanism",
    "patient_association": "other",
    "patient_identification": "other",
    "patient_phenotype": "clinical",
    "penetrance": "inheritance",
    "pharmacogenetics": "clinical",
    "phenotype": "clinical",
    "phenotype_association": "clinical",
    "phenotypic": "clinical",
    "population": "frequency",
    "population_association": "frequency",
    "prediction": "pathogenicity",
    "prenatal_features": "other",
    "presence": "other",
    "prognosis": "clinical",
    "protein_change": "other",
    "protein_consequence": "function",
    "screening": "clinical",
    "screening_method": "clinical",
    "segregation": "inheritance",
    "severity": "clinical",
    "significance": "pathogenicity",
    "technique": "other",
    "timing": "other",
    "treatment": "clinical",
    "treatment_association": "other",
    "treatment_response": "clinical",
    "variant_identification": "other",
    "variant_type": "function",
    "vital_status": "other",
    "zygosity": "inheritance",
    "zygosity_comparison": "other",
}


def normalize_question_type(qt):
    """Return the main category for a raw question_type, or None if no rule matches."""
    if qt in MAIN_TYPES:
        return qt
    if isinstance(qt, str) and ("|" in qt or "/" in qt):
        for token in re.split(r"[|/]", qt):
            if token in MAIN_TYPES:
                return token
    return LONG_TAIL_MAP.get(qt)


def main():
    parser = argparse.ArgumentParser(description="Normalize question_type to the 8 main categories")
    parser.add_argument("--input", type=Path, default=INPUT_FILE,
                        help="correct_gene_assignments.py output (before normalization)")
    parser.add_argument("--output", type=Path, default=OUTPUT_FILE)
    args = parser.parse_args()

    if args.input.resolve() == args.output.resolve():
        sys.exit("ERROR: --input and --output must differ")

    total = changed = 0
    unmapped = Counter()
    final_counts = Counter()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.input) as f_in, open(args.output, "w") as f_out:
        for line in f_in:
            record = json.loads(line)
            total += 1
            qt = record.get("question_type")
            new_qt = normalize_question_type(qt)
            if new_qt is None:
                unmapped[qt] += 1
            elif new_qt != qt:
                record["question_type"] = new_qt
                record["question_type_original"] = qt
                changed += 1
            final_counts[record.get("question_type")] += 1
            f_out.write(json.dumps(record) + "\n")

    print(f"Records: {total:,}")
    print(f"Changed: {changed:,} ({changed / max(total, 1) * 100:.1f}%)")
    for qt, n in final_counts.most_common():
        print(f"  {qt}: {n:,}")
    if unmapped:
        print(f"WARNING: {sum(unmapped.values()):,} records with unmapped question_type left unchanged:")
        for qt, n in unmapped.most_common(20):
            print(f"  {qt!r}: {n:,}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
