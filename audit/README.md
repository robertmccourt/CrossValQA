# audit/

LLM-free deterministic claim auditor for MutQA. For each QA row it derives ground-truth facts
from the mutation itself, regex-extracts factual claims from both the phase-1 answer (`qagen`) and
the phase-2 answer (`llm2`), and labels each claim `verified` or `contradicted`. No LLM or API calls;
only the two UniProt fetch scripts use the network.

- **Protein side** (`verifier/mutation_verifier_dual.py`): facts from `ref_sequence` + the HGVS `p.`
  descriptor (wild-type/mutant residue and position, mutation kind, protein length, charge/polarity/
  hydropathy/size change, N-/C-terminal region, special residues) plus UniProt feature annotation
  (domain, site, modified residue, glycosylation, disulfide, transmembrane, topology, signal peptide).
  Feature checks only emit `contradicted` on a direct mismatch at the position, never on missing annotation.
- **DNA side** (`verifier/mutation_verifier_dna.py`): for DNA-notation mutations (`c.`, `g.`, or bare forms such as `89T>C`), maps the variant
  onto GRCh38 through a UniProt-to-hg38 CDS map and checks ref/alt base, mutation kind,
  transition/transversion, predicted wild-type/mutant amino acid, premature stop and splice-region
  claims; the protein-side residue-chemistry checks are reused for exonic non-synonymous substitutions.
- Row label (`qagen_label`, `llm2_label`): `contradicted` if any check is contradicted, else `verified`
  if any check fired, else `not_applicable`. Hedged, hypothetical and negated phrasings are skipped.

## Run order (from the repo root; outputs go to `audit/output/`, which is gitignored)

```bash
python audit/fetch_uniprot_sequences.py      # -> audit/output/uniprot_sequences.json (UniProt REST)
python audit/fetch_uniprot_features.py       # -> audit/output/uniprot_features.json  (UniProt REST)
python audit/build_gold_rows_v2.py           # v2 + sequences -> audit/output/_gold_rows_v2.jsonl

# protein labels, v2 run (no UniProt feature checks) and v3 run (with them)
python audit/verifier/mutation_verifier_dual.py --input audit/output/_gold_rows_v2.jsonl \
    --out_dir audit/output/protein_v2 --no_features
python audit/verifier/mutation_verifier_dual.py --input audit/output/_gold_rows_v2.jsonl \
    --out_dir audit/output/protein_v3

# DNA labels, on the working set written by release/build_v2_verified_answerable.py
python audit/verifier/mutation_verifier_dna.py \
    --input release/output/m2tqa_gene_corrected_v2_verified_answerable.jsonl --out_dir audit/output/dna_v3
```

Then `release/build_v2_verified_answerable.py --protein-labels audit/output/protein_v2/mutation_verification.jsonl`
and `release/build_v3_filtered.py --protein-labels audit/output/protein_v3/mutation_verification.jsonl
--dna-labels audit/output/dna_v3/mutation_verification_dna.jsonl` (see `release/README.md`). Both drop a
row when `qagen_label == "contradicted"` or `llm2_label == "contradicted"`, keyed on `row_id`: the
0-based line index of the v2 file (protein) or of the working set (DNA). Every script takes `--help`.

## Outputs

| file | rows | fields used by `release/` |
|---|---|---|
| `protein_*/mutation_verification.jsonl` | one per input row (unparseable rows: `parse_failed: true`) | `row_id`, `qagen_label`, `llm2_label` |
| `dna_v3/mutation_verification_dna.jsonl` | one per DNA-notation row | `row_id`, `qagen_label`, `llm2_label` |

Each row also carries `facts` and the per-claim `qagen_checks` / `llm2_checks`
(`{claim_type, claimed, ref, status}`); each run also writes a `*_summary.md` with
coverage and verified/contradicted rates by question bucket, split and claim type.

On the 423,415-row v2 file, the v3 protein run could check 79,639 phase-1 answers (149,811 claims, 95.85% verified)
and 77,394 phase-2 answers (146,405 claims, 96.34% verified); 246,291 mutations do not parse as HGVS `p.`.
On the 232,311 DNA-notation rows of the working set, the DNA run could check 174,698 phase-1 answers
(454,804 claims, 98.19% verified) and 175,944 phase-2 answers (465,203 claims, 98.03% verified).

## External data (not distributed)

The DNA verifier needs `pysam` and these files under `$CROSSVALQA_REF_DIR` (default `reference_data/` at the repo root, shared with `release/get_mutation_dna.py`).
It loads them through `release/get_mutation_dna.py`, the same helper that is shipped for dataset users.

```
reference_genome/GRCh38.primary_assembly.genome.fa   (+ .fai; samtools faidx)
uniprot/hg38_uniprot_sequences.csv                   gene -> UniProt accession, chrom, strand, isMain
uniprot/uniprotGtfHg38_swiss.csv                     CDS exon blocks per Swiss-Prot entry (GTF)
uniprot/uniprotGtfHg38_trembl.csv                    CDS exon blocks per TrEMBL entry (GTF)
```

- GRCh38: the GENCODE primary-assembly FASTA (`GRCh38.primary_assembly.genome.fa.gz` from any GENCODE
  human release; chr-prefixed names, 194 sequences).
- UniProt-to-hg38 tables (tab-separated despite the `.csv` suffix): exports of the UCSC Genome Browser
  hg38 "UniProt" track. The GTF source column reads `hg38_unipAliSwissprot` / `hg38_unipAliTrembl`. The
  exports date from 2024–2025; the exact release was not recorded.

The protein verifier needs only the two UniProt caches written by the fetch scripts.

## Credit

The deterministic auditor was written by Oladimeji Macaulay.
Copyright (c) 2026 Oladimeji Macaulay; released under the MIT license in [`LICENSE`](../LICENSE).
