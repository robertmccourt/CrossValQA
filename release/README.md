# release/

Scripts that turn the pipeline output into the released MutQA dataset (published on
Hugging Face). Nothing here calls an LLM. Working files are written to
`release/output/` and external inputs are read from `release/input/` (data files are gitignored). Every
script has `--help`, and its paths can be overridden on the command line. Run from the repo root.

## Run order

```bash
# 0. upstream: uniprot/scripts/normalize_question_types.py -> uniprot/output/m2tqa_gene_corrected.jsonl (416,188)
#    expansion: see claim_extractor/expansion/README.md      -> uniprot/output/m2tqa_gene_corrected_v2.jsonl (423,415)

# 1. Working set: drop deterministic-verifier contradictions (labels from audit/) and rows with phase2_answerable=False
python release/build_v2_verified_answerable.py --protein-labels <v2 protein-verifier mutation_verification.jsonl>  # 411,060
python release/build_v3_filtered.py --protein-labels <v3 protein labels> --dna-labels <v3 DNA labels>              # 401,406

# 2. Enriched v2 release file (sequences, structures, isoforms); needs only the v2 dataset
python release/prepare_release_v2.py                                  # projection + ref_sequence (UniProt REST)
python release/enrich_release_v2.py                                   # PDB/AlphaFold ids, positions, features
python release/merge_isoform_enrichment.py --isoform-tsv <TSV>        # matched_isoform_* columns
python release/fix_mangled_pmids.py                                   # 37250416_ARPP21 -> 37250416
#    -> release/output/crossvalqa_release_v2.jsonl (423,415 rows, 32 columns)

# 3. 39-column release, splits, HF staging
python release/build_release_v3.py                     # -> crossvalqa_release_v3.jsonl (401,406)
python release/create_splits.py both                   # -> homology_split/, random_split/
python release/stage_hf.py                             # -> hf_staging/{homology,random,full}/ (Parquet + jsonl.gz)
```

## Release filter and split policy

- Splits use only rows with `grounding == "cross-grounded"` AND `uniprot_primary` AND a non-null
  `ref_sequence`. That is 330,387 of the 401,406 rows, covering 5,437 UniProt accessions, 27,442 PMIDs,
  6,620 genes and 88,587 distinct mutation strings. The `full` config holds all 401,406 rows, unfiltered.
- `homology`: MMseqs2 `easy-cluster --min-seq-id 0.4 -c 0.8 --cov-mode 0` on the wild-type sequences
  (4,160 clusters). Whole clusters go to train or test in an 80/20 ratio (seeded greedy fill on
  record counts, seed 42, no validation split). A GraphPart-style audit then discards every test
  UniProt with >= 40% identity to any train UniProt (15 UniProts, 1,651 records).
  Result: train 264,390 / test 64,346 records (4,240 / 1,182 UniProts), max test-vs-train identity 39.7%.
- `random`: the same filter at the UniProt level (each UniProt is its own cluster), 80/20, seed 42.
  Result: train 263,175 / test 67,212 (4,310 / 1,127 UniProts).
- Each split directory has a `manifest.json` with the parameters, counts and audit results.

## External inputs

- Deterministic verifier labels (protein-side v2 and v3 runs, DNA-side v3 run): these come from
  the verifier under `audit/`. The protein labels are row-aligned to the full v2 file; the DNA labels
  are row-aligned to the 411,060-row working set.
- `release/input/crossvalqa_release.jsonl`: the v1 (416,188-row) enriched release, plus
  `protein_3d_ids.json` and `mutation_locations.json`. `prepare_release_v2.py` and `enrich_release_v2.py`
  carry v1 sequences and structure fields forward from these. They were built outside this repo.
- Isoform TSV for `merge_isoform_enrichment.py`: comes from a collaborator's isoform-matching
  pipeline and is not distributed. The required columns are listed in the script docstring.
- Network: UniProt REST (`rest.uniprot.org`) for the expansion accessions; skip it with `--skip-fetch`.
- MMseqs2 (tested with 18.8cc5c) on `PATH` or `--mmseqs`; `pyarrow` for `stage_hf.py`.
- `get_mutation_dna.py` is a helper for dataset users, not a build step. It needs `pysam` and GRCh38
  plus UniProt CDS tables under `$CROSSVALQA_REF_DIR` (default `reference_data/`). In
  `test_get_mutation_dna.py`, the parser tests run without these files; the genome tests are skipped.

## Reproducing the published splits

The published Hugging Face splits are canonical. Rerunning `create_splits.py homology` does not
reproduce them exactly: MMseqs2 writes the same 4,160 clusters in a different order on each run,
and the seeded shuffle runs over that order. A rerun here gave 264,653 / 65,608 records and dropped
7 test UniProts in the audit (published: 15). To get the exact published assignment, use each
config's `uniprot_manifest.json` from the Hugging Face repo, which lists the UniProt accessions
assigned to train and to test. The `random` split does reproduce exactly.

## Known issue

`build_release_v3.py` joins on `pmid`, but the working set still has 1,637 rows with mangled
`<PMID>_<GENE>` pmids (213 distinct, all from the `reprocess_filtered` expansion batch). The enriched
file was already fixed, so those rows find no match: their pmid stays mangled and their sequence and
structure fields are null. As a result they are missing from both splits (1,486 of them are
cross-grounded) and appear with the mangled pmid in the `full` config. The released files are
reproduced as built, without fixing this.
