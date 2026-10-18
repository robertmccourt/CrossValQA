# claim_extractor/expansion/

These scripts add articles for genes missing from the pre-expansion dataset
`uniprot/output/m2tqa_gene_corrected.jsonl` (416,188 rows, question types already normalized by
`uniprot/scripts/normalize_question_types.py`). Each batch runs through Phases 0-3 and is then
merged into `uniprot/output/m2tqa_gene_corrected_v2.jsonl` (423,415 rows; +7,227 records, +1,041 genes).
Run the commands from this directory; outputs go to `output/` (gitignored). `run_expansion.py`
calls OpenRouter, so run a `--limit` pilot first. The fetchers need `NCBI_EMAIL` (and optionally
`NCBI_API_KEY`).

```bash
# batch1: ClinVar-cited PMIDs (needs data/variant_summary.txt.gz + data/var_citations.txt from the ClinVar FTP)
python extract_clinvar_targets.py && python fetch_pmc_articles.py && python run_expansion.py --batch batch1
# batch2: remaining protein-coding genes, PMC search
python extract_batch2_targets.py
python fetch_pmc_articles.py --targets output/batch2_targets.jsonl --output output/batch2_articles.json
python run_expansion.py --batch batch2
# batch2_round2: re-fetch genes that produced no QA, 5 papers per gene, excluding already-fetched PMC ids
python find_failed_genes.py --batch batch2
python fetch_pmc_articles.py --targets output/batch2_failed_genes.jsonl --output output/batch2_round2_articles.json \
    --papers-per-gene 5 --exclude-from <earlier articles JSON or PMC-id list>
python run_expansion.py --batch batch2_round2
# humsavar: UniProt humsavar (data/humsavar.txt) + EBI Proteins API evidence
python extract_humsavar_targets.py && python fetch_humsavar_articles.py
python run_expansion.py --batch humsavar --articles output/humsavar_articles.json
# reprocess_filtered: python run_expansion.py --batch reprocess_filtered  (input output/reprocess_filtered.json, see below)
python merge_expansion.py      # -> uniprot/output/m2tqa_gene_corrected_v2.jsonl + expansion_merge_stats.json
```

Gap: the script that built `output/reprocess_filtered.json` (707 articles keyed `<PMID>_<GENE>`,
with an `original_pmid` field) is not part of this repo. Those keys are the source of the mangled
pmids that `release/fix_mangled_pmids.py` removes later.
