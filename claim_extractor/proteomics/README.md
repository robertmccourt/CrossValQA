# Proteomics adaptation (proof of concept)

A small adaptation of the CrossValQA pipeline from mutations to proteins. The
scripts here are forked copies of the main `claim_extractor/phase_*` scripts with
protein-specific prompts and entity handling. They reuse
`claim_extractor/phase_0/prefilter.py` (sentence splitting) and, for the Phase 0
verifier and Phase 1, the shared `claim_extractor/openrouter.py` client.

## Requirements

- `.env` at the repo root with `OPENROUTER_API_KEY` (phases 0-verify, 1, 2, 3)
  and `NCBI_EMAIL` (`get_articles.py`). The phase scripts read
  `OPENROUTER_API_KEY` from the environment or `--api-key`.
- Network access to NCBI E-utilities, PubTator3 and UniProt (steps 1-2).

## Run order

Run from any directory; all paths are relative to this folder.

| Step | Command | Input | Output |
|---|---|---|---|
| 1 | `python get_articles.py` | PMID list in the script | `proteomics_dataset/all_articles.json` |
| 2 | `python phase_0/extract_proteins.py` | `all_articles.json` | `phase_0/output/extracted_proteins.json` |
| 3 | `python phase_0/prefilter_proteomics.py` | `extracted_proteins.json` | `phase_0/output/prefiltered_proteins.json` |
| 4 | `python phase_0/verify_proteomics.py` | `prefiltered_proteins.json` | `phase_0/output/verified_proteins.json`, `verification_stats.json` |
| 5 | `python phase_1/phase1_generator_qa.py` | `verified_proteins.json`, `all_articles.json` | `phase_1/output/phase1_qa_pairs.jsonl` |
| 6 | `python phase_2/phase2_generator.py` | `phase1_qa_pairs.jsonl`, `all_articles.json` | `phase_2/output/phase2_qa_pairs.jsonl` |
| 7 | `python phase_3/phase3_judge.py` | `phase2_qa_pairs.jsonl` | `phase_3/output/phase3_validated.jsonl`, `stats.json` |

Steps 1-3 call no LLM. Steps 4-7 call OpenRouter; default models are
`google/gemini-2.0-flash-001` (4), `google/gemini-2.5-flash-lite` (5),
`deepseek/deepseek-chat` (6) and `meta-llama/llama-3.1-8b-instruct` (7), each
overridable with `--model`. Phases 1-3 support `--resume` and `--limit`.
Output data files are not tracked in git.

Phase 0 (steps 1-4) is described in [docs/phase0_protein_extraction.md](docs/phase0_protein_extraction.md).

## Results

| Stage | Count |
|---|---|
| PMIDs requested (5 categories) | 26 |
| Full-text articles downloaded | 22 (4 skipped, no PMC full text) |
| Articles with at least one verified protein | 10 |
| Verified proteins | 71 |
| QA pairs (Phase 2 / Phase 3) | 210 |
| Cross-grounded | 187 (89.0%) |
| Grounded only | 10 (4.8%) |
| Ungrounded | 13 (6.2%) |
| Grounded or cross-grounded | 197 (93.8%) |
