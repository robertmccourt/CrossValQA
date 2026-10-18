# CrossValQA

Code for **MutQA: A Cross-Validated Q&A Dataset for Genetic Mutations** (NeurIPS 2026).

CrossValQA generates question–answer pairs from scientific articles and keeps only the pairs that
two independent LLMs can each justify from the source text. MutQA is the dataset built by applying
it to PubMed literature on genetic variants.

**Dataset:** [huggingface.co/datasets/random987654321/MutQA](https://huggingface.co/datasets/random987654321/MutQA) ·
**Prompts:** [`prompts/README.md`](prompts/README.md) ·
**License:** code MIT, dataset CC BY 4.0

## How it works

```
 article ──► Phase 0 ──► Phase 1 ──────────► Phase 2 ──────────► Phase 3 ──────────────► Phase 4 + grounding
            find and     QA generation       independent answer   cross-validation       gene anchoring,
            verify       Q, A₁, R₁           A₂, R₂ (never sees   does R₁ support A₂?    UniProt mapping,
            mutations    (Gemini 2.5 Flash)  A₁ or R₁; DeepSeek)  does R₂ support A₁?    evidence sentences
                                                                  (Llama 3.1 8B)
```

- **Information isolation.** The second model answers the question from the article alone; it
  never sees the first model's answer or rationale.
- **Bidirectional cross-validation.** A third model is shown each rationale paired with the
  *other* model's answer, without knowing which model wrote what. A pair is **cross-grounded**
  when both directions pass, **grounded** when one passes, and **ungrounded** otherwise.
- **Citations.** Every answer and rationale cites sentence IDs in the source article; the cited
  sentences are stored with the record.

## The dataset

The released train/test splits contain only cross-grounded records with a UniProt accession and a
reference protein sequence:

| | |
|:--|--:|
| QA pairs | 330,387 |
| PubMed articles | 27,442 |
| Genes | 6,620 |
| Distinct mutation strings | 88,587 |
| UniProt accessions | 5,437 |

| Config | Split | Records | Notes |
|:--|:--|--:|:--|
| `homology` (default) | train / test | 264,390 / 64,346 | MMseqs2 clusters at 40% identity; no test protein ≥ 40% identical to a train protein |
| `random` | train / test | 263,175 / 67,212 | random 80/20 at the protein level; ablation baseline |
| `full` | full | 401,406 | all grounding tiers, unfiltered; not for training |

```python
from datasets import load_dataset

ds = load_dataset("random987654321/MutQA")             # homology split
full = load_dataset("random987654321/MutQA", "full")   # every record, all grounding tiers
```

Each record keeps both models' answers, rationales and evidence IDs, both cross-validation verdicts
with the judge's reasoning, the cited sentences, and protein context (UniProt accession, reference
and mutated sequence, position, AlphaFold/PDB identifiers). The dataset card describes all 39 fields.

## Repository layout

```
claim_extractor/          the CrossValQA pipeline
  config.py                 all default paths and models (relative to the repo root)
  openrouter.py             shared OpenRouter client
  phase_0/                  mutation extraction (regex) and LLM verification
  phase_1/                  QA generation
  phase_2/                  independent answering
  phase_3/                  cross-validation judge (+ recovery pass for failed calls)
  phase_4/                  final record assembly, cited sentences, gene recovery
    anchoring/              PubTator-based gene anchoring
  expansion/                extra articles for genes missing from the first run (ClinVar, humsavar)
  proteomics/               proof-of-concept adaptation to proteomics papers
  docs/                     Phase 0 technical notes
uniprot/scripts/          gene → UniProt grounding, gene-assignment correction, question-type normalization
audit/                    deterministic (LLM-free) claim auditor: checks answers against the sequence and genome
release/                  audit filter, release projection, splits, Parquet export
prompts/                  every LLM prompt used to build and evaluate the dataset
```

## Reproducing the pipeline

Most users only need the dataset above. Rerunning the pipeline costs LLM API credits and needs
several large external inputs (listed below).

### Setup

```bash
git clone https://github.com/robertmccourt/CrossValQA.git && cd CrossValQA
pip install -r requirements.txt          # Python 3.10+
cp .env.example .env                     # add OPENROUTER_API_KEY (and NCBI_EMAIL for the fetchers)
```

Every script has `--help`; defaults come from `claim_extractor/config.py` and can be overridden on the
command line. Data files are written next to the scripts (`*/output/`, `database/`) and are gitignored.
Use `--limit` for a small pilot before any full LLM run.

### Steps

```bash
# Phase 0: mutations
python claim_extractor/phase_0/batch_extractor.py        # regex extraction over the article corpus
python claim_extractor/phase_0/build_article_set.py      # drop outliers / short / mutation-free articles
python claim_extractor/phase_0/verify_mutations.py       # LLM check that each mutation is discussed

# Phases 1–3: generate, answer independently, cross-validate
python claim_extractor/phase_1/phase1_generator_qa.py --workers 50
python claim_extractor/phase_2/phase2_generator.py --workers 50
python claim_extractor/phase_3/phase3_judge.py --workers 50 --once
python claim_extractor/phase_3/phase3_recovery.py        # retry judge calls that failed
python claim_extractor/phase_3/merge_recovery.py

# Phase 4: assemble records, attach cited sentences, recover genes, anchor genes
python claim_extractor/phase_4/build_final.py --gene-map <pmid_mutation_gene_map.json>
python claim_extractor/phase_4/extract_cited_sentences.py
python claim_extractor/phase_4/run_gene_recovery.py
#   then the five anchoring steps in claim_extractor/phase_4/anchoring/README.md

# UniProt grounding
python uniprot/scripts/master_grounding_pipeline.py
python uniprot/scripts/correct_gene_assignments.py --pmid-resource <ClinVar PMID resource JSON>
python uniprot/scripts/normalize_question_types.py       # -> 416,188 records

# Expansion, deterministic audit, release and splits
#   claim_extractor/expansion/README.md                  # -> 423,415 records
#   audit/README.md                                      # verified / contradicted label per claim
#   release/README.md                                    # -> 401,406 records, then the splits
```

### External inputs not in this repository

| Input | Used by | Notes |
|:--|:--|:--|
| Article corpus `database/pubmed_to_proteins_with_articles.backfilled.json` | Phase 0, anchoring | `{pmid: {"article": text, ...}}`; PMC full text (abstract fallback) for a MIRIAD-derived PMID list |
| `pmid_mutation_gene_map.json` | `build_final.py` | `{"<pmid>\|<mutation>": gene}`; a noisy prior that later steps override |
| PubTator3 annotations | anchoring | fetched by the anchoring scripts from the public API |
| NCBI `Homo_sapiens.gene_info`, UniProt `*_idmapping.dat.gz` | grounding | public downloads; see the script docstrings |
| ClinVar `variant_summary.txt.gz`, `var_citations.txt`; UniProt `humsavar.txt` | expansion | public downloads |
| GRCh38 primary assembly, UniProt ↔ hg38 CDS tables | `audit/` (DNA side) | see `audit/README.md` |
| Isoform-matching TSV, v1 enriched release | `release/` | not distributed; see `release/README.md` |
| MMseqs2 | `release/create_splits.py` | tested with 18.8cc5c |

### Known limitations

- Rerunning `release/create_splits.py homology` does not reproduce the published homology split
  exactly (MMseqs2 cluster order varies between runs). Use the `uniprot_manifest.json` in each
  Hugging Face config for the exact assignment.
- 1,637 records in the `full` config keep a `<PMID>_<GENE>` pmid from one expansion batch and
  lack sequence fields, so they are excluded from the splits. See `release/README.md`.
- LLM outputs are not deterministic, so a rerun of Phases 0–3 gives a similar but not identical set.

## Proteomics adaptation

[`claim_extractor/proteomics/`](claim_extractor/proteomics/README.md) applies the same pipeline to
proteomics papers as a proof of concept: 22 articles, 71 verified proteins, 210 QA pairs, 89.0%
cross-grounded.

## Citation

```bibtex
@inproceedings{mutqa2026,
  title     = {MutQA: A Cross-Validated Q\&A Dataset for Genetic Mutations},
  author    = {McCourt, Robert and Macaulay, Oladimeji and Arredondo, David and
               Tafoya, Luis E. and Edeh, Oluoma and Virupakshappa, Kushal and Hu, Yue and
               Nguyen Bach, My and Almeida Ruiz, Enrique and Poshiya, Shrey and
               Hudgens, Debjani P. and Kundnani, Deepali and Sahu, Avinash},
  booktitle = {Advances in Neural Information Processing Systems},
  year      = {2026}
}
```

## Acknowledgments

The deterministic auditor in `audit/` was written by Oladimeji Macaulay.

## License

Code: [MIT](LICENSE). Dataset: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
