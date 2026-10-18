# Phase 4b: PubTator / NER gene anchoring

Anchors each QA record's mutation to an NCBI Gene ID. The input is
`phase_4/output/m2tqa_gene_recovered.jsonl` (written by `run_gene_recovery.py`). The output,
`m2tqa_tiered_anchored_v3.jsonl`, is the input of `uniprot/scripts/master_grounding_pipeline.py`.

All outputs go to `phase_4/anchoring/output/` by default. Every script takes `--help`, and its
paths can be overridden on the command line.

## Run order

| Step | Script | Reads | Writes | Network |
|------|--------|-------|--------|---------|
| 1 | `make_pmid_list.py` | `m2tqa_gene_recovered.jsonl` | `pmid_list.txt` (unique PMIDs, sorted) | no |
| 2 | `fetch_pubtator_bioc.py` | `pmid_list.txt` | `pubtator_annotations.jsonl`, `failed_pmids.txt` | PubTator3 API |
| 3 | `prepare_ner_input.py` | article corpus + `pmid_list.txt` | `ner_input/batch_NNNN.txt` (1,000 articles each, PubTator format) | no |
| 4 | `fetch_pubtator_batches.py` | `ner_input/batch_*.txt` (PMIDs only) | `ner_output/batch_NNNN.out` | PubTator3 API |
| 5 | `tiered_anchor_v3.py` | `m2tqa_gene_recovered.jsonl`, `pubtator_annotations.jsonl`, `ner_output/`, article corpus | `m2tqa_tiered_anchored_v3.jsonl`, `tiered_stats_v3.json` | no |

Steps 2 and 4 are resumable. Step 2 skips PMIDs that are already fetched, and step 4 skips
batches whose `.out` file is already non-empty. Both steps use the public PubTator3 API, which
needs no API key.

## External inputs

- **Article corpus** (`config.SOURCE_ARTICLES_FILE`, the same file Phase 0 starts from, 4.7 GB).
  It is a JSON object `{pmid: {"article": <full text>, "proteins": [...]}}` that must contain
  every PMID in `pmid_list.txt`. It was built outside this repository; no script here produces
  it. The PMID set comes from a MIRIAD-derived list of PubMed papers (MIRIAD is a public medical
  QA corpus grounded in the literature). The article texts were then backfilled from PMC: about
  26.6K papers got PMC full text and about 2.4K fell back to the PubMed abstract, according to
  the backfill logs kept with the corpus. The texts are needed for two reasons: to find where the mutation appears
  (Tiers 1 and 2), and to write the NER batch files. Step 5 streams the corpus with `ijson`;
  step 3 loads all of it into memory.
- **Tier 1 NER annotations** (`ner_output/`, given to step 5 with `--ner-dir`). These are
  PubTator text-format files named `batch_NNNN.out` (also accepted: `batch_NNNN.txt.out`,
  `*.PubTator`). The script checks for all 31 batches `batch_0000`–`batch_0030`; a missing
  batch only triggers a warning. Each annotation is one line:
  `PMID<TAB>Start<TAB>End<TAB>Mention<TAB>Type<TAB>ConceptID`.
  Only `Type == "Gene"` rows are used. A ConceptID of the form `taxon:gene` sets
  `ncbi_taxon_id`; a bare gene ID leaves the Tier 1 taxon as `"unknown"`.
  - **What the released dataset used:** step 4 produced these files. They are PubTator3
    annotations of each paper's **title and abstract**, not full-text NER.
    The full-text GNormPlus run was started but ran out of memory (Java `OutOfMemoryError`).
    The `.out` files that were finally used match `fetch_pubtator_batches.py` output: real
    PubMed titles, abstract-length passages, bare NCBI Gene IDs, and PMIDs sorted within
    each file.
  - **Full-text alternative (not used for the release):** run GNormPlus on `ner_input/` instead
    of step 4. The original run used the `lanax/gnormplus` Docker image with 10 GB of heap per
    container, `FocusSpecies = All`, `GeneIDMatch = True` and `Normalization2Protein = False`:
    `java -Xmx10G -jar GNormPlus.jar input_data/batch_NNNN.txt output_data/batch_NNNN.txt.out setup.txt`.
    Before that run, the batch files were cleaned: U+00A0 became a space, U+00AD became `-`,
    and titles shorter than 2 characters became `Title for <PMID>`.
    GNormPlus output will give different Tier 1 anchors from the released ones.

## Anchoring logic (`tiered_anchor_v3.py`)

Decisions are made per record, from `(pmid, mutation)` and that PMID's annotations:

1. If the PMID has no PubTator record → `no_pubtator`. If it has no gene annotations (after a
   blacklist of tool names, ACMG terms and similar) → `no_genes`.
2. If PubTator lists exactly one gene → `single_gene`, confidence `high`.
3. **Tier 1** (`tier1_fulltext`): find the mutation in the article text, then take the closest
   Tier 1 NER gene mention within 2,000 characters.
4. **Tier 2** (`tier2_abstract`): find the mutation in the abstract extracted from the article,
   then take the closest PubTator gene mention within 2,000 characters.
5. **Tier 3** (`tier3_most_mentioned`): the PubTator gene with the most mentions, confidence
   `low`. If no gene qualifies → `anchor_failed`.

In Tiers 1 and 2, confidence comes from the distance: `high` up to 200 characters, `medium`
up to 500, `low` beyond that.

## Fields added to each record

Fields are appended in this order:

| Field | Description |
|-------|-------------|
| `anchor_gene_id` | NCBI Gene ID (string), or `null` |
| `anchor_gene_name` | Gene mention or symbol from PubTator or NER, or `null` |
| `anchor_tier` | `single_gene`, `tier1_fulltext`, `tier2_abstract`, `tier3_most_mentioned`, `no_pubtator`, `no_genes` or `anchor_failed` |
| `anchor_confidence` | `high`, `medium`, `low` or `none` |
| `ncbi_taxon_id` | Taxon of the anchor: the paper-level PubTator species (`9606` if present), the NER taxon, `"unknown"`, or `null`. Species filtering happens downstream. |
| `anchor_distance` | Characters between the mutation and the gene mention. Only present for Tier 1 and Tier 2. |
| `vep_input_string` | rsID, genomic string, `NEEDS_TRANSCRIPT:c.…` or `PROTEIN_LEVEL:p.…`. Omitted when the notation is not recognized. |

All other fields pass through unchanged, and records stay in input order.
`tiered_stats_v3.json` counts records per tier, per confidence level and per taxon, so it
reflects only the records that were processed.
