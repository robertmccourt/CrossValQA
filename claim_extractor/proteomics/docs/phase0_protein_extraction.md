# Proteomics Phase 0: Protein Extraction and Verification

The proteomics proof of concept replaces the mutation regex of the main
pipeline with PubTator3 gene annotations. Phase 0 has four steps (steps 1-4
of `../README.md`): download articles, extract and map proteins, prefilter by
functional language, and verify with an LLM. Steps 1-3 call no LLM.

All paths below are relative to `claim_extractor/proteomics/`.

| Step | Script | Output |
|---|---|---|
| 1 | `get_articles.py` | `proteomics_dataset/all_articles.json`, `proteomics_dataset/<category>/articles.json` |
| 2 | `phase_0/extract_proteins.py` | `phase_0/output/extracted_proteins.json` |
| 3 | `phase_0/prefilter_proteomics.py` | `phase_0/output/prefiltered_proteins.json` |
| 4 | `phase_0/verify_proteomics.py` | `phase_0/output/verified_proteins.json`, `verification_stats.json`, `failed_batches.jsonl` |

## Step 1: articles (`get_articles.py`)

`PROTEOMICS_DATA` lists PMIDs under five categories
(`enzymology_and_kinetics`, `structural_biology_and_motifs`,
`protein_protein_interactions`, `post_translational_modifications`,
`clinical_proteomics_and_biomarkers`); 26 unique PMIDs (40545624 is in two
categories). For each PMID, `fetch_article` gets PubMed metadata via
E-utilities `efetch`, looks up a PMC ID via `elink`, and fetches the PMC XML
full text. Articles without PMC full text are skipped. Requests are spaced by
0.35-0.4 s. The contact e-mail comes from `NCBI_EMAIL` (`get_ncbi_email`).

Recorded run: 22 full-text articles, 4 skipped (40545624, 41381349, 41524014,
41569177), per `all_articles.json` metadata.

## Step 2: PubTator3 extraction and UniProt mapping (`extract_proteins.py`)

1. `fetch_pubtator_annotations` requests BioC JSON from
   `PUBTATOR_BIOC_URL` (`https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson`)
   in batches of 10 PMIDs, with `REQUEST_DELAY = 0.5` s between requests.
2. `extract_gene_annotations` keeps annotations whose `infons.type` is `gene`
   and whose `identifier` (or `ncbi_gene`) is numeric (comma-separated IDs are
   split). Annotations are merged per gene ID into a `GeneAnnotation`
   (`gene_id`, `symbol`, `name`, `taxon_id`, `mentions`, `offsets`). The
   symbol is `infons.name` when present.
3. `map_gene_ids_to_uniprot` submits UniProt ID-mapping jobs
   (`from: GeneID`, `to: UniProtKB`, 100 IDs per job) with
   `submit_uniprot_mapping`, polls them with `poll_uniprot_mapping` (up to 20
   attempts, 1 s apart), and keeps the first hit per gene (`primaryAccession`,
   `uniProtkbId`, recommended full name).
4. Genes still unmapped go to `batch_lookup_uniprot`, which searches
   `(xref:geneid-<id>)` (5 results) and prefers a reviewed (Swiss-Prot) entry.

Offsets are kept in `GeneAnnotation` but are not written to the output.

`extracted_proteins.json` has `metadata`, `protein_mappings` (gene ID ->
`ProteinMapping`: `gene_id`, `gene_symbol`, `taxon_id`, `uniprot_accession`,
`uniprot_entry_name`, `protein_name`) and `articles` (PMID -> `pmid`, `title`,
`gene_annotations`, `article_text`, `abstract`). Each gene annotation carries
the UniProt fields.

Recorded run (`metadata`): 22 articles, all with PubTator data; 67 gene
annotations, 67 unique genes, 65 mapped to UniProt.

## Step 3: prefilter (`prefilter_proteomics.py`)

Candidate names per article (`get_protein_names_from_annotations`) are each
gene's symbol, every mention string, and the UniProt protein name. One gene
can therefore produce several candidates, which is why "discussed" counts can
exceed the number of genes.

`prefilter_proteins(article_text, proteins, window=1)`:

1. `is_substantive_protein` drops names shorter than 3 characters, names in a
   generic-term set (`protein`, `enzyme`, `receptor`, `factor`, `kinase`,
   `antibody`, `antigen`, `marker`, `biomarker`, `target`, `substrate`,
   `inhibitor`, `activator`), and names occurring fewer than 2 times
   (case-insensitive) in the text. These go to `filtered_generic`.
2. `split_into_sentences` (this file's own version: protects initials,
   `et al.`, common abbreviations including `e.g.`/`i.e.`, and decimals).
3. `find_protein_sentence_indices` finds sentences containing the name
   (case-insensitive; a space in the name also matches a hyphen).
4. `check_discussion_proximity` checks the protein sentences and one sentence
   on either side. A sentence counts as evidence if it matches
   `PROTEOMICS_REGEX` (built from `PROTEOMICS_PATTERNS`: characterization,
   enzymatic activity and kinetics, binding and interactions, PTMs,
   structure, localization and expression, functional roles,
   biomarker/clinical, stability and physical properties,
   recombinant/variant) and does not match `JUST_MENTIONED_REGEX`
   (`JUST_MENTIONED_PATTERNS`: protocol/kit/antibody usage, `see`/`cf.`/
   `reviewed in`, `previously reported`/`well-known`, `marker|control|standard for|of`).

Names with evidence go to `likely_discussed` (`protein`, `protein_sentences`,
`evidence_sentences`, up to three `evidence_text`); the rest go to
`likely_just_mentioned`. Discussed entries are enriched with `gene_id`,
`taxon_id`, `uniprot_accession`, `uniprot_entry_name` and
`protein_name_canonical` by looking the name up in
`build_gene_identifier_map`, which is keyed by symbol and mentions only; a
candidate that came from the UniProt protein name gets no identifiers.

`prefiltered_proteins.json` has `metadata`, `protein_mappings` (carried
forward) and `articles` (PMID -> `title`, `likely_discussed`,
`likely_just_mentioned`, `filtered_generic`, `gene_annotations`, `stats`,
`sentences`).

Recorded run (`metadata`): 74 discussed (59 with a UniProt accession),
7 just mentioned, across 22 articles.

## Step 4: LLM verification (`verify_proteomics.py`)

`process_article` sends each article's `likely_discussed` list in batches of
`--batch-size` (default 5) to `verify_protein_batch`. The prompt
(`SYSTEM_PROMPT`, `USER_PROMPT_TEMPLATE`) lists the protein names and the
numbered sentences from `build_numbered_sentences` (protein and evidence
sentence indices expanded by `context_window=2`, capped at 100, formatted
`[idx] sentence`). The model returns per protein `protein`, `is_discussed`,
`consequence_summary`, `evidence_sentences` and `confidence`, and is told to
state only what the cited sentences say and to mark methods/control mentions
as not discussed. `parse_llm_response` keeps results whose `protein` is in the
batch; identifiers from step 3 are attached as `ProteinIdentifiers`.

| Setting | Value |
|---|---|
| Model | `--model`, default `google/gemini-2.0-flash-001` |
| Client | shared `OpenRouterClient` in `claim_extractor/openrouter.py` |
| `temperature` / `max_tokens` | 0.1 / 2000 |
| Concurrency | `--workers`, default 10 |
| Timeout | 120 s per request |
| Retries | 3; on HTTP 429 wait 5 s x attempt, on 5xx wait 2 s x attempt, other 4xx not retried |

There is no checkpointing or resume. `failed_batches.jsonl` (`pmid`,
`proteins`, `error`, `timestamp`) is deleted at the start of each run.

`verified_proteins.json` maps PMID to `verified_proteins`,
`unverified_proteins` (each entry: the five LLM fields plus `gene_id`,
`taxon_id`, `uniprot_accession`, `uniprot_entry_name`,
`protein_name_canonical`), `just_mentioned`, `all_identifiers`
(`{"gene_annotations": [...]}`) and `error`.

Recorded run (`verification_stats.json`, model `google/gemini-2.0-flash-001`):

| Field | Value |
|---|---|
| `total_articles` | 22 |
| `articles_with_verified_proteins` | 10 |
| `total_verified_proteins` | 71 |
| `verified_with_identifiers` (gene ID or UniProt accession) | 57 |
| `total_unverified_by_llm` | 3 |
| `total_just_mentioned` | 7 |
| `total_errors` | 0 |

Of the 71 verified proteins in `verified_proteins.json`, 57 have a
`gene_id` and 56 have a `uniprot_accession`.

## Differences from the mutation Phase 0

| | Mutations (`claim_extractor/phase_0/`) | Proteins (this folder) |
|---|---|---|
| Entity detection | `mutation_regex.py` | PubTator3 gene annotations |
| Database IDs | added after Phase 3 (Phase 4 gene recovery, then UniProt grounding) | NCBI Gene and UniProt in step 2 |
| Prefilter evidence | `CONSEQUENCE_PATTERNS`, +-1 sentence | `PROTEOMICS_PATTERNS` minus `JUST_MENTIONED_PATTERNS`, +-1 sentence |
| LLM context | +-2 sentences, max 100 | +-2 sentences, max 100 |
| Batch size / concurrency | 10 / 50 | 5 / 10 |
| Checkpoint / resume | yes | no |
