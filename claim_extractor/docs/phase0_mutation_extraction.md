# Phase 0: Mutation Extraction and Verification

Phase 0 finds mutation strings in full-text articles with regular expressions,
keeps the articles that contain them, and asks an LLM which of those mutations
the paper discusses functionally. Its output, `database/verified_mutations.json`,
is the input to Phase 1.

Code: `claim_extractor/phase_0/`. Paths come from `claim_extractor/config.py`.

## Run order

| Step | Script | Input | Output |
|---|---|---|---|
| - | (external, not in this repo) | PubMed/PMC download | `database/pubmed_to_proteins_with_articles.backfilled.json` (`SOURCE_ARTICLES_FILE`) |
| 1 | `batch_extractor.py` | `SOURCE_ARTICLES_FILE` | `database/pubmed_to_proteins_with_articles.mutations_extracted.json` (`EXTRACTED_ARTICLES_FILE`), `database/extraction_stats.json` |
| 2 | `build_article_set.py` | `EXTRACTED_ARTICLES_FILE` | `database/articles_with_mutations.json` (`ARTICLES_FILE`) |
| 3 | `verify_mutations.py` | `ARTICLES_FILE` | `database/verified_mutations.json`, `database/verification_stats.json`, `database/failed_batches.jsonl`, `database/checkpoints/` |

Steps 1 and 2 call no LLM. `mutation_regex.py` and `prefilter.py` are libraries
used by steps 1 and 3.

The source file maps PMID to an entry with an `article` field (full text).
`batch_extractor.py` overwrites each entry's `proteins` field with the
extracted mutation list; downstream scripts read mutations from `proteins`.

## Step 1: regex extraction (`batch_extractor.py`, `mutation_regex.py`)

`batch_extractor.py` runs `process_single_article` over every entry with a
`ProcessPoolExecutor` (`--workers`, default `os.cpu_count()`; `executor.map`
with `chunksize=256`). Each article goes through
`extract_all_variants(text, apply_section_filter=True)` and then
`get_mutations_for_llm(results)["mutations"]`.

### Section filtering

`filter_article_sections` splits the text on headings matched by
`SECTION_SPLIT_RE` (a line in ALL CAPS or Title Case between newlines).
Headings are normalized with `normalize_heading`, and:

- headings longer than 60 characters are skipped;
- sections in `DROP_SECTIONS` (references, acknowledgments, funding, author
  contributions, competing interests / conflicts of interest, supplementary
  references/material/methods, extended methods, materials) are dropped;
- sections in `KEEP_SECTIONS` (abstract, introduction, background, methods,
  results, discussion, conclusion(s), objectives, findings, interpretation,
  measurements and main results) are kept, along with any text before the
  first heading;
- sections in neither set are dropped.

The full original text is used instead if no heading is found, no
`KEEP_SECTIONS` heading is found, or the filtered text is shorter than 30% of
the original. `strip_emails_and_urls` then removes e-mail addresses and URLs.

### Non-HGVS patterns (`extract_nonhgvs_variants`)

Preprocessing replaces `\u0001` with a space and splits packed tokens where a
`3` sits between two letters (`C292T3R98W` -> `C292T R98W`). The gene prefix
is `(?:[A-Za-z]{2,10}\d{0,5})` with an optional trailing `p`.

| Key | Example | Notes |
|---|---|---|
| `gene_aa` | `PEX26R98W`, `Pex26pL44P`, `VWF-N528S`, `BRAF V600E` | `pat_gene_aa` (optional hyphen) plus `pat_gene_space_aa` (gene of 3+ uppercase chars, a space, then the change); `BLOCKLIST` removes cell-line names (`HEK293T`, `COS7`, `NIH3T3`, ...) |
| `gene_int` | `PEX26intG231T` | |
| `gene_indel` | `PEX26G255insT`, `PEX26C861del` | |
| `gene_compound` | `PEX26L153V/C861del` | |
| `nt_sub` | `C292T` | not gated |
| `nt_indel` | `C861del` | kept if the suffix of a `gene_indel` hit, seen 2+ times, or inside a compound |
| `compound` | `L153V/C861del` | not gated |
| `bare_sub_arrow` | `1234G>C` | not gated |
| `bare_aa` | `R98W` | kept only if it is the suffix of a `gene_aa` hit |
| `bare_int` | `intG231T` | kept only if it is the suffix of a `gene_int` hit |

Gene-anchored tokens are normalized by `normalize_gene_token` (remove
whitespace, drop the protein `p` before the change, uppercase), so
`Pex26pR98W` becomes `PEX26R98W`.

### HGVS patterns (`extract_hgvs_variants`)

`þ` is replaced with `+` first. `c.` and `p.` must not follow a letter or digit.

| Key | Examples |
|---|---|
| `c_splice` | `c.768+358C>T`, `c.859-9T>C`, `c.6729+5_+19del` |
| `c_sub` | `c.1040C>T`, `c.1040C→T`, and the non-standard `c.G2578A` |
| `c_del` / `c_dup` / `c_ins` | `c.4248_4250del`, `c.247_250dup`, `c.123_124insACGT` |
| `p_stop` | `p.W855*`, `p.Trp855*`, `p.R408X` |
| `p_complex` | `p.[L541P;A1038V]` |
| `p_3letter` | `p.Gly12Asp`, `p.Arg97Ter`, `p.Phe508del`, `p.Arg97Glyfs*12` |
| `p_1letter` | `p.G12D`, `p.S47fs`, `p.K814del`, `p.R98fsX12`, `p.delF508` |
| `rsids` | `rs12345` |

### Final list (`get_mutations_for_llm`)

The `mutations` list is the sorted union of:

- gene-anchored hits (`gene_aa`, `gene_int`, `gene_indel`, `gene_compound`);
- all HGVS `c.` and `p.` hits;
- `nt_sub` hits, except those that are the position part of a `gene_int` hit;
- `bare_sub_arrow` hits;
- rsIDs.

`nt_indel`, `compound`, `bare_aa` and `bare_int` are extracted but are not
added to this list.

`batch_extractor.py` writes `extraction_stats.json` with the totals, the
per-article count distribution, up to 1,000 zero-mutation PMIDs and all error
PMIDs. In the released run it reported 114,831 articles and 1,888,872
extracted mutations.

## Step 2: article set (`build_article_set.py`)

Filters run in this order:

1. drop articles with more than `MAX_MUTATIONS = 200` extracted mutations;
2. drop articles whose `article` text is shorter than `MIN_ARTICLE_CHARS = 5000`;
3. drop articles with zero extracted mutations (optionally saved with
   `--no-mutations-output`).

Released run (the script's console log):

| Stage | Articles |
|---|---|
| Input | 114,831 |
| After step 1 (875 removed) | 113,956 |
| After step 2 (9,931 removed) | 104,025 |
| After step 3 (17,708 removed) | 86,317 |

The resulting `articles_with_mutations.json` has 86,317 articles and
1,539,316 mutations in their `proteins` lists (counted from the file).

## Step 3: LLM verification (`verify_mutations.py`, `prefilter.py`)

### Prefilter (`prefilter_mutations`, no LLM)

For each article, the full `article` text (not the section-filtered text) and
its `proteins` list are passed to `prefilter_mutations(article_text, mutations, window=1)`:

1. `is_likely_mutation` splits the list into mutations and
   `filtered_gene_names`. Strings matching a mutation pattern (`c.`, `p.`,
   `rs`, `>`, `fs`, `del`/`ins`/`dup`, `IVS`, hyphenated gene-anchored,
   3-letter/1-letter amino-acid changes) are kept; strings matching
   `GENE_NAME_PATTERNS` or the all-caps heuristic are filtered.
2. `split_into_sentences` splits the text with a regex that protects
   initials, `et al.`, `Fig.`/`vs.`-style abbreviations and digits followed
   by a period.
3. `find_mutation_sentence_indices` finds the sentences containing each
   mutation (case-insensitive literal match).
4. `check_consequence_proximity` looks for `CONSEQUENCE_PATTERNS` (16
   patterns, e.g. `loss[- ]of[- ]function`, `dominant[- ]negative`,
   `(is|was|were|classified\s+as)\s+(likely\s+)?(pathogenic|benign|deleterious)`)
   in the mutation sentences and the sentence on either side.

A mutation with at least one such sentence goes to `likely_discussed` (with
`mutation_sentences`, `evidence_sentences`, and up to three `evidence_text`
strings); the rest go to `likely_just_mentioned`. Only `likely_discussed` is
sent to the LLM.

### LLM call

`process_article` sends `likely_discussed` in batches of `--batch-size`
(default 10) to `verify_mutation_batch`. The prompt (`SYSTEM_PROMPT`,
`USER_PROMPT_TEMPLATE`) lists the mutations and the numbered sentences from
`build_numbered_sentences`: every mutation and evidence sentence index expanded
by `context_window=2` on each side, capped at the first 100 indices, formatted
as `[idx] sentence`. The model returns, per mutation, `mutation`,
`is_discussed`, `consequence_summary`, `evidence_sentences` and `confidence`,
and is told to state only what the cited sentences say and to mark
`is_discussed` false when in doubt. `parse_llm_response` keeps only results
whose `mutation` is in the batch.

If an entry has a `gene` key (used by `claim_extractor/expansion/`), the
`USER_PROMPT_TEMPLATE_GENE` prompt is used instead and the result also carries
`gene_match`.

| Setting | Value |
|---|---|
| Model | `--model`, default `DEFAULT_PHASE0_VERIFY_MODEL = "google/gemini-2.5-flash-lite"` (used for the released run) |
| Endpoint | `https://openrouter.ai/api/v1/chat/completions` (`OpenRouterClient` in this file) |
| `temperature` / `max_tokens` | 0.1 / 2000 |
| Concurrency | `--workers`, default 50 (asyncio semaphore) |
| Articles per task chunk | 2000 (`chunk_size`) |
| Timeout | 120 s per request |
| Retries | 3; on HTTP 429 wait 5 s x attempt, on 5xx wait 2 s x attempt, other 4xx not retried |
| Checkpoint | every `--checkpoint-interval` articles (default 1000) to `database/checkpoints/checkpoint_<n>.json` |
| Resume | `--resume` skips PMIDs already in the output file and merges with it |
| Testing | `--limit N` processes the first N articles |

LLM errors and unparseable responses are appended to `failed_batches.jsonl`
(`pmid`, `mutations`, `error`, `timestamp`); that file is deleted at the
start of a run without `--resume`.

### Output

`verified_mutations.json` maps PMID to:

```json
{
  "verified_mutations":   [{"mutation": "...", "is_discussed": true,  "consequence_summary": "...", "evidence_sentences": [34, 35], "confidence": "high", "gene_match": null}],
  "unverified_mutations": [{"mutation": "...", "is_discussed": false, "consequence_summary": null,  "evidence_sentences": [],       "confidence": "low",  "gene_match": null}],
  "filtered_gene_names": ["..."],
  "just_mentioned": ["..."],
  "error": null
}
```

`evidence_sentences` are indices into the `split_into_sentences` output for
that article's full text.

Released run (`verification_stats.json`, model `google/gemini-2.5-flash-lite`):

| Field | Value |
|---|---|
| `total_articles` | 86,317 |
| `articles_with_verified_mutations` | 33,901 |
| `total_verified_mutations` | 234,148 |
| `total_unverified_by_llm` | 29,954 |
| `total_filtered_gene_names` | 150,859 |
| `total_just_mentioned` | 1,115,275 |
| `total_errors` | 0 |
| `failed_batches` | 1,147 |
