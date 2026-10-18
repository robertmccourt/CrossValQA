"""
Central configuration for the CrossValQA claim extraction pipeline.

All paths are relative to PIPELINE_ROOT, which is auto-detected from this file's location.
Override any path by setting environment variables or editing this file.
"""

from pathlib import Path
import os

# ============================================================================
# ROOT PATHS
# ============================================================================

# Auto-detect: this file is at claim_extractor/config.py
CLAIM_EXTRACTOR_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = CLAIM_EXTRACTOR_DIR.parent

# ============================================================================
# DATABASE PATHS
# ============================================================================

DATABASE_DIR = PIPELINE_ROOT / "database"

# Phase 0 inputs/intermediates (see phase_0/build_article_set.py for the run order).
# SOURCE_ARTICLES_FILE is an external PubMed/PMC article collection built outside this repo.
SOURCE_ARTICLES_FILE = DATABASE_DIR / "pubmed_to_proteins_with_articles.backfilled.json"
EXTRACTED_ARTICLES_FILE = DATABASE_DIR / "pubmed_to_proteins_with_articles.mutations_extracted.json"
ARTICLES_FILE = DATABASE_DIR / "articles_with_mutations.json"
ARTICLES_NO_MUTATIONS_FILE = DATABASE_DIR / "articles_no_mutations_extracted.json"
VERIFIED_MUTATIONS_FILE = DATABASE_DIR / "verified_mutations.json"
EXTRACTION_STATS_FILE = DATABASE_DIR / "extraction_stats.json"
VERIFICATION_STATS_FILE = DATABASE_DIR / "verification_stats.json"
FAILED_BATCHES_FILE = DATABASE_DIR / "failed_batches.jsonl"
MISMATCHED_PMIDS_FILE = DATABASE_DIR / "mismatched_pmids.json"

# ============================================================================
# PHASE OUTPUT PATHS
# ============================================================================

PHASE0_OUTPUT_DIR = CLAIM_EXTRACTOR_DIR / "phase_0" / "output"
PHASE1_OUTPUT_DIR = CLAIM_EXTRACTOR_DIR / "phase_1" / "output"
PHASE2_OUTPUT_DIR = CLAIM_EXTRACTOR_DIR / "phase_2" / "output"
PHASE3_OUTPUT_DIR = CLAIM_EXTRACTOR_DIR / "phase_3" / "output"

PHASE1_QA_FILE = PHASE1_OUTPUT_DIR / "phase1_qa_pairs.jsonl"
PHASE2_QA_FILE = PHASE2_OUTPUT_DIR / "phase2_qa_pairs.jsonl"
PHASE3_VALIDATED_FILE = PHASE3_OUTPUT_DIR / "phase3_validated.jsonl"
PHASE3_RECOVERY_DIR = PHASE3_OUTPUT_DIR / "recovery"
PHASE3_RECOVERED_FILE = PHASE3_RECOVERY_DIR / "phase3_recovered.jsonl"
PHASE3_MERGED_FILE = PHASE3_OUTPUT_DIR / "phase3_validated_merged.jsonl"

# Phase 4 (gene recovery)
PHASE4_DIR = CLAIM_EXTRACTOR_DIR / "phase_4"
PHASE4_OUTPUT_DIR = PHASE4_DIR / "output"
PHASE4_FINAL_FILE = PHASE4_OUTPUT_DIR / "m2tqa_final.jsonl"
PHASE4_WITH_CITED_FILE = PHASE4_OUTPUT_DIR / "m2tqa_with_cited.jsonl"
PHASE4_GENE_RECOVERED_FILE = PHASE4_OUTPUT_DIR / "m2tqa_gene_recovered.jsonl"

# ============================================================================
# REFERENCE DATA
# ============================================================================

MAPPINGS_DIR = PIPELINE_ROOT / "mappings"
GENE_INFO_FILE = MAPPINGS_DIR / "Homo_sapiens.gene_info"
HGNC_GENES_FILE = PHASE4_DIR / "hgnc_genes.txt"

# ============================================================================
# UNIPROT GROUNDING (final pipeline output)
# ============================================================================
# This is the last file produced by this pipeline. The released MutQA splits are
# built from it downstream (release step) and published on Hugging Face.

UNIPROT_DIR = PIPELINE_ROOT / "uniprot"
UNIPROT_INPUT_DIR = UNIPROT_DIR / "input"
UNIPROT_OUTPUT_DIR = UNIPROT_DIR / "output"
UNIPROT_GROUNDED_FILE = UNIPROT_OUTPUT_DIR / "m2tqa_grounded_final.jsonl"
UNIPROT_GROUNDING_STATS_FILE = UNIPROT_OUTPUT_DIR / "grounding_stats_final.json"
UNIPROT_GENE_CORRECTED_FILE = UNIPROT_OUTPUT_DIR / "m2tqa_gene_corrected.jsonl"
PIPELINE_OUTPUT_FILE = UNIPROT_GENE_CORRECTED_FILE

# ============================================================================
# LLM DEFAULTS
# ============================================================================

DEFAULT_PHASE0_VERIFY_MODEL = "google/gemini-2.5-flash-lite"  # model used for the released run
DEFAULT_PHASE1_MODEL = "google/gemini-2.5-flash"
DEFAULT_PHASE2_MODEL = "deepseek/deepseek-chat"
DEFAULT_PHASE3_MODEL = "meta-llama/llama-3.1-8b-instruct"

# ============================================================================
# API
# ============================================================================

def get_openrouter_api_key() -> str:
    """Get OpenRouter API key from environment."""
    key = os.environ.get("OPENROUTER_API_KEY", "")
    if not key:
        from dotenv import load_dotenv
        load_dotenv(PIPELINE_ROOT / ".env")
        key = os.environ.get("OPENROUTER_API_KEY", "")
    if not key:
        raise ValueError("OPENROUTER_API_KEY not set. Add to .env or export it.")
    return key


def get_ncbi_api_key() -> str:
    """Get NCBI API key from environment."""
    key = os.environ.get("NCBI_API_KEY", "")
    if not key:
        from dotenv import load_dotenv
        load_dotenv(PIPELINE_ROOT / ".env")
        key = os.environ.get("NCBI_API_KEY", "")
    return key  # optional, not required
