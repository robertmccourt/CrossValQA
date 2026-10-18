#!/usr/bin/env python3
"""
Expansion Pipeline Runner

Runs Phases 0-3 on one expansion article batch (batch1, batch2, batch2_round2,
reprocess_filtered, humsavar, or any --articles file).

Phase 0 is handled internally (regex extraction + gene-aware LLM verification,
keeping only gene_match=True mutations). Phases 1-3 are called as subprocesses
using the existing phase scripts with path overrides. Makes OpenRouter calls
(OPENROUTER_API_KEY) — run a --limit pilot first.

Input:  expansion/output/<batch articles>.json  ({pmid: {article, gene, ...}})
Output: expansion/output/<batch>/{extracted_mutations,verified_mutations}.json,
        <batch>/phase{1,2,3}/..., notably <batch>/phase3/phase3_validated.jsonl

Usage:
    # Full run on Batch 1
    python run_expansion.py --batch batch1

    # Full run on Batch 2
    python run_expansion.py --batch batch2

    # Phase 0 only (extraction + verification)
    python run_expansion.py --batch batch1 --phase 0

    # Phases 1-3 only (assumes Phase 0 already done)
    python run_expansion.py --batch batch1 --phase 1 2 3

    # Pilot (10 articles)
    python run_expansion.py --batch batch1 --limit 10

    # Resume interrupted run
    python run_expansion.py --batch batch1 --resume
"""

import json
import os
import sys
import argparse
import subprocess
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# Setup imports
EXPANSION_DIR = Path(__file__).resolve().parent
CLAIM_EXTRACTOR_DIR = EXPANSION_DIR.parent
PIPELINE_ROOT = CLAIM_EXTRACTOR_DIR.parent
sys.path.insert(0, str(CLAIM_EXTRACTOR_DIR))
sys.path.insert(0, str(CLAIM_EXTRACTOR_DIR / "phase_0"))

# Load .env early so subprocesses inherit API keys
from dotenv import load_dotenv
load_dotenv(PIPELINE_ROOT / ".env")

from config import get_openrouter_api_key, DEFAULT_PHASE1_MODEL, DEFAULT_PHASE2_MODEL, DEFAULT_PHASE3_MODEL


# ============================================================================
# PATHS
# ============================================================================

def get_batch_paths(batch_name: str, articles_override: str | None = None) -> dict:
    """Get all input/output paths for a batch."""
    output_dir = EXPANSION_DIR / "output" / batch_name
    output_dir.mkdir(parents=True, exist_ok=True)

    if articles_override:
        articles_file = Path(articles_override)
    elif batch_name == "batch1":
        articles_file = EXPANSION_DIR / "output" / "expansion_articles.json"
    elif batch_name == "batch2":
        articles_file = EXPANSION_DIR / "output" / "batch2_articles.json"
    elif batch_name == "batch2_round2":
        articles_file = EXPANSION_DIR / "output" / "batch2_round2_articles.json"
    elif batch_name == "reprocess":
        articles_file = EXPANSION_DIR / "output" / "reprocess_articles.json"
    elif batch_name == "reprocess_filtered":
        articles_file = EXPANSION_DIR / "output" / "reprocess_filtered.json"
    else:
        raise ValueError(f"Unknown batch: {batch_name}. Use --articles to override.")

    return {
        "articles": str(articles_file),
        "output_dir": str(output_dir),
        # Phase 0 outputs
        "extracted": str(output_dir / "extracted_mutations.json"),
        "verified": str(output_dir / "verified_mutations.json"),
        "verification_stats": str(output_dir / "verification_stats.json"),
        # Phase 1 outputs
        "phase1_dir": str(output_dir / "phase1"),
        "phase1_qa": str(output_dir / "phase1" / "phase1_qa_pairs.jsonl"),
        # Phase 2 outputs
        "phase2_dir": str(output_dir / "phase2"),
        "phase2_qa": str(output_dir / "phase2" / "phase2_qa_pairs.jsonl"),
        # Phase 3 outputs
        "phase3_dir": str(output_dir / "phase3"),
        "phase3_validated": str(output_dir / "phase3" / "phase3_validated.jsonl"),
    }


# ============================================================================
# PHASE 0: EXTRACTION + GENE-AWARE VERIFICATION
# ============================================================================

def _extract_one(args):
    """Extract mutations from a single article (module-level for pickling)."""
    from mutation_regex import extract_all_variants, get_mutations_for_llm
    pmid, entry = args
    text = entry.get("article", "")
    if not text:
        return pmid, [], 0
    try:
        results = extract_all_variants(text, apply_section_filter=True)
        llm_data = get_mutations_for_llm(results)
        mutations = llm_data["mutations"]
        return pmid, mutations, len(mutations)
    except Exception:
        return pmid, [], 0


def run_extraction(articles_path: str, output_path: str, limit: int = 0):
    """Run regex mutation extraction on expansion articles."""
    print(f"Loading articles from {articles_path}...")
    with open(articles_path) as f:
        articles = json.load(f)

    if limit > 0:
        articles = dict(list(articles.items())[:limit])

    print(f"Extracting mutations from {len(articles):,} articles...")
    start = time.time()

    total_mutations = 0
    workers = min(os.cpu_count() or 4, len(articles))
    work_items = list(articles.items())

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_extract_one, item): item[0] for item in work_items}
        for future in as_completed(futures):
            pmid, mutations, count = future.result()
            articles[pmid]["proteins"] = mutations
            total_mutations += count

    elapsed = time.time() - start
    print(f"Extracted {total_mutations:,} mutations from {len(articles):,} articles in {elapsed:.1f}s")

    # Save
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(articles, f)
    print(f"Saved to {output_path}")
    return articles


def run_verification(articles_path: str, output_path: str, stats_path: str,
                     model: str, workers: int, resume: bool):
    """Run gene-aware LLM verification."""
    api_key = get_openrouter_api_key()
    cmd = [
        sys.executable, str(CLAIM_EXTRACTOR_DIR / "phase_0" / "verify_mutations.py"),
        "--input", articles_path,
        "--output", output_path,
        "--stats-output", stats_path,
        "--model", model,
        "--workers", str(workers),
        "--api-key", api_key,
    ]
    if resume:
        cmd.append("--resume")

    print(f"\nRunning verification: --model {model} --workers {workers}")
    result = subprocess.run(cmd, cwd=str(CLAIM_EXTRACTOR_DIR))
    if result.returncode != 0:
        print(f"ERROR: Verification failed with code {result.returncode}")
        sys.exit(1)


def filter_gene_matches(verified_path: str):
    """Filter verified mutations to only gene_match=True entries."""
    with open(verified_path) as f:
        data = json.load(f)

    total_verified = 0
    total_gene_match = 0
    total_gene_mismatch = 0
    total_no_gene_field = 0

    for pmid, content in data.items():
        original = content.get("verified_mutations", [])
        filtered = []
        for mut in original:
            total_verified += 1
            gm = mut.get("gene_match")
            if gm is True:
                filtered.append(mut)
                total_gene_match += 1
            elif gm is False:
                total_gene_mismatch += 1
            else:
                # gene_match=None means original pipeline (no target gene) — keep
                filtered.append(mut)
                total_no_gene_field += 1
        content["verified_mutations"] = filtered

    print(f"\n=== GENE MATCH FILTER ===")
    print(f"Total verified mutations: {total_verified:,}")
    print(f"Gene match (kept):       {total_gene_match:,}")
    print(f"Gene mismatch (removed): {total_gene_mismatch:,}")
    print(f"No gene field (kept):    {total_no_gene_field:,}")

    with open(verified_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Updated {verified_path}")


# ============================================================================
# PHASES 1-3: SUBPROCESS CALLS
# ============================================================================

def run_phase1(paths: dict, model: str, workers: int, resume: bool, limit: int = 0):
    """Run Phase 1 QA generation."""
    Path(paths["phase1_dir"]).mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(CLAIM_EXTRACTOR_DIR / "phase_1" / "phase1_generator_qa.py"),
        "--model", model,
        "--workers", str(workers),
        "--input-verified", paths["verified"],
        "--input-articles", paths["articles"],
        "--output-dir", paths["phase1_dir"],
    ]
    if resume:
        cmd.append("--resume")
    if limit > 0:
        cmd.extend(["--limit", str(limit)])

    print(f"\n{'='*60}")
    print(f"PHASE 1: QA Generation")
    print(f"{'='*60}")
    result = subprocess.run(cmd, cwd=str(CLAIM_EXTRACTOR_DIR))
    if result.returncode != 0:
        print(f"ERROR: Phase 1 failed with code {result.returncode}")
        sys.exit(1)


def run_phase2(paths: dict, model: str, workers: int, resume: bool, limit: int = 0):
    """Run Phase 2 independent answer generation."""
    Path(paths["phase2_dir"]).mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(CLAIM_EXTRACTOR_DIR / "phase_2" / "phase2_generator.py"),
        "--model", model,
        "--workers", str(workers),
        "--input-phase1", paths["phase1_qa"],
        "--input-articles", paths["articles"],
        "--output-dir", paths["phase2_dir"],
    ]
    if resume:
        cmd.append("--resume")
    if limit > 0:
        cmd.extend(["--limit", str(limit)])

    print(f"\n{'='*60}")
    print(f"PHASE 2: Independent Answers")
    print(f"{'='*60}")
    result = subprocess.run(cmd, cwd=str(CLAIM_EXTRACTOR_DIR))
    if result.returncode != 0:
        print(f"ERROR: Phase 2 failed with code {result.returncode}")
        sys.exit(1)


def run_phase3(paths: dict, model: str, workers: int, resume: bool):
    """Run Phase 3 cross-validation."""
    Path(paths["phase3_dir"]).mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(CLAIM_EXTRACTOR_DIR / "phase_3" / "phase3_judge.py"),
        "--model", model,
        "--workers", str(workers),
        "--input-phase2", paths["phase2_qa"],
        "--output-dir", paths["phase3_dir"],
        "--once",  # Don't poll, just process and exit
    ]
    if resume:
        cmd.append("--resume")

    print(f"\n{'='*60}")
    print(f"PHASE 3: Cross-Validation")
    print(f"{'='*60}")
    result = subprocess.run(cmd, cwd=str(CLAIM_EXTRACTOR_DIR))
    if result.returncode != 0:
        print(f"ERROR: Phase 3 failed with code {result.returncode}")
        sys.exit(1)


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Run expansion pipeline (Phases 0-3)")
    parser.add_argument("--articles", default=None,
                        help="Override path to articles JSON (enables arbitrary batch names)")
    parser.add_argument("--batch", required=True,
                        help="Which batch to process")
    parser.add_argument("--phase", nargs="*", type=int, default=None,
                        help="Run specific phases (default: all 0-3)")
    parser.add_argument("--limit", type=int, default=0,
                        help="Limit articles/mutations (for testing)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume interrupted run")
    parser.add_argument("--workers", type=int, default=50,
                        help="Max concurrent API requests")
    parser.add_argument("--verification-model", default="google/gemini-2.5-flash",
                        help="Model for Phase 0 verification")
    parser.add_argument("--phase1-model", default=DEFAULT_PHASE1_MODEL,
                        help="Model for Phase 1 QA generation")
    parser.add_argument("--phase2-model", default=DEFAULT_PHASE2_MODEL,
                        help="Model for Phase 2 answers")
    parser.add_argument("--phase3-model", default=DEFAULT_PHASE3_MODEL,
                        help="Model for Phase 3 judging")
    args = parser.parse_args()

    phases = set(args.phase) if args.phase is not None else {0, 1, 2, 3}
    paths = get_batch_paths(args.batch, articles_override=args.articles)

    print(f"=== EXPANSION PIPELINE: {args.batch} ===")
    print(f"Articles: {paths['articles']}")
    print(f"Output:   {paths['output_dir']}")
    print(f"Phases:   {sorted(phases)}")
    if args.limit:
        print(f"Limit:    {args.limit}")
    print()

    # Verify articles exist
    if not Path(paths["articles"]).exists():
        print(f"ERROR: Articles file not found: {paths['articles']}")
        sys.exit(1)

    # Phase 0: Extraction + Verification
    if 0 in phases:
        print(f"{'='*60}")
        print(f"PHASE 0: Extraction + Gene-Aware Verification")
        print(f"{'='*60}")

        # Step 1: Regex extraction
        if not args.resume or not Path(paths["extracted"]).exists():
            run_extraction(paths["articles"], paths["extracted"], args.limit)
        else:
            print(f"Skipping extraction (resume mode, {paths['extracted']} exists)")

        # Step 2: LLM verification with gene_match
        run_verification(
            paths["extracted"], paths["verified"], paths["verification_stats"],
            args.verification_model, args.workers, args.resume,
        )

        # Step 3: Filter to gene_match=True only
        filter_gene_matches(paths["verified"])

    # Phase 1: QA Generation
    if 1 in phases:
        run_phase1(paths, args.phase1_model, args.workers, args.resume, args.limit)

    # Phase 2: Independent Answers
    if 2 in phases:
        run_phase2(paths, args.phase2_model, args.workers, args.resume, args.limit)

    # Phase 3: Cross-Validation
    if 3 in phases:
        run_phase3(paths, args.phase3_model, args.workers, args.resume)

    # Summary
    print(f"\n{'='*60}")
    print(f"EXPANSION PIPELINE COMPLETE: {args.batch}")
    print(f"{'='*60}")
    if Path(paths["phase3_validated"]).exists():
        count = sum(1 for _ in open(paths["phase3_validated"]))
        print(f"Phase 3 output: {count:,} validated QA pairs")
    print(f"Output directory: {paths['output_dir']}")


if __name__ == "__main__":
    main()
