#!/usr/bin/env python3
"""
Fetch PMC full-text articles for humsavar-derived PMIDs.

Input: expansion/output/humsavar_targets.jsonl — one record per (gene, pmid) with curated variants.
Output: expansion/output/humsavar_articles.json — {pmid: {article, gene, pmc_id, variants}}
        (+ humsavar_articles_stats.json, humsavar_articles_checkpoint.json)
NCBI ID Converter + E-utilities (network); requires NCBI_EMAIL (env or repo-root .env).

Usage:
  python fetch_humsavar_articles.py [--limit N] [--resume] [--workers 8]
"""

import json
import time
import argparse
import logging
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import get_ncbi_api_key

# Reuse JATS XML parser from fetch_pmc_articles
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fetch_pmc_articles import extract_text_from_jats_xml, fetch_pmc_fulltext, get_ncbi_email

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)

EXPANSION_DIR = Path(__file__).resolve().parent
TARGETS = EXPANSION_DIR / "output" / "humsavar_targets.jsonl"
ARTICLES = EXPANSION_DIR / "output" / "humsavar_articles.json"
STATS = EXPANSION_DIR / "output" / "humsavar_articles_stats.json"
CHECKPOINT = EXPANSION_DIR / "output" / "humsavar_articles_checkpoint.json"

IDCONV_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"


def load_targets() -> dict[str, dict]:
    """Return {pmid: {gene, uniprot_acc, variants}}. If multiple genes share a pmid, merge."""
    by_pmid: dict[str, dict] = {}
    with open(TARGETS) as f:
        for line in f:
            r = json.loads(line)
            pmid = r["pmid"]
            if r.get("already_fetched"):
                continue
            if pmid not in by_pmid:
                by_pmid[pmid] = {
                    "gene": r["gene"],  # primary gene
                    "uniprot_acc": r["uniprot_acc"],
                    "additional_genes": [],
                    "variants": list(r["variants"]),
                }
            else:
                entry = by_pmid[pmid]
                if r["gene"] != entry["gene"] and r["gene"] not in entry["additional_genes"]:
                    entry["additional_genes"].append(r["gene"])
                entry["variants"].extend(r["variants"])
    return by_pmid


def pmid_to_pmcid_batch(pmids: list[str], api_key: str) -> dict[str, str]:
    """Convert PMIDs → PMCIDs using NCBI ID Converter. Returns {pmid: pmcid} (only for ones with PMC)."""
    mapping: dict[str, str] = {}
    batch = 200
    for i in range(0, len(pmids), batch):
        chunk = pmids[i:i+batch]
        params = {
            "ids": ",".join(chunk),
            "format": "json",
            "tool": "crossvalqa",
            }
        params["email"] = get_ncbi_email()
        if api_key:
            params["api_key"] = api_key

        for attempt in range(3):
            try:
                resp = requests.get(IDCONV_URL, params=params, timeout=30)
                if resp.status_code == 200:
                    data = resp.json()
                    for rec in data.get("records", []):
                        pmcid = rec.get("pmcid")
                        pmid = rec.get("pmid") or rec.get("requested-id")
                        if pmcid and pmid:
                            mapping[str(pmid)] = pmcid
                    break
                elif resp.status_code == 429:
                    time.sleep(2 * (attempt + 1))
                else:
                    break
            except Exception:
                time.sleep(1 * (attempt + 1))
        time.sleep(0.15)
    return mapping


def load_checkpoint() -> set[str]:
    if not CHECKPOINT.exists():
        return set()
    with open(CHECKPOINT) as f:
        return set(json.load(f).get("done_pmids", []))


def save_checkpoint(done: set[str]):
    tmp = CHECKPOINT.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump({"done_pmids": sorted(done)}, f)
    tmp.replace(CHECKPOINT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()

    ARTICLES.parent.mkdir(parents=True, exist_ok=True)
    get_ncbi_email()  # fail fast if NCBI_EMAIL is missing
    api_key = get_ncbi_api_key()
    if api_key:
        logger.info("Using NCBI API key")

    targets = load_targets()
    logger.info(f"Unique new PMIDs to process: {len(targets):,}")

    # Load existing articles + checkpoint
    articles: dict = {}
    if args.resume and ARTICLES.exists():
        with open(ARTICLES) as f:
            articles = json.load(f)
        logger.info(f"Resuming: {len(articles):,} articles already fetched")
    done = load_checkpoint() if args.resume else set()

    # Filter PMIDs to process
    pmids_to_process = [p for p in targets if p not in done and p not in articles]
    if args.limit:
        pmids_to_process = pmids_to_process[:args.limit]
    logger.info(f"PMIDs to process: {len(pmids_to_process):,}")

    # Step 1: batch-resolve PMID → PMCID
    logger.info("Step 1: PMID → PMCID batch conversion...")
    raw_pmc_map = pmid_to_pmcid_batch(pmids_to_process, api_key)
    # Filter to PMIDs we actually have in targets (API may normalize PMIDs)
    pmc_map = {p: pmcid for p, pmcid in raw_pmc_map.items() if p in targets}
    logger.info(f"  PMC availability: {len(pmc_map):,}/{len(pmids_to_process):,} ({100*len(pmc_map)/max(1,len(pmids_to_process)):.1f}%)")

    # Mark non-PMC PMIDs as done (nothing to fetch)
    for pmid in pmids_to_process:
        if pmid not in pmc_map:
            done.add(pmid)

    # Step 2: fetch full text for PMC-available PMIDs
    logger.info(f"Step 2: Fetching full text for {len(pmc_map):,} articles...")
    fetched = 0
    no_fulltext = 0
    errors = 0
    save_every = 200

    def fetch_one(pmid: str, pmcid: str) -> tuple[str, str | None]:
        pmc_num = pmcid.replace("PMC", "")
        text = fetch_pmc_fulltext(pmc_num, api_key)
        return pmid, text

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(fetch_one, pmid, pmcid): (pmid, pmcid)
                for pmid, pmcid in pmc_map.items()}

        for fut in as_completed(futs):
            pmid, pmcid = futs[fut]
            try:
                _, text = fut.result()
            except Exception as e:
                text = None
                errors += 1

            if text:
                entry = targets[pmid]
                articles[pmid] = {
                    "article": text,
                    "gene": entry["gene"],
                    "pmc_id": pmcid,
                    "uniprot_acc": entry["uniprot_acc"],
                    "additional_genes": entry["additional_genes"],
                    "variants": entry["variants"],
                }
                fetched += 1
            else:
                no_fulltext += 1

            done.add(pmid)
            total = fetched + no_fulltext + errors
            if total % save_every == 0:
                with open(ARTICLES, "w") as f:
                    json.dump(articles, f)
                save_checkpoint(done)
                logger.info(f"  [{total}/{len(pmc_map)}] fetched={fetched:,} no_fulltext={no_fulltext:,} errors={errors:,}")
            time.sleep(0.12)  # be polite

    # Final save
    with open(ARTICLES, "w") as f:
        json.dump(articles, f)
    save_checkpoint(done)

    stats = {
        "total_pmids": len(pmids_to_process),
        "pmc_available": len(pmc_map),
        "fulltext_fetched": fetched,
        "no_fulltext": no_fulltext,
        "errors": errors,
        "total_articles": len(articles),
    }
    with open(STATS, "w") as f:
        json.dump(stats, f, indent=2)

    logger.info(f"\n=== DONE ===")
    for k, v in stats.items():
        logger.info(f"  {k}: {v:,}")
    logger.info(f"Output: {ARTICLES}")


if __name__ == "__main__":
    main()
