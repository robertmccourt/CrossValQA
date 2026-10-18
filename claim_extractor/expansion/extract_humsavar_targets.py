#!/usr/bin/env python3
"""
Humsavar targets: (gene, UniProt, variant, PMID) for genes missing from the
pre-expansion MutQA dataset, from UniProt humsavar + the EBI Proteins API.

Input:
  - uniprot/output/m2tqa_gene_corrected.jsonl (pre-expansion dataset; genes to skip)
  - expansion/data/humsavar.txt (UniProt curated human missense variants)
  - EBI Proteins API /variation (per-variant PubMed evidence; network)
  - optional expansion/output/exclude_pmids.json (JSON list of PMIDs fetched in
    earlier batches; matching targets are flagged already_fetched)

Output:
  - expansion/output/humsavar_targets.jsonl: one record per (gene, pmid), with list of
    (variant, ACMG category, disease) tuples attached
  - expansion/output/humsavar_checkpoint.json (resume state)

Usage:
  python extract_humsavar_targets.py [--limit N] [--workers 10] [--resume]
"""

import json
import time
import argparse
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

import requests

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import UNIPROT_GENE_CORRECTED_FILE

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)

EXPANSION_DIR = Path(__file__).resolve().parent
HUMSAVAR = EXPANSION_DIR / "data" / "humsavar.txt"
OUTPUT = EXPANSION_DIR / "output" / "humsavar_targets.jsonl"
CHECKPOINT = EXPANSION_DIR / "output" / "humsavar_checkpoint.json"
EXCLUDE_PMIDS_FILE = EXPANSION_DIR / "output" / "exclude_pmids.json"

UNIPROT_API = "https://www.ebi.ac.uk/proteins/api/variation"


def load_current_genes() -> set[str]:
    genes = set()
    with open(UNIPROT_GENE_CORRECTED_FILE) as f:
        for line in f:
            r = json.loads(line)
            g = r.get("corrected_gene") or r.get("recovered_gene") or r.get("anchor_gene_name", "")
            if g:
                genes.add(g)
    return genes


def parse_humsavar() -> dict[str, dict]:
    """Return {gene: {acc: ..., variants: [(featureId, aa_change, category, dbsnp, disease)]}}"""
    data: dict[str, dict] = {}
    with open(HUMSAVAR) as f:
        in_data = False
        for line in f:
            if line.startswith("_________"):
                in_data = True
                continue
            if not in_data:
                continue
            parts = line.split(None, 6)
            if len(parts) < 6:
                continue
            gene, acc, ftid, aa_change, category, dbsnp = parts[:6]
            disease = parts[6].strip() if len(parts) > 6 else ""
            if gene not in data:
                data[gene] = {"acc": acc, "variants": []}
            data[gene]["variants"].append({
                "feature_id": ftid,
                "aa_change": aa_change,
                "acmg": category,
                "dbsnp": dbsnp if dbsnp != "-" else "",
                "disease": disease if disease != "-" else "",
            })
    return data


def fetch_uniprot_variants(acc: str, retries: int = 3) -> list[dict]:
    """Fetch per-variant PubMed evidences from UniProt. Returns list of features (UniProt-curated only)."""
    url = f"{UNIPROT_API}/{acc}"
    for attempt in range(retries):
        try:
            r = requests.get(url, headers={"Accept": "application/json"}, timeout=20)
            if r.status_code == 200:
                data = r.json()
                return [f for f in data.get("features", []) if f.get("sourceType") in ("uniprot", "mixed")]
            elif r.status_code == 429:
                time.sleep(2 * (attempt + 1))
            else:
                return []
        except Exception:
            time.sleep(1 * (attempt + 1))
    return []


def extract_gene_targets(gene: str, acc: str, humsavar_vars: list[dict]) -> dict:
    """Fetch UniProt variants for a gene and produce target records keyed by PMID."""
    features = fetch_uniprot_variants(acc)

    # Build humsavar lookup by feature_id (VAR_xxxxxx)
    hs_by_ftid = {v["feature_id"]: v for v in humsavar_vars}

    # For each UniProt-curated variant, pull PubMed IDs, match to humsavar metadata
    pmid_to_vars: dict[str, list[dict]] = defaultdict(list)
    for f in features:
        ftid = None
        for xref in f.get("xrefs", []):
            if xref.get("name") == "UniProt" and xref.get("id", "").startswith("VAR_"):
                ftid = xref["id"]
                break

        hs_meta = hs_by_ftid.get(ftid, {})
        # Fallback: use the API's own aa change
        aa_change = hs_meta.get("aa_change") or f"p.{f.get('wildType','')}{f.get('begin','')}{f.get('mutatedType','')}"

        pmids = {ev["source"]["id"] for ev in f.get("evidences", [])
                 if ev.get("source", {}).get("name") == "pubmed"}

        for pmid in pmids:
            pmid_to_vars[pmid].append({
                "variant": aa_change,
                "feature_id": ftid or "",
                "acmg": hs_meta.get("acmg", ""),
                "dbsnp": hs_meta.get("dbsnp", ""),
                "disease": hs_meta.get("disease", ""),
            })

    return {"gene": gene, "uniprot_acc": acc, "pmid_to_vars": dict(pmid_to_vars)}


def load_checkpoint() -> set[str]:
    if not CHECKPOINT.exists():
        return set()
    with open(CHECKPOINT) as f:
        return set(json.load(f).get("done_genes", []))


def save_checkpoint(done: set[str]):
    tmp = CHECKPOINT.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump({"done_genes": sorted(done)}, f)
    tmp.replace(CHECKPOINT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="Limit genes (0=all)")
    ap.add_argument("--workers", type=int, default=8, help="Concurrent UniProt requests")
    ap.add_argument("--resume", action="store_true", help="Resume from checkpoint")
    ap.add_argument("--include-benign", action="store_true",
                    help="Include genes with only LB/B variants (default: skip)")
    args = ap.parse_args()

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    current_genes = load_current_genes()
    logger.info(f"Current genes: {len(current_genes):,}")

    # Optional exclude: PMIDs already fetched in prior batches
    exclude_pmids: set[str] = set()
    if EXCLUDE_PMIDS_FILE.exists():
        with open(EXCLUDE_PMIDS_FILE) as f:
            exclude_pmids = set(json.load(f))
        logger.info(f"Excluding {len(exclude_pmids):,} previously fetched PMIDs")

    humsavar = parse_humsavar()
    missing = {g: d for g, d in humsavar.items() if g not in current_genes}
    logger.info(f"Humsavar missing genes: {len(missing):,}")

    # Filter to genes with at least one LP/P or US variant (skip benign-only)
    if not args.include_benign:
        filtered = {}
        for g, d in missing.items():
            if any(v["acmg"] in ("LP/P", "US") for v in d["variants"]):
                filtered[g] = d
        logger.info(f"After LP/P+US filter: {len(filtered):,} genes (dropped {len(missing)-len(filtered):,} benign-only)")
        missing = filtered

    done = load_checkpoint() if args.resume else set()
    if done:
        logger.info(f"Resuming: {len(done):,} genes already processed")

    todo = sorted(g for g in missing if g not in done)
    if args.limit:
        todo = todo[:args.limit]
    logger.info(f"To process: {len(todo):,}")

    # Append mode for output (compatible with resume)
    out_mode = "a" if args.resume and OUTPUT.exists() else "w"
    total_records = 0
    total_pmids = set()
    new_pmids = set()
    processed = 0

    with open(OUTPUT, out_mode) as out, \
         ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(extract_gene_targets, g, missing[g]["acc"], missing[g]["variants"]): g
                   for g in todo}

        for fut in as_completed(futures):
            gene = futures[fut]
            try:
                result = fut.result()
            except Exception as e:
                logger.warning(f"Failed {gene}: {e}")
                done.add(gene)
                continue

            for pmid, vlist in result["pmid_to_vars"].items():
                rec = {
                    "gene": result["gene"],
                    "uniprot_acc": result["uniprot_acc"],
                    "pmid": pmid,
                    "variants": vlist,
                    "already_fetched": pmid in exclude_pmids,
                }
                out.write(json.dumps(rec) + "\n")
                total_records += 1
                total_pmids.add(pmid)
                if pmid not in exclude_pmids:
                    new_pmids.add(pmid)

            done.add(gene)
            processed += 1
            if processed % 200 == 0:
                save_checkpoint(done)
                logger.info(f"  [{processed}/{len(todo)}] records={total_records:,} unique_pmids={len(total_pmids):,} new={len(new_pmids):,}")

    save_checkpoint(done)
    logger.info(f"\n=== DONE ===")
    logger.info(f"Processed genes: {processed:,}")
    logger.info(f"Records written: {total_records:,}")
    logger.info(f"Unique PMIDs: {len(total_pmids):,}")
    logger.info(f"New PMIDs (not in exclude): {len(new_pmids):,}")
    logger.info(f"Output: {OUTPUT}")


if __name__ == "__main__":
    main()
