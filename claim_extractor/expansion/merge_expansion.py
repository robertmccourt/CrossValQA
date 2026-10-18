#!/usr/bin/env python3
"""
Merge expansion phase-3 outputs into the main MutQA dataset with UniProt identifiers.

For each expansion batch, we:
  1. Load <batch>/phase3/phase3_validated.jsonl
  2. Join to the per-batch articles file by pmid to recover the target gene
  3. Map gene symbol -> NCBI gene_id -> UniProt acc (human 9606)
  4. For humsavar, also attach the UniProt acc + curated variant metadata
     (feature_id, dbSNP, ACMG category, disease) from humsavar_targets.jsonl
Then dedupe on (pmid, mutation, question), drop records already in the main
dataset, and append the novel ones after the unchanged main records.

Input:  uniprot/output/m2tqa_gene_corrected.jsonl (416,188 rows, question types normalized)
        expansion/output/{batch1,batch2,batch2_round2,reprocess_filtered,humsavar}/phase3/
        expansion/output/*_articles.json, humsavar_targets.jsonl
        mappings/Homo_sapiens.gene_info, uniprot/input/HUMAN_9606_idmapping.dat.gz
Output: uniprot/output/m2tqa_gene_corrected_v2.jsonl (423,415 rows)
        uniprot/output/expansion_merge_stats.json

Usage:
    python merge_expansion.py
    python merge_expansion.py --main IN.jsonl --output OUT.jsonl --stats STATS.json
"""

import argparse
import csv
import gzip
import json
import sys
from collections import defaultdict, Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import (
    GENE_INFO_FILE, UNIPROT_INPUT_DIR, UNIPROT_OUTPUT_DIR, UNIPROT_GENE_CORRECTED_FILE,
)

EXP_OUT = Path(__file__).resolve().parent / "output"
NCBI_GENE_INFO = GENE_INFO_FILE
UNIPROT_HUMAN = UNIPROT_INPUT_DIR / "HUMAN_9606_idmapping.dat.gz"
MAIN_DATASET = UNIPROT_GENE_CORRECTED_FILE
MERGED_OUT = UNIPROT_OUTPUT_DIR / "m2tqa_gene_corrected_v2.jsonl"
STATS_OUT = UNIPROT_OUTPUT_DIR / "expansion_merge_stats.json"

BATCHES = [
    ("batch1",             EXP_OUT / "expansion_articles.json"),
    ("batch2",             EXP_OUT / "batch2_articles.json"),
    ("batch2_round2",      EXP_OUT / "batch2_round2_articles.json"),
    ("reprocess_filtered", EXP_OUT / "reprocess_filtered.json"),
    ("humsavar",           EXP_OUT / "humsavar_articles.json"),
]


def load_ncbi_symbol_to_geneid():
    symbol_to_id: dict[str, tuple[str, str]] = {}
    with open(NCBI_GENE_INFO) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid = row["GeneID"]
            official = row["Symbol"]
            symbol_to_id[official.upper()] = (gid, official)
            syns = row["Synonyms"].split("|") if row["Synonyms"] != "-" else []
            for syn in syns:
                if syn and syn.upper() not in symbol_to_id:
                    symbol_to_id[syn.upper()] = (gid, official)
    return symbol_to_id


def load_uniprot_human():
    gene_to_uniprot: dict[str, list[str]] = defaultdict(list)
    with gzip.open(UNIPROT_HUMAN, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 3 and parts[1] == "GeneID":
                gene_to_uniprot[parts[2]].append(parts[0])
    return dict(gene_to_uniprot)


def load_humsavar_pmid_variants():
    """pmid -> list of {feature_id, variant, acmg, dbsnp, disease}, plus pmid -> uniprot_acc."""
    path = EXP_OUT / "humsavar_targets.jsonl"
    pmid_to_variants: dict[str, list[dict]] = defaultdict(list)
    pmid_to_uniprot: dict[str, str] = {}
    if not path.exists():
        return pmid_to_variants, pmid_to_uniprot
    with open(path) as f:
        for line in f:
            r = json.loads(line)
            pmid = r["pmid"]
            pmid_to_uniprot.setdefault(pmid, r.get("uniprot_acc", ""))
            pmid_to_variants[pmid].extend(r.get("variants", []))
    return pmid_to_variants, pmid_to_uniprot


def process_batch(batch_name: str, articles_path: Path,
                  symbol_to_id: dict, human_gene_to_uniprot: dict,
                  humsavar_vars: dict, humsavar_uniprot: dict) -> tuple[list[dict], Counter]:
    phase3_file = EXP_OUT / batch_name / "phase3/phase3_validated.jsonl"
    if not phase3_file.exists():
        print(f"[{batch_name}] no phase3 file, skipping")
        return [], Counter()

    with open(articles_path) as f:
        arts = json.load(f)
    pmid_to_gene = {pmid: info.get("gene") for pmid, info in arts.items()}
    pmid_to_addl = {pmid: info.get("additional_genes", []) for pmid, info in arts.items()}

    out: list[dict] = []
    stats: Counter = Counter()

    with open(phase3_file) as f:
        for line in f:
            r = json.loads(line)
            stats["records"] += 1
            pmid = r["pmid"]
            gene = pmid_to_gene.get(pmid)
            if not gene:
                stats["no_gene"] += 1
                continue

            # Symbol -> gene_id
            gid_pair = symbol_to_id.get(gene.upper())
            if gid_pair:
                gene_id, official = gid_pair
                uniprots = human_gene_to_uniprot.get(gene_id, [])
            else:
                gene_id, official, uniprots = "", gene, []

            rec = dict(r)
            rec["gene"] = gene
            rec["anchor_gene_name"] = official
            rec["anchor_gene_id"] = gene_id
            rec["anchor_tier"] = "expansion_target"
            rec["anchor_confidence"] = "high"
            rec["ncbi_taxon_id"] = "9606"
            rec["corrected_gene"] = official
            rec["gene_confidence"] = "gold" if gid_pair else "low"
            rec["gene_correction_method"] = "expansion_direct"
            rec["expansion_batch"] = batch_name
            if pmid_to_addl.get(pmid):
                rec["additional_genes"] = pmid_to_addl[pmid]

            # UniProt accession
            if batch_name == "humsavar" and humsavar_uniprot.get(pmid):
                primary = humsavar_uniprot[pmid]
                all_ids = [primary] + [u for u in uniprots if u != primary]
                rec["uniprot_ids"] = all_ids
                rec["uniprot_primary"] = primary
                rec["grounding_method"] = "humsavar_curated"
                stats["uniprot_humsavar"] += 1
            elif uniprots:
                rec["uniprot_ids"] = uniprots
                rec["uniprot_primary"] = uniprots[0]
                rec["grounding_method"] = "expansion_symbol_lookup"
                stats["uniprot_symbol"] += 1
            else:
                rec["uniprot_ids"] = []
                rec["uniprot_primary"] = ""
                rec["grounding_method"] = "ungrounded"
                stats["uniprot_missing"] += 1

            # Humsavar curated variant metadata (match by consequence if possible)
            if batch_name == "humsavar" and humsavar_vars.get(pmid):
                mutation = rec.get("mutation", "")
                matches = [v for v in humsavar_vars[pmid]
                           if v.get("variant") and (
                               v["variant"] == mutation or
                               v["variant"].lstrip("p.") == mutation.lstrip("p.")
                           )]
                rec["humsavar_variants_all"] = humsavar_vars[pmid]
                if matches:
                    m = matches[0]
                    rec["humsavar_feature_id"] = m.get("feature_id", "")
                    rec["humsavar_acmg"] = m.get("acmg", "")
                    rec["humsavar_dbsnp"] = m.get("dbsnp", "")
                    rec["humsavar_disease"] = m.get("disease", "")
                    stats["humsavar_matched"] += 1

            out.append(rec)
            stats["kept"] += 1

    print(f"[{batch_name}] records={stats['records']} kept={stats['kept']} "
          f"no_gene={stats['no_gene']} uniprot_missing={stats['uniprot_missing']}")
    return out, stats


def main():
    global MAIN_DATASET, MERGED_OUT, STATS_OUT
    ap = argparse.ArgumentParser(description="Merge expansion batches into the main dataset (-> v2).")
    ap.add_argument("--main", type=Path, default=MAIN_DATASET, help="pre-expansion dataset (default: uniprot/output/m2tqa_gene_corrected.jsonl)")
    ap.add_argument("--output", type=Path, default=MERGED_OUT, help="merged output (default: uniprot/output/m2tqa_gene_corrected_v2.jsonl)")
    ap.add_argument("--stats", type=Path, default=STATS_OUT, help="merge stats JSON (default: uniprot/output/expansion_merge_stats.json)")
    args = ap.parse_args()
    MAIN_DATASET, MERGED_OUT, STATS_OUT = args.main, args.output, args.stats

    print("Loading NCBI gene info...")
    symbol_to_id = load_ncbi_symbol_to_geneid()
    print(f"  {len(symbol_to_id):,} symbols/aliases")

    print("Loading UniProt human mappings...")
    human_gene_to_uniprot = load_uniprot_human()
    print(f"  {len(human_gene_to_uniprot):,} gene IDs with UniProt")

    print("Loading humsavar variant metadata...")
    humsavar_vars, humsavar_uniprot = load_humsavar_pmid_variants()
    print(f"  {len(humsavar_vars):,} PMIDs with curated variants")

    all_expansion: list[dict] = []
    per_batch_stats: dict[str, dict] = {}
    for batch_name, path in BATCHES:
        if not path.exists():
            print(f"[{batch_name}] articles file missing: {path}")
            continue
        recs, stats = process_batch(batch_name, path, symbol_to_id,
                                     human_gene_to_uniprot, humsavar_vars, humsavar_uniprot)
        all_expansion.extend(recs)
        per_batch_stats[batch_name] = dict(stats)

    print(f"\nTotal expansion records ready: {len(all_expansion):,}")

    # Dedupe: (pmid, mutation, question) — same QA generated in multiple batches is rare but possible
    seen: set[tuple] = set()
    deduped = []
    for r in all_expansion:
        key = (r["pmid"], r["mutation"], r["question"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(r)
    dupes_removed = len(all_expansion) - len(deduped)
    print(f"Duplicate QA pairs removed: {dupes_removed:,}")

    # Collision check: do any expansion (pmid, mutation, question) already exist in main dataset?
    print("\nChecking collisions with main dataset...")
    main_keys = set()
    with open(MAIN_DATASET) as f:
        for line in f:
            r = json.loads(line)
            main_keys.add((r["pmid"], r["mutation"], r["question"]))
    collisions = sum(1 for r in deduped if (r["pmid"], r["mutation"], r["question"]) in main_keys)
    print(f"  Main dataset records: {len(main_keys):,}")
    print(f"  Expansion records already present: {collisions:,}")
    novel = [r for r in deduped if (r["pmid"], r["mutation"], r["question"]) not in main_keys]
    print(f"  Novel expansion records to append: {len(novel):,}")

    # Write merged file: main + novel expansion
    print(f"\nWriting merged dataset: {MERGED_OUT}")
    MERGED_OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(MAIN_DATASET) as src, open(MERGED_OUT, "w") as dst:
        total = 0
        for line in src:
            dst.write(line)
            total += 1
        for r in novel:
            dst.write(json.dumps(r) + "\n")
            total += 1
    print(f"  Total records written: {total:,}")

    # Summary
    genes = Counter()
    uniprots = Counter()
    pmids = set()
    for r in novel:
        if r.get("anchor_gene_name"):
            genes[r["anchor_gene_name"]] += 1
        if r.get("uniprot_primary"):
            uniprots[r["uniprot_primary"]] += 1
        pmids.add(r["pmid"])

    main_genes = set()
    main_pmids = set()
    with open(MAIN_DATASET) as f:
        for line in f:
            r = json.loads(line)
            g = r.get("corrected_gene") or r.get("recovered_gene") or r.get("anchor_gene_name", "")
            if g:
                main_genes.add(g)
            main_pmids.add(r["pmid"])

    new_genes = set(genes) - main_genes
    new_pmids = set(pmids) - main_pmids

    stats_summary = {
        "per_batch": per_batch_stats,
        "total_expansion_records_pre_dedup": len(all_expansion),
        "duplicates_removed": dupes_removed,
        "collisions_with_main": collisions,
        "novel_records_appended": len(novel),
        "main_dataset_records": len(main_keys),
        "merged_dataset_records": total,
        "expansion_unique_genes": len(genes),
        "expansion_new_genes_vs_main": len(new_genes),
        "expansion_unique_pmids": len(pmids),
        "expansion_new_pmids_vs_main": len(new_pmids),
        "expansion_unique_uniprots": len(uniprots),
        "top_genes": genes.most_common(20),
    }
    with open(STATS_OUT, "w") as f:
        json.dump(stats_summary, f, indent=2)

    print("\n=== MERGE SUMMARY ===")
    print(f"Main dataset records:          {len(main_keys):,}")
    print(f"Novel expansion records:       {len(novel):,}")
    print(f"Merged total:                  {total:,}")
    print(f"New genes added to dataset:    {len(new_genes):,}")
    print(f"New PMIDs added to dataset:    {len(new_pmids):,}")
    print(f"Unique UniProt accs (new):     {len(uniprots):,}")
    print(f"\nOutput: {MERGED_OUT}")
    print(f"Stats:  {STATS_OUT}")


if __name__ == "__main__":
    main()
