#!/usr/bin/env python3
"""Stage the MutQA release for Hugging Face as per-config Parquet + gzipped JSONL.

Configs:
  homology  homology_split/{train,test}.jsonl  -> <out>/homology/{train,test}.{parquet,jsonl.gz} + manifest.json
  random    random_split/{train,test}.jsonl    -> <out>/random/{train,test}.{parquet,jsonl.gz} + manifest.json
  full      crossvalqa_release_v3.jsonl (all rows, unfiltered) -> <out>/full/all.{parquet,jsonl.gz}

Parquet rows are normalized for a stable Arrow schema: phase1/phase2_evidence_ids
become list[str], and dict entries of mutation_locations become strings
("<type>: <description> (<start>-<end>)"). The .jsonl.gz files are byte copies
of the source JSONL. A stage_sizes.json summary is written to <out>.
This script only writes local files; uploading is done separately.

Input:  release/output/{homology_split,random_split}/, release/output/crossvalqa_release_v3.jsonl
Output: release/output/hf_staging/
Requires: pyarrow

Example:
    python release/stage_hf.py
    python release/stage_hf.py --configs full --full-source release/output/crossvalqa_release_v3.jsonl
"""

import argparse
import gzip
import json
import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "release/output"
STAGING = OUTPUT_DIR / "hf_staging"
SPLIT_DIRS = {
    "homology": OUTPUT_DIR / "homology_split",
    "random": OUTPUT_DIR / "random_split",
}
FULL_SOURCE = OUTPUT_DIR / "crossvalqa_release_v3.jsonl"
SPLITS = ("train", "test")
ALL_CONFIGS = ("homology", "random", "full")

EVIDENCE_FIELDS = ("phase1_evidence_ids", "phase2_evidence_ids")


def load_jsonl(path):
    with open(path) as f:
        for line in f:
            yield json.loads(line)


def gzip_file(src, dst):
    with open(src, "rb") as fin, gzip.open(dst, "wb", compresslevel=6) as fout:
        shutil.copyfileobj(fin, fout, length=1 << 20)


def _fmt_location(x):
    if isinstance(x, str):
        return x
    if isinstance(x, dict):
        t = x.get("type") or ""
        desc = x.get("description") or ""
        start = x.get("start")
        end = x.get("end")
        head = f"{t}: {desc}" if t and desc else (t or desc)
        if start is not None and end is not None:
            return f"{head} ({start}-{end})" if head else f"({start}-{end})"
        return head or ""
    return str(x)


def _normalize_row(row):
    for k in EVIDENCE_FIELDS:
        v = row.get(k)
        if isinstance(v, list):
            row[k] = [str(x) for x in v]
    ml = row.get("mutation_locations")
    if isinstance(ml, list):
        row["mutation_locations"] = [_fmt_location(x) for x in ml]
    return row


def jsonl_to_parquet(src, dst):
    rows = [_normalize_row(r) for r in load_jsonl(src)]
    table = pa.Table.from_pylist(rows)
    pq.write_table(table, dst, compression="zstd", compression_level=9)
    return len(rows)


def stage_file(src, cfg_out, stem, cfg, sizes):
    pq_out = cfg_out / f"{stem}.parquet"
    gz_out = cfg_out / f"{stem}.jsonl.gz"
    n = jsonl_to_parquet(src, pq_out)
    gzip_file(src, gz_out)
    sizes[stem] = {
        "records": n,
        "parquet_bytes": pq_out.stat().st_size,
        "jsonl_gz_bytes": gz_out.stat().st_size,
    }
    print(
        f"  {cfg:9s} {stem:5s}  {n:>7,} rec  "
        f"parquet={pq_out.stat().st_size/1e6:7.1f} MB  "
        f"gz={gz_out.stat().st_size/1e6:7.1f} MB"
    )


def main():
    ap = argparse.ArgumentParser(description="Stage MutQA configs as Parquet + jsonl.gz for Hugging Face (local only).")
    ap.add_argument("--configs", nargs="+", choices=ALL_CONFIGS, default=list(ALL_CONFIGS))
    ap.add_argument("--homology-dir", type=Path, default=SPLIT_DIRS["homology"], help="homology split dir (default: release/output/homology_split)")
    ap.add_argument("--random-dir", type=Path, default=SPLIT_DIRS["random"], help="random split dir (default: release/output/random_split)")
    ap.add_argument("--full-source", type=Path, default=FULL_SOURCE, help="release JSONL for the full config (default: release/output/crossvalqa_release_v3.jsonl)")
    ap.add_argument("--out", type=Path, default=STAGING, help="staging dir (default: release/output/hf_staging)")
    args = ap.parse_args()
    split_dirs = {"homology": args.homology_dir, "random": args.random_dir}

    args.out.mkdir(parents=True, exist_ok=True)
    sizes = {}
    for cfg in args.configs:
        cfg_out = args.out / cfg
        if cfg_out.exists():
            shutil.rmtree(cfg_out)
        cfg_out.mkdir(parents=True)
        sizes[cfg] = {}
        if cfg == "full":
            stage_file(args.full_source, cfg_out, "all", cfg, sizes[cfg])
            continue
        src_dir = split_dirs[cfg]
        for split in SPLITS:
            stage_file(src_dir / f"{split}.jsonl", cfg_out, split, cfg, sizes[cfg])
        shutil.copy2(src_dir / "manifest.json", cfg_out / "manifest.json")

    (args.out / "stage_sizes.json").write_text(json.dumps(sizes, indent=2))
    print(f"\nDone. Staging dir: {args.out}")


if __name__ == "__main__":
    main()
