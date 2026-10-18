"""Fetch UniProt protein sequences for every unique uniprot_primary in v2.

Input:  the v2 dataset JSONL (field `uniprot_primary`).
Output: audit/output/uniprot_sequences.json as {acc: sequence}. Re-running
        only fetches accessions missing from the cache.
Network: UniProt REST (rest.uniprot.org/uniprotkb/accessions, FASTA).

Example (from the repo root):
    python audit/fetch_uniprot_sequences.py \
        --v2 uniprot/output/m2tqa_gene_corrected_v2.jsonl
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import urllib.request
import urllib.parse
import urllib.error

REPO_ROOT = Path(__file__).resolve().parents[1]
V2 = REPO_ROOT / "uniprot/output/m2tqa_gene_corrected_v2.jsonl"
CACHE = REPO_ROOT / "audit/output/uniprot_sequences.json"
BATCH_SIZE = 100  # accessions per request
RETRY_DELAY = 5
MAX_RETRIES = 3


def collect_accessions() -> set[str]:
    accs: set[str] = set()
    with V2.open() as f:
        for line in f:
            d = json.loads(line)
            acc = d.get("uniprot_primary")
            if acc:
                accs.add(acc)
    return accs


def parse_fasta(text: str) -> dict[str, str]:
    """UniProt FASTA: > sp|ACC|NAME ... \n SEQ\n"""
    out: dict[str, str] = {}
    cur_acc = None
    cur_seq: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur_acc is not None:
                out[cur_acc] = "".join(cur_seq)
            parts = line[1:].split("|")
            cur_acc = parts[1] if len(parts) >= 2 else None
            cur_seq = []
        elif line:
            cur_seq.append(line.strip())
    if cur_acc is not None:
        out[cur_acc] = "".join(cur_seq)
    return out


def fetch_batch(accs: list[str]) -> dict[str, str]:
    url = (
        "https://rest.uniprot.org/uniprotkb/accessions?"
        + urllib.parse.urlencode({
            "accessions": ",".join(accs),
            "format": "fasta",
            "size": str(len(accs)),
        })
    )
    last_err = None
    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return parse_fasta(resp.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError) as e:
            last_err = e
            time.sleep(RETRY_DELAY * (attempt + 1))
    raise RuntimeError(f"failed after {MAX_RETRIES} attempts: {last_err}")


def main() -> None:
    global V2, CACHE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--v2", type=Path, default=V2,
                    help=f"dataset JSONL with uniprot_primary (default: {V2.relative_to(REPO_ROOT)})")
    ap.add_argument("--cache", type=Path, default=CACHE,
                    help=f"sequence cache JSON, read and updated in place (default: {CACHE.relative_to(REPO_ROOT)})")
    args = ap.parse_args()
    V2, CACHE = args.v2, args.cache
    CACHE.parent.mkdir(parents=True, exist_ok=True)

    cache: dict[str, str] = {}
    if CACHE.exists():
        cache = json.loads(CACHE.read_text())
        print(f"Loaded {len(cache):,} cached sequences from {CACHE.name}")

    print("Scanning v2 for unique uniprot_primary accessions...")
    needed = collect_accessions()
    print(f"Unique accessions in v2: {len(needed):,}")

    todo = sorted(needed - set(cache))
    print(f"Already cached: {len(needed) - len(todo):,}")
    print(f"To fetch: {len(todo):,}")
    if not todo:
        print("Nothing to fetch.")
        return

    fetched = 0
    failed = 0
    started = time.time()
    for i in range(0, len(todo), BATCH_SIZE):
        batch = todo[i : i + BATCH_SIZE]
        try:
            seqs = fetch_batch(batch)
        except Exception as e:
            print(f"  batch {i}: {e}", flush=True)
            failed += len(batch)
            continue

        cache.update(seqs)
        fetched += len(seqs)
        # Periodic checkpoint every ~10 batches.
        if (i // BATCH_SIZE) % 10 == 9:
            CACHE.write_text(json.dumps(cache))
            elapsed = time.time() - started
            rate = fetched / max(elapsed, 1)
            print(
                f"  fetched {fetched:,}/{len(todo):,} ({rate:.0f}/s, "
                f"checkpoint -> {CACHE.name})",
                flush=True,
            )

    CACHE.write_text(json.dumps(cache))
    print()
    print(f"Total fetched this run: {fetched:,}")
    print(f"Failed (no FASTA returned in batch): {failed:,}")
    print(f"Cache size now: {len(cache):,} sequences")
    not_returned = [a for a in todo if a not in cache]
    print(f"Accessions requested but not returned by UniProt: {len(not_returned):,}")
    if not_returned[:10]:
        print(f"  examples: {not_returned[:10]}")


if __name__ == "__main__":
    main()
