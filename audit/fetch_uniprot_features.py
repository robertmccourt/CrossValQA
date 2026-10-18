"""Fetch UniProt feature annotations for every accession in v2.

Input:  the v2 dataset JSONL (field `uniprot_primary`).
Output: audit/output/uniprot_features.json as
        {acc: {"length": int, "features": [...]}}. Incremental, like
        fetch_uniprot_sequences.py.
Network: UniProt REST (rest.uniprot.org/uniprotkb/accessions, JSON).

Example (from the repo root):
    python audit/fetch_uniprot_features.py \
        --v2 uniprot/output/m2tqa_gene_corrected_v2.jsonl

Filters to the ~15 feature types the verifier checks against; drops noise
types (Natural variant, Helix, Beta strand, Mutagenesis, Sequence conflict,
Compositional bias, etc.) which are not deterministically verifiable claims
the verifier can use.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import urllib.request
import urllib.parse
import urllib.error

REPO_ROOT = Path(__file__).resolve().parents[1]
V2 = REPO_ROOT / "uniprot/output/m2tqa_gene_corrected_v2.jsonl"
CACHE = REPO_ROOT / "audit/output/uniprot_features.json"
BATCH_SIZE = 100  # accessions per request
RETRY_DELAY = 5
MAX_RETRIES = 3

# Feature types we keep. Everything else (variants, secondary structure,
# mutagenesis, sequence conflicts, compositional bias) is dropped at fetch
# time — those are not claim families the verifier can check.
KEEP_TYPES = frozenset({
    "Domain",
    "Region",
    "Motif",
    "Repeat",
    "Zinc finger",
    "Active site",
    "Binding site",
    "Site",
    "Modified residue",
    "Glycosylation",
    "Disulfide bond",
    "Cross-link",
    "Lipidation",
    "Transmembrane",
    "Intramembrane",
    "Topological domain",
    "Signal peptide",
    "Propeptide",
    "Transit peptide",
    "Chain",
    "Peptide",
    "Coiled coil",
    "DNA binding",
})


def collect_accessions() -> set[str]:
    accs: set[str] = set()
    with V2.open() as f:
        for line in f:
            d = json.loads(line)
            acc = d.get("uniprot_primary")
            if acc:
                accs.add(acc)
    return accs


def slim_feature(f: dict) -> dict:
    """Strip a UniProt feature down to what the verifier actually uses."""
    loc = f.get("location", {})
    start = (loc.get("start") or {}).get("value")
    end = (loc.get("end") or {}).get("value")
    out = {
        "type": f["type"],
        "start": start,
        "end": end,
        "description": f.get("description", ""),
    }
    # Ligand info on Binding sites is the most useful detail beyond position.
    lig = f.get("ligand")
    if lig:
        out["ligand"] = lig.get("name", "")
    # Disulfide partner: encoded as a second location.
    if f["type"] == "Disulfide bond":
        # UniProt encodes the partner as the end position; for inter-chain
        # disulfides there's a `featureCrossReferences` entry. Keep the basic
        # range — interpretation done in the loader.
        pass
    # Evidence codes (compact form). We keep just the highest-confidence code
    # per feature so we can filter computational-only features later.
    codes = []
    for ev in f.get("evidences", []) or []:
        c = ev.get("evidenceCode")
        if c and c not in codes:
            codes.append(c)
    if codes:
        out["evidence"] = codes
    return out


def parse_response(data: dict) -> dict[str, dict]:
    """UniProt JSON → {acc: {length, features: [slim]}}."""
    out: dict[str, dict] = {}
    for entry in data.get("results", []):
        acc = entry["primaryAccession"]
        length = entry.get("sequence", {}).get("length", 0)
        feats = []
        for f in entry.get("features", []) or []:
            t = f.get("type")
            if t in KEEP_TYPES:
                feats.append(slim_feature(f))
        out[acc] = {"length": length, "features": feats}
    return out


def fetch_batch(accs: list[str]) -> dict[str, dict]:
    url = (
        "https://rest.uniprot.org/uniprotkb/accessions?"
        + urllib.parse.urlencode({
            "accessions": ",".join(accs),
            "format": "json",
            "size": str(len(accs)),
        })
    )
    last_err = None
    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                return parse_response(json.loads(resp.read()))
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
                    help=f"feature cache JSON, read and updated in place (default: {CACHE.relative_to(REPO_ROOT)})")
    args = ap.parse_args()
    V2, CACHE = args.v2, args.cache
    CACHE.parent.mkdir(parents=True, exist_ok=True)

    cache: dict[str, dict] = {}
    if CACHE.exists():
        cache = json.loads(CACHE.read_text())
        print(f"Loaded {len(cache):,} cached entries from {CACHE.name}")

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
            entries = fetch_batch(batch)
        except Exception as e:
            print(f"  batch {i}: {e}", flush=True)
            failed += len(batch)
            continue

        cache.update(entries)
        fetched += len(entries)
        if (i // BATCH_SIZE) % 5 == 4:
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
    print(f"Failed (no entry returned in batch): {failed:,}")
    print(f"Cache size now: {len(cache):,} entries")
    print(f"Cache file size: {CACHE.stat().st_size / 1024 / 1024:.1f} MB")
    not_returned = [a for a in todo if a not in cache]
    print(f"Accessions requested but not returned: {len(not_returned):,}")
    if not_returned[:10]:
        print(f"  examples: {not_returned[:10]}")

    # Summary stats
    total_features = sum(len(e["features"]) for e in cache.values())
    print(f"\nTotal features across cache: {total_features:,}")
    print(f"  avg per accession: {total_features / max(len(cache), 1):.1f}")
    type_counts: dict[str, int] = {}
    for e in cache.values():
        for f in e["features"]:
            type_counts[f["type"]] = type_counts.get(f["type"], 0) + 1
    print("\nFeature type distribution:")
    for t, c in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"  {t}: {c:,}")


if __name__ == "__main__":
    main()
