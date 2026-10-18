#!/usr/bin/env python3
"""
Anchoring step 2: fetch gene annotations from PubTator Central for all PMIDs.
Uses the PubTator3 BioC-JSON export endpoint with batched requests (network).

Input:
- anchoring/output/pmid_list.txt (make_pmid_list.py)

Output:
- anchoring/output/pubtator_annotations.jsonl, one record per PMID:
  {"pmid", "genes": [{"ncbi_gene_id", "text", "name", "offsets": [[offset, length], ...],
   "passage_type"}], "species": [taxon ids], "is_human"}
  PMIDs missing from the API response get an empty record with "no_annotations": true.
- anchoring/output/failed_pmids.txt (PMIDs whose batch failed after retries)

Resumable: PMIDs already present in the output file are skipped.
Used by tiered_anchor_v3.py for the single-gene, Tier 2 and Tier 3 decisions.
"""

import argparse
import json
import sys
import time
import requests
from pathlib import Path
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import config  # noqa: E402

# Configuration
ANCHORING_OUTPUT_DIR = config.PHASE4_DIR / "anchoring" / "output"
PMID_FILE = ANCHORING_OUTPUT_DIR / "pmid_list.txt"
OUTPUT_FILE = ANCHORING_OUTPUT_DIR / "pubtator_annotations.jsonl"
FAILED_FILE = ANCHORING_OUTPUT_DIR / "failed_pmids.txt"

BATCH_SIZE = 100  # PMIDs per request
DELAY = 0.34      # seconds between requests
MAX_RETRIES = 3

API_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson"

def load_pmids(pmid_file):
    """Load PMIDs from file."""
    with open(pmid_file) as f:
        return [line.strip() for line in f if line.strip()]

def load_completed(output_file):
    """Load already-fetched PMIDs to enable resume."""
    completed = set()
    if output_file.exists():
        with open(output_file) as f:
            for line in f:
                try:
                    data = json.loads(line)
                    if data.get("pmid"):
                        completed.add(data["pmid"])
                except:
                    pass
    return completed

def fetch_batch(pmids: list[str], session: requests.Session) -> dict:
    """Fetch annotations for a batch of PMIDs."""
    pmid_str = ",".join(pmids)

    for attempt in range(MAX_RETRIES):
        try:
            resp = session.get(API_URL, params={"pmids": pmid_str}, timeout=60)
            if resp.status_code == 200:
                return resp.json()
            elif resp.status_code == 429:  # Rate limited
                wait = (attempt + 1) * 5
                print(f"\nRate limited, waiting {wait}s...")
                time.sleep(wait)
            else:
                print(f"\nHTTP {resp.status_code} for batch starting {pmids[0]}")
                if attempt < MAX_RETRIES - 1:
                    time.sleep(2)
        except requests.exceptions.Timeout:
            print(f"\nTimeout for batch starting {pmids[0]}, retry {attempt+1}")
            time.sleep(2)
        except Exception as e:
            print(f"\nError: {e}")
            time.sleep(2)

    return None

def extract_genes(pubtator_doc: dict) -> dict:
    """Extract genes (all species) and paper-level species from a PubTator document."""
    # Handle PubTator3 format: id is at top level
    pmid = pubtator_doc.get("id", "").split("|")[0]  # Handle "12345|None" format

    genes = []
    species = set()

    # Process passages (title, abstract, full text sections)
    for passage in pubtator_doc.get("passages", []):
        passage_type = passage.get("infons", {}).get("type", "")
        passage_offset = passage.get("offset", 0)

        for annotation in passage.get("annotations", []):
            infons = annotation.get("infons", {})
            ann_type = infons.get("type", "")

            # Collect species - can be string or int
            if ann_type == "Species":
                identifier = infons.get("identifier")
                if identifier:
                    species.add(str(identifier))
                # Also check normalized array
                normalized = infons.get("normalized", [])
                for n in normalized:
                    species.add(str(n))

            # Collect genes
            elif ann_type == "Gene":
                gene_id = infons.get("identifier", "")  # NCBI Gene ID
                gene_text = annotation.get("text", "")
                gene_name = infons.get("name", "")  # Normalized gene name

                # Get character offsets
                locations = annotation.get("locations", [])
                offsets = [(loc.get("offset", 0), loc.get("length", 0)) for loc in locations]

                genes.append({
                    "ncbi_gene_id": str(gene_id) if gene_id else "",
                    "text": gene_text,
                    "name": gene_name,
                    "offsets": offsets,
                    "passage_type": passage_type
                })

    return {
        "pmid": pmid,
        "genes": genes,
        "species": list(species),
        "is_human": "9606" in species
    }

def main():
    parser = argparse.ArgumentParser(description="Fetch PubTator3 BioC-JSON gene annotations for all PMIDs")
    parser.add_argument("--pmid-file", type=Path, default=PMID_FILE, help="pmid_list.txt (make_pmid_list.py)")
    parser.add_argument("--output", type=Path, default=OUTPUT_FILE)
    parser.add_argument("--failed-output", type=Path, default=FAILED_FILE)
    args = parser.parse_args()

    print("Loading PMIDs...")
    all_pmids = load_pmids(args.pmid_file)
    completed = load_completed(args.output)

    remaining = [p for p in all_pmids if p not in completed]
    print(f"Total PMIDs: {len(all_pmids):,}")
    print(f"Already fetched: {len(completed):,}")
    print(f"Remaining: {len(remaining):,}")

    if not remaining:
        print("All PMIDs already fetched!")
        return

    # Create batches
    batches = [remaining[i:i+BATCH_SIZE] for i in range(0, len(remaining), BATCH_SIZE)]
    print(f"Batches to process: {len(batches)}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    failed_pmids = []

    session = requests.Session()
    session.headers.update({"Accept": "application/json"})

    with open(args.output, "a") as out_f, open(args.failed_output, "a") as fail_f:
        for batch in tqdm(batches, desc="Fetching"):
            result = fetch_batch(batch, session)

            if result:
                # PubTator3 returns {"PubTator3": [...]} format
                docs = result.get("PubTator3", [])
                if not docs and isinstance(result, list):
                    docs = result

                fetched_pmids = set()

                for doc in docs:
                    extracted = extract_genes(doc)
                    if extracted["pmid"]:  # Only write if we got a valid PMID
                        out_f.write(json.dumps(extracted) + "\n")
                        fetched_pmids.add(extracted["pmid"])

                # Track PMIDs that weren't in response (no annotations or not found)
                for pmid in batch:
                    if pmid not in fetched_pmids:
                        # Write empty record for PMIDs with no annotations
                        out_f.write(json.dumps({
                            "pmid": pmid,
                            "genes": [],
                            "species": [],
                            "is_human": False,
                            "no_annotations": True
                        }) + "\n")
            else:
                # Batch failed after retries
                for pmid in batch:
                    failed_pmids.append(pmid)
                    fail_f.write(pmid + "\n")

            out_f.flush()
            time.sleep(DELAY)

    print(f"\nDone! Failed PMIDs: {len(failed_pmids)}")
    print(f"Output: {args.output}")

if __name__ == "__main__":
    main()
