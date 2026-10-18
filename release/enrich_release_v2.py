#!/usr/bin/env python3
"""Enrich the v2 release file with structure metadata (step 2 of 4).

Strategy:
  - For the 416,188 v2 records that exist in the v1 release (matched by
    (pmid, mutation, question)), carry forward v1's enrichment wholesale:
    pdb_ids/pdb_urls/alphafold_id/alphafold_url, mutation_position,
    mutation_locations, wt_aa/mut_aa, mutated_sequence,
    has_structure/structure_evaluable.
  - For the 7,227 expansion records, compute enrichment fresh:
      * PDB/AlphaFold IDs + features for the ~1,000 new uniprots via the
        UniProt REST /uniprotkb/accessions?accessions=... JSON endpoint.
      * wt_aa/mut_aa/mutation_position from parse_protein_aa() on
        p.-prefix mutations; left null for c./bare mutations.
      * mutated_sequence by substituting the WT residue at mutation_position
        in ref_sequence (only when all three are available).
      * mutation_locations by cross-referencing mutation_position against
        the feature annotations for that uniprot.

Input:  release/output/crossvalqa_release_v2.jsonl   (from prepare_release_v2.py)
        release/input/crossvalqa_release.jsonl       (v1 enriched release)
        release/input/protein_3d_ids.json            (v1 uniprot -> PDB/AlphaFold ids)
        release/input/mutation_locations.json        (v1 uniprot -> feature list)
        UniProt REST unless --skip-fetch
Output: overwrites release/output/crossvalqa_release_v2.jsonl in place
        (with --limit, a pilot file <release>.jsonl.tmp is kept instead)
        release/output/protein_3d_ids_v2_new.json, features_v2_new.json (fetch caches)

Example:
    python release/enrich_release_v2.py --limit 1000 --skip-fetch   # pilot
    python release/enrich_release_v2.py
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = REPO_ROOT / "release/input"
OUTPUT_DIR = REPO_ROOT / "release/output"
V1_RELEASE = INPUT_DIR / "crossvalqa_release.jsonl"
V2_RELEASE = OUTPUT_DIR / "crossvalqa_release_v2.jsonl"
PROTEIN_3D_IDS = INPUT_DIR / "protein_3d_ids.json"
MUTATION_LOCATIONS = INPUT_DIR / "mutation_locations.json"
PROTEIN_3D_V2_EXTENSION = OUTPUT_DIR / "protein_3d_ids_v2_new.json"
FEATURES_V2_EXTENSION = OUTPUT_DIR / "features_v2_new.json"

AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Stop": "*",
}
RE_1LETTER = re.compile(r"^([A-Z])(\d+)([A-Z*])$")
RE_3LETTER = re.compile(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)$")

# Feature types the v1 cache kept (subset of all UniProt feature types).
V1_FEATURE_TYPES = {
    "Binding site", "Region", "Compositional bias", "Disulfide bond",
    "Domain", "Glycosylation", "Transmembrane", "Repeat", "Motif",
    "Site", "Active site", "Zinc finger", "Signal", "Coiled coil",
    "Lipidation", "DNA binding",
}


def parse_protein_aa(mutation: str, position_resolved: bool = False) -> tuple[str | None, str | None, int | None]:
    """Return (wt_aa, mut_aa, position) if mutation is a parseable protein substitution.

    Bare 1-letter form (V600E) is accepted only when position_resolved=True,
    since strings like T89C are also valid cDNA notation.
    """
    if not mutation:
        return None, None, None
    s = mutation.strip()
    if ":" in s:
        s = s.split(":")[-1]
    explicit_protein = s.startswith("p.")
    if explicit_protein:
        s = s[2:]
    m = RE_3LETTER.match(s)
    if m:
        wt = AA_3TO1.get(m.group(1))
        mut = AA_3TO1.get(m.group(3), "*")
        pos = int(m.group(2)) if wt else None
        if wt:
            return wt, mut, pos
    m = RE_1LETTER.match(s)
    if m and (explicit_protein or position_resolved):
        return m.group(1), m.group(3), int(m.group(2))
    return None, None, None


def fetch_uniprot_entries_json(accessions: list[str], batch: int = 100) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for i in range(0, len(accessions), batch):
        chunk = accessions[i:i + batch]
        url = (
            "https://rest.uniprot.org/uniprotkb/accessions"
            f"?accessions={','.join(chunk)}&format=json&size=500"
        )
        req = urllib.request.Request(url, headers={"User-Agent": "CrossValQA-release-prep/1.0"})
        try:
            with urllib.request.urlopen(req, timeout=90) as resp:
                payload = json.loads(resp.read())
        except Exception as exc:
            print(f"  batch {i}: ERROR {exc}", file=sys.stderr)
            continue
        for entry in payload.get("results", []):
            acc = entry.get("primaryAccession")
            if acc:
                out[acc] = entry
        print(f"  {i + len(chunk):>5}/{len(accessions):>5}  ({len(out):>5} entries cum)")
        if i + batch < len(accessions):
            time.sleep(0.3)
    return out


def extract_pdb_and_af(entry: dict) -> dict:
    pdb_ids: list[str] = []
    af_id: str | None = None
    for xref in entry.get("uniProtKBCrossReferences", []):
        db = xref.get("database")
        acc = xref.get("id")
        if db == "PDB" and acc:
            pdb_ids.append(acc)
        elif db == "AlphaFoldDB" and acc and af_id is None:
            af_id = f"AF-{acc}-F1"
    return {"pdb_ids": pdb_ids, "alphafold_id": af_id}


def extract_features(entry: dict) -> list[dict]:
    feats: list[dict] = []
    for f in entry.get("features", []):
        t = f.get("type")
        if t not in V1_FEATURE_TYPES:
            continue
        loc = f.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue
        feats.append({
            "type": t,
            "start": start,
            "end": end,
            "description": f.get("description", "") or "",
        })
    return feats


def locations_for_position(features: list[dict], position: int) -> list[dict] | None:
    if position is None or not features:
        return None
    hits = [f for f in features if f["start"] <= position <= f["end"]]
    return hits if hits else None


def pdb_url(pdb_id: str) -> str:
    return f"https://files.rcsb.org/download/{pdb_id}.pdb"


def af_url(af_id: str) -> str:
    return f"https://alphafold.ebi.ac.uk/files/{af_id}-model_v4.pdb"


def apply_mutation(ref_seq: str, position: int, wt: str, mut: str) -> str | None:
    if ref_seq is None or position is None or position < 1 or position > len(ref_seq):
        return None
    if ref_seq[position - 1] != wt:
        return None
    if mut == "*":
        return ref_seq[:position - 1]  # truncate at stop
    return ref_seq[:position - 1] + mut + ref_seq[position:]


def load_v1_enrichment_lookup() -> dict[tuple, dict]:
    """Build (pmid,mutation,question) -> enrichment-dict from v1 release."""
    enrichment_keys = [
        "pdb_ids", "pdb_urls", "alphafold_id", "alphafold_url",
        "mutation_position", "mutation_locations",
        "wt_aa", "mut_aa", "mutated_sequence",
        "has_structure", "structure_evaluable",
    ]
    lookup: dict[tuple, dict] = {}
    with open(V1_RELEASE) as f:
        for line in f:
            r = json.loads(line)
            key = (r.get("pmid"), r.get("mutation"), r.get("question"))
            lookup[key] = {k: r.get(k) for k in enrichment_keys}
    return lookup


def load_existing_3d_cache() -> dict[str, dict]:
    """Load protein_3d_ids.json (v1) and any prior v2 extension."""
    cache: dict[str, dict] = {}
    if PROTEIN_3D_IDS.exists():
        d = json.loads(PROTEIN_3D_IDS.read_text())
        cache.update(d.get("mappings", {}))
    if PROTEIN_3D_V2_EXTENSION.exists():
        ext = json.loads(PROTEIN_3D_V2_EXTENSION.read_text())
        cache.update(ext)
    return cache


def load_existing_features_cache() -> dict[str, list[dict]]:
    cache: dict[str, list[dict]] = {}
    if MUTATION_LOCATIONS.exists():
        d = json.loads(MUTATION_LOCATIONS.read_text())
        cache.update(d.get("features_by_uniprot", {}))
    if FEATURES_V2_EXTENSION.exists():
        ext = json.loads(FEATURES_V2_EXTENSION.read_text())
        cache.update(ext)
    return cache


def main(argv: list[str] | None = None) -> int:
    global V1_RELEASE, V2_RELEASE, PROTEIN_3D_IDS, MUTATION_LOCATIONS
    global PROTEIN_3D_V2_EXTENSION, FEATURES_V2_EXTENSION
    ap = argparse.ArgumentParser(description="Fill structure/position fields of the v2 release in place.")
    ap.add_argument("--limit", type=int, default=None, help="Pilot: process only the first N records (output kept as .tmp)")
    ap.add_argument("--skip-fetch", action="store_true", help="Skip UniProt REST fetch; use only local caches")
    ap.add_argument("--release", type=Path, default=V2_RELEASE, help=f"v2 release file, rewritten in place (default: {V2_RELEASE.relative_to(REPO_ROOT)})")
    ap.add_argument("--v1-release", type=Path, default=V1_RELEASE, help=f"v1 enriched release (default: {V1_RELEASE.relative_to(REPO_ROOT)})")
    ap.add_argument("--protein-3d-ids", type=Path, default=PROTEIN_3D_IDS, help=f"v1 structure-id cache (default: {PROTEIN_3D_IDS.relative_to(REPO_ROOT)})")
    ap.add_argument("--mutation-locations", type=Path, default=MUTATION_LOCATIONS, help=f"v1 feature cache (default: {MUTATION_LOCATIONS.relative_to(REPO_ROOT)})")
    ap.add_argument("--cache-dir", type=Path, default=OUTPUT_DIR, help=f"where the v2 fetch caches are read/written (default: {OUTPUT_DIR.relative_to(REPO_ROOT)})")
    args = ap.parse_args(argv)
    V1_RELEASE, V2_RELEASE = args.v1_release, args.release
    PROTEIN_3D_IDS, MUTATION_LOCATIONS = args.protein_3d_ids, args.mutation_locations
    PROTEIN_3D_V2_EXTENSION = args.cache_dir / PROTEIN_3D_V2_EXTENSION.name
    FEATURES_V2_EXTENSION = args.cache_dir / FEATURES_V2_EXTENSION.name

    print("[1/5] Loading v1 enrichment lookup + existing caches")
    v1_lookup = load_v1_enrichment_lookup()
    structure_cache = load_existing_3d_cache()
    features_cache = load_existing_features_cache()
    print(f"  v1 lookup: {len(v1_lookup):,} keys")
    print(f"  structure cache: {len(structure_cache):,} uniprots")
    print(f"  features cache: {len(features_cache):,} uniprots")

    print("[2/5] Identifying uniprots needing fetch")
    all_v2_uniprots: set[str] = set()
    with open(V2_RELEASE) as f:
        for line in f:
            r = json.loads(line)
            up = r.get("uniprot_primary")
            if up:
                all_v2_uniprots.add(up)
    missing_struct = sorted(all_v2_uniprots - set(structure_cache.keys()))
    missing_feats = sorted(all_v2_uniprots - set(features_cache.keys()))
    needed = sorted(set(missing_struct) | set(missing_feats))
    print(f"  uniprots in v2: {len(all_v2_uniprots):,}")
    print(f"  missing structure xrefs: {len(missing_struct):,}")
    print(f"  missing features: {len(missing_feats):,}")
    print(f"  to fetch (union): {len(needed):,}")

    if needed and not args.skip_fetch:
        print("[3/5] Fetching UniProt JSON entries")
        entries = fetch_uniprot_entries_json(needed)
        new_struct: dict[str, dict] = {}
        new_feats: dict[str, list[dict]] = {}
        for acc, entry in entries.items():
            sx = extract_pdb_and_af(entry)
            new_struct[acc] = sx
            structure_cache[acc] = sx
            fs = extract_features(entry)
            new_feats[acc] = fs
            features_cache[acc] = fs
        # Persist extensions
        PROTEIN_3D_V2_EXTENSION.parent.mkdir(parents=True, exist_ok=True)
        prior_struct = {}
        if PROTEIN_3D_V2_EXTENSION.exists():
            prior_struct = json.loads(PROTEIN_3D_V2_EXTENSION.read_text())
        prior_struct.update(new_struct)
        PROTEIN_3D_V2_EXTENSION.write_text(json.dumps(prior_struct))
        prior_feats = {}
        if FEATURES_V2_EXTENSION.exists():
            prior_feats = json.loads(FEATURES_V2_EXTENSION.read_text())
        prior_feats.update(new_feats)
        FEATURES_V2_EXTENSION.write_text(json.dumps(prior_feats))
        print(f"  fetched {len(entries):,} entries")
        print(f"  persisted extensions: {PROTEIN_3D_V2_EXTENSION.name}, {FEATURES_V2_EXTENSION.name}")
    elif args.skip_fetch:
        print("[3/5] --skip-fetch: using caches only")

    print("[4/5] Streaming v2 release, enriching in place")
    tmp_out = V2_RELEASE.with_suffix(".jsonl.tmp")
    n_carry = 0
    n_expansion = 0
    n_protein_hgvs_parsed = 0
    n_mutated_seq = 0
    n_structure = 0
    n_evaluable = 0
    total = 0
    with open(V2_RELEASE) as fin, open(tmp_out, "w") as fout:
        for line in fin:
            if args.limit is not None and total >= args.limit:
                break
            r = json.loads(line)
            total += 1
            key = (r.get("pmid"), r.get("mutation"), r.get("question"))
            up = r.get("uniprot_primary")
            ref_seq = r.get("ref_sequence")
            mutation = r.get("mutation")

            if key in v1_lookup:
                # Carry forward v1 enrichment (authoritative for overlap records)
                v1e = v1_lookup[key]
                for k, v in v1e.items():
                    r[k] = v
                n_carry += 1
            else:
                # Expansion record: enrich fresh
                n_expansion += 1
                sx = structure_cache.get(up, {}) if up else {}
                pdb_ids = sx.get("pdb_ids", []) or []
                af_id = sx.get("alphafold_id")
                r["pdb_ids"] = pdb_ids
                r["pdb_urls"] = [pdb_url(p) for p in pdb_ids]
                r["alphafold_id"] = af_id
                r["alphafold_url"] = af_url(af_id) if af_id else None

                wt, mut, pos = parse_protein_aa(mutation or "", position_resolved=False)
                if wt is not None and pos is not None:
                    n_protein_hgvs_parsed += 1
                r["wt_aa"] = wt
                r["mut_aa"] = mut
                r["mutation_position"] = pos

                feats = features_cache.get(up, []) if up else []
                r["mutation_locations"] = locations_for_position(feats, pos)

                if ref_seq and wt and mut and pos:
                    mseq = apply_mutation(ref_seq, pos, wt, mut)
                    r["mutated_sequence"] = mseq
                    if mseq is not None:
                        n_mutated_seq += 1
                else:
                    r["mutated_sequence"] = None

                r["has_structure"] = bool(pdb_ids or af_id)
                r["structure_evaluable"] = r["has_structure"] and r["mutation_position"] is not None

            if r.get("has_structure"):
                n_structure += 1
            if r.get("structure_evaluable"):
                n_evaluable += 1
            fout.write(json.dumps(r) + "\n")

    print(f"[5/5] Summary")
    print(f"  total written: {total:,}")
    print(f"  carry-forward from v1: {n_carry:,}")
    print(f"  expansion enriched fresh: {n_expansion:,}")
    print(f"  expansion records with protein HGVS parsed: {n_protein_hgvs_parsed:,}")
    print(f"  expansion records with mutated_sequence computed: {n_mutated_seq:,}")
    print(f"  has_structure (all records): {n_structure:,} ({100*n_structure/total:.1f}%)")
    print(f"  structure_evaluable (all records): {n_evaluable:,} ({100*n_evaluable/total:.1f}%)")

    if args.limit is None:
        tmp_out.replace(V2_RELEASE)
        print(f"  wrote {V2_RELEASE}")
    else:
        print(f"  pilot output kept at {tmp_out} (not swapped in)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
