#!/usr/bin/env python3
"""Merge an isoform-aware p. enrichment table into the v2 release (step 3 of 4).

For every record whose mutation is p.-notation, looks up (GENE, protein_hgvs)
in an externally produced TSV and:
  - sets matched_isoform / matched_isoform_status / matched_isoform_reason
    (null for unmatched and non-p. records),
  - fills mutation_position (from start_pos) and mutated_sequence only where
    they are still empty,
  - recomputes structure_evaluable = has_structure AND mutation_position set.

The TSV comes from a collaborator's isoform-matching pipeline that is not part
of this repository and is not distributed. Required tab-separated columns:
  gene, protein_hgvs (e.g. "BRAF:p.V600E" or "p.V600E"), chosen_isoid,
  status, reason, start_pos, mutated_sequence

Input:  release/output/crossvalqa_release_v2.jsonl + --isoform-tsv
Output: release/output/crossvalqa_release_v2.jsonl rewritten in place;
        the pre-merge file is kept as crossvalqa_release_v2_pre_isoform.jsonl

Example:
    python release/merge_isoform_enrichment.py --isoform-tsv release/input/protein_hgvs_with_sequences.tsv
"""

import argparse
import csv
import json
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
RELEASE = REPO_ROOT / "release/output/crossvalqa_release_v2.jsonl"


def load_isoform_map(path):
    d = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = (row["gene"].strip().upper(), row["protein_hgvs"].strip())
            d[key] = row
    return d


def parse_int(s):
    if not s:
        return None
    try:
        return int(float(s))
    except (ValueError, TypeError):
        return None


def main():
    ap = argparse.ArgumentParser(description="Merge isoform-aware p. enrichment into the v2 release (in place).")
    ap.add_argument("--isoform-tsv", type=Path, required=True, help="external isoform-matching TSV (see module docstring)")
    ap.add_argument("--release", type=Path, default=RELEASE, help=f"v2 release file, rewritten in place (default: {RELEASE.relative_to(REPO_ROOT)})")
    args = ap.parse_args()
    release = args.release
    BACKUP = release.with_name(release.stem + "_pre_isoform.jsonl")
    TMP_OUT = release.with_name(release.stem + ".tmp.jsonl")

    print(f"Loading isoform TSV from {args.isoform_tsv}...")
    isoforms = load_isoform_map(args.isoform_tsv)
    print(f"  {len(isoforms):,} rows")

    if not BACKUP.exists():
        print(f"Backing up {release} → {BACKUP}")
        shutil.copy2(release, BACKUP)
    else:
        print(f"Backup already exists at {BACKUP}; leaving it alone")

    n_total = 0
    n_p = 0
    n_matched = 0
    n_filled_pos = 0
    n_filled_mut = 0
    n_status_set = 0
    n_struct_eval_true = 0

    print(f"Streaming {BACKUP} → {TMP_OUT}")
    with open(BACKUP) as fin, open(TMP_OUT, "w") as fout:
        for line in fin:
            r = json.loads(line)
            n_total += 1

            mut = (r.get("mutation") or "").strip()
            core = mut.split(":")[-1] if ":" in mut else mut

            if core.startswith("p."):
                n_p += 1
                gene = (r.get("gene") or "").strip().upper()
                hgvs_with_gene = f"{gene}:{core}" if ":" not in mut else mut
                drow = isoforms.get((gene, hgvs_with_gene)) or isoforms.get((gene, core))

                if drow is not None:
                    n_matched += 1

                    chosen = (drow.get("chosen_isoid") or "").strip() or None
                    status = (drow.get("status") or "").strip() or None
                    reason = (drow.get("reason") or "").strip() or None
                    r["matched_isoform"] = chosen
                    r["matched_isoform_status"] = status
                    r["matched_isoform_reason"] = reason
                    if status:
                        n_status_set += 1

                    his_pos = parse_int(drow.get("start_pos"))
                    if r.get("mutation_position") is None and his_pos is not None:
                        r["mutation_position"] = his_pos
                        n_filled_pos += 1

                    his_mut = (drow.get("mutated_sequence") or "").strip()
                    his_mut_valid = his_mut and len(his_mut) > 1
                    if not r.get("mutated_sequence") and his_mut_valid:
                        r["mutated_sequence"] = his_mut
                        n_filled_mut += 1
                else:
                    r.setdefault("matched_isoform", None)
                    r.setdefault("matched_isoform_status", None)
                    r.setdefault("matched_isoform_reason", None)
            else:
                r.setdefault("matched_isoform", None)
                r.setdefault("matched_isoform_status", None)
                r.setdefault("matched_isoform_reason", None)

            r["structure_evaluable"] = bool(r.get("has_structure")) and r.get("mutation_position") is not None
            if r["structure_evaluable"]:
                n_struct_eval_true += 1

            fout.write(json.dumps(r) + "\n")

    print(f"\n  total records:              {n_total:,}")
    print(f"  p. records:                 {n_p:,}")
    print(f"  matched to isoform TSV:     {n_matched:,}")
    print(f"  mutation_position filled:   {n_filled_pos:,}")
    print(f"  mutated_sequence filled:    {n_filled_mut:,}")
    print(f"  matched_isoform_status set: {n_status_set:,}")
    print(f"  structure_evaluable=true:   {n_struct_eval_true:,}")

    print(f"\nSwapping {TMP_OUT} → {release}")
    TMP_OUT.replace(release)
    print("Done.")


if __name__ == "__main__":
    main()
