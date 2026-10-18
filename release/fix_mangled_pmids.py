#!/usr/bin/env python3
"""Strip `_<suffix>` from mangled pmid values in the v2 release (step 4 of 4).

Some expansion-batch records (the `reprocess_filtered` batch, whose article
keys are `<PMID>_<GENE>`) carry pmid values like `37250416_ARPP21`. The
prefix is the real PMID and the gene is already in the `gene` field, so the
suffix is stripped. `PMC<digits>_<suffix>` is handled the same way.

Input:  release/output/crossvalqa_release_v2.jsonl
        (optionally a splits directory with train/val/test.jsonl)
Output: the same files rewritten in place. The original release is kept as
        <release>_pre_pmid_fix.jsonl (and <splits>/pre_pmid_fix/); if that
        backup already exists it is used as the source, so reruns are safe.

Example:
    python release/fix_mangled_pmids.py
    python release/fix_mangled_pmids.py --release R.jsonl --splits-dir S/
"""

import argparse
import json
import re
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
RELEASE = REPO_ROOT / "release/output/crossvalqa_release_v2.jsonl"
SPLIT_NAMES = ["train", "val", "test"]

PMID_SUFFIX_RE = re.compile(r"^([0-9]+)_.+$")
PMCID_SUFFIX_RE = re.compile(r"^(PMC[0-9]+)_.+$")


def clean_pmid(pmid):
    if not pmid:
        return pmid, False
    pmid = str(pmid)
    m = PMID_SUFFIX_RE.match(pmid) or PMCID_SUFFIX_RE.match(pmid)
    if m:
        return m.group(1), True
    return pmid, False


def fix_file(src, dst, label):
    n_total = 0
    n_fixed = 0
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            r = json.loads(line)
            cleaned, was_mangled = clean_pmid(r.get("pmid"))
            if was_mangled:
                r["pmid"] = cleaned
                n_fixed += 1
            fout.write(json.dumps(r) + "\n")
            n_total += 1
    print(f"  {label}: {n_total:,} records, {n_fixed:,} fixed")
    return n_total, n_fixed


def main():
    ap = argparse.ArgumentParser(description="Strip _<GENE> suffixes from mangled pmids (in place, with backup).")
    ap.add_argument("--release", type=Path, default=RELEASE,
                    help=f"release JSONL to fix (default: {RELEASE.relative_to(REPO_ROOT)})")
    ap.add_argument("--splits-dir", type=Path, default=None,
                    help="optional directory with train/val/test.jsonl to fix as well")
    args = ap.parse_args()

    release = args.release
    release_backup = release.with_name(release.stem + "_pre_pmid_fix.jsonl")
    targets = [("release", release)]

    # ── release ──
    if not release_backup.exists():
        print(f"Backing up {release} → {release_backup}")
        shutil.copy2(release, release_backup)
    else:
        print(f"Backup {release_backup} already exists; using it as source")
    tmp = release.with_suffix(".tmp.jsonl")
    fix_file(release_backup, tmp, "release")
    tmp.replace(release)

    # ── splits ──
    if args.splits_dir is not None:
        splits_backup_dir = args.splits_dir / "pre_pmid_fix"
        splits_backup_dir.mkdir(exist_ok=True)
        for split in SPLIT_NAMES:
            src_file = args.splits_dir / f"{split}.jsonl"
            if not src_file.exists():
                continue
            backup_file = splits_backup_dir / f"{split}.jsonl"
            if not backup_file.exists():
                shutil.copy2(src_file, backup_file)
            tmp = src_file.with_suffix(".tmp.jsonl")
            fix_file(backup_file, tmp, f"splits/{split}")
            tmp.replace(src_file)
            targets.append((f"splits/{split}", src_file))

    # ── verify ──
    print("\nVerifying no mangled pmids remain...")
    for label, path in targets:
        bad = 0
        with open(path) as f:
            for line in f:
                p = str(json.loads(line).get("pmid", ""))
                if PMID_SUFFIX_RE.match(p) or PMCID_SUFFIX_RE.match(p):
                    bad += 1
        print(f"  {label}: {bad} mangled pmids remaining")


if __name__ == "__main__":
    main()
