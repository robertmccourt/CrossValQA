"""Homology-aware and random train/test split generator for MutQA.

Prevents sequence-level leakage by clustering wild-type protein sequences
with MMseqs2 and assigning whole clusters to a single split. To close the
representative-vs-pairwise gap inherent to MMseqs2 clustering, we audit the
split with a GraphPart-style criterion (Teufel et al., 2023): every test
sequence with >= identity_threshold to any training sequence is discarded
from the test split. Train is left untouched.

Release config:
    - Source: crossvalqa_release_v3.jsonl (401,406 rows)
    - Filter: grounding == 'cross-grounded' AND uniprot_primary AND ref_sequence
      (330,387 rows, 5,437 uniprots)
    - MMseqs2: --min-seq-id 0.4 -c 0.8 --cov-mode 0
    - Ratios: 80/20 train/test (no validation split), seed 42,
      whole clusters assigned by seeded greedy capacity-fill on record counts
    - GraphPart-style audit: discard test uniprots with >= 40% identity to
      any train uniprot (test side only)
    - Random baseline (build_random_split): same filter, uniprot-level
      singleton clusters, same 80/20 and seed

Input:  release/output/crossvalqa_release_v3.jsonl; `mmseqs` on PATH (or --mmseqs)
Output: release/output/homology_split/{train,test}.jsonl + manifest.json
        release/output/random_split/{train,test}.jsonl + manifest.json

Run directly::

    python release/create_splits.py both            # or: homology | random
    python release/create_splits.py homology --input R.jsonl --mmseqs /path/to/mmseqs

Or import (with release/ on sys.path)::

    from create_splits import build_split, verify_split
    split = build_split("crossvalqa_release_v3.jsonl")
    split.write("splits/")
"""

from __future__ import annotations

import argparse
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JSONL = str(REPO_ROOT / "release/output/crossvalqa_release_v3.jsonl")
DEFAULT_OUTPUT_DIR = str(REPO_ROOT / "release/output/homology_split")
DEFAULT_RANDOM_OUTPUT_DIR = str(REPO_ROOT / "release/output/random_split")
DEFAULT_MMSEQS = "mmseqs"


@dataclass
class Split:
    uniprot_to_split: dict[str, str]
    cluster_to_split: dict[str, str]
    cluster_members: dict[str, list[str]]
    source_jsonl: str
    discarded_test_uniprots: set[str] = field(default_factory=set)
    manifest: dict = field(default_factory=dict)

    def write(self, output_dir: str) -> None:
        os.makedirs(output_dir, exist_ok=True)
        active_splits = sorted(set(self.uniprot_to_split.values()))
        paths = {name: os.path.join(output_dir, f"{name}.jsonl") for name in active_splits}
        handles = {name: open(path, "w") for name, path in paths.items()}

        per_split_counts = Counter()
        per_split_qtype = {name: Counter() for name in active_splits}
        per_split_uniprots = {name: set() for name in active_splits}
        dropped_filter = 0
        dropped_graphpart = 0

        try:
            with open(self.source_jsonl) as f:
                for line in f:
                    if not line.strip():
                        continue
                    rec = json.loads(line)
                    if rec.get("grounding") != "cross-grounded":
                        dropped_filter += 1
                        continue
                    up = rec.get("uniprot_primary")
                    if up is None or rec.get("ref_sequence") is None:
                        dropped_filter += 1
                        continue
                    if up in self.discarded_test_uniprots:
                        dropped_graphpart += 1
                        continue
                    split = self.uniprot_to_split.get(up)
                    if split is None:
                        dropped_filter += 1
                        continue
                    rec["split"] = split
                    handles[split].write(json.dumps(rec) + "\n")
                    per_split_counts[split] += 1
                    per_split_qtype[split][rec.get("question_type", "unknown")] += 1
                    per_split_uniprots[split].add(up)
        finally:
            for h in handles.values():
                h.close()

        for name, path in paths.items():
            print(f"Wrote {path} ({per_split_counts[name]:,} records, {len(per_split_uniprots[name]):,} uniprots)")
        print(f"GraphPart discard: {len(self.discarded_test_uniprots):,} test uniprots, {dropped_graphpart:,} records")

        self.manifest["split_record_counts"] = dict(per_split_counts)
        self.manifest["split_uniprot_counts"] = {k: len(v) for k, v in per_split_uniprots.items()}
        self.manifest["split_question_types"] = {
            name: dict(ctr) for name, ctr in per_split_qtype.items()
        }
        self.manifest["dropped_records"] = dropped_filter
        self.manifest["graphpart_discarded_records"] = dropped_graphpart
        self.manifest["graphpart_discarded_test_uniprots"] = sorted(self.discarded_test_uniprots)

        manifest_path = os.path.join(output_dir, "manifest.json")
        with open(manifest_path, "w") as mf:
            json.dump(self.manifest, mf, indent=2, sort_keys=True)
        print(f"Wrote {manifest_path}")


def _scan_source(jsonl_path: str) -> tuple[dict[str, str], Counter, Counter, int, int]:
    """First pass: collect uniprot→ref_sequence, record counts per uniprot,
    global qtype distribution. Returns (uniprot_to_seq, uniprot_counts,
    qtype_counts, total_records_used, total_records_scanned)."""
    uniprot_to_seq: dict[str, str] = {}
    uniprot_counts: Counter = Counter()
    qtype_counts: Counter = Counter()
    total_used = 0
    total_scanned = 0
    seq_conflicts = 0

    with open(jsonl_path) as f:
        for line in f:
            if not line.strip():
                continue
            total_scanned += 1
            rec = json.loads(line)
            if rec.get("grounding") != "cross-grounded":
                continue
            up = rec.get("uniprot_primary")
            seq = rec.get("ref_sequence")
            if up is None or seq is None:
                continue
            total_used += 1
            uniprot_counts[up] += 1
            qtype_counts[rec.get("question_type", "unknown")] += 1
            if up in uniprot_to_seq:
                if uniprot_to_seq[up] != seq:
                    seq_conflicts += 1
            else:
                uniprot_to_seq[up] = seq

    if seq_conflicts:
        print(f"WARN: {seq_conflicts} records had a ref_sequence differing from the first seen for the same uniprot (kept first).", file=sys.stderr)

    return uniprot_to_seq, uniprot_counts, qtype_counts, total_used, total_scanned


def _write_fasta(uniprot_to_seq: dict[str, str], fasta_path: str) -> None:
    with open(fasta_path, "w") as f:
        for up, seq in uniprot_to_seq.items():
            f.write(f">{up}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


def _run_mmseqs_cluster(
    mmseqs: str,
    fasta_path: str,
    out_prefix: str,
    tmp_dir: str,
    identity: float,
    coverage: float,
) -> str:
    cmd = [
        mmseqs, "easy-cluster",
        fasta_path, out_prefix, tmp_dir,
        "--min-seq-id", f"{identity}",
        "-c", f"{coverage}",
        "--cov-mode", "0",
    ]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"mmseqs easy-cluster failed ({result.returncode})")
    cluster_tsv = f"{out_prefix}_cluster.tsv"
    if not os.path.exists(cluster_tsv):
        raise RuntimeError(f"Expected {cluster_tsv} but not found")
    return cluster_tsv


def _parse_cluster_tsv(path: str) -> dict[str, list[str]]:
    clusters: dict[str, list[str]] = defaultdict(list)
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                continue
            rep, member = parts
            clusters[rep].append(member)
    return dict(clusters)


def _find_violating_test_uniprots(
    mmseqs: str,
    uniprot_to_seq: dict[str, str],
    uniprot_to_split: dict[str, str],
    identity_threshold: float,
    coverage: float,
    tmp_root: str,
) -> tuple[set[str], dict[str, float]]:
    """GraphPart-style audit: return the set of TEST uniprots with at least
    one TRAIN hit at >= identity_threshold (fraction in [0, 1]). Also returns
    a {test_uniprot: max_observed_pid_pct} map for telemetry.

    Implementation: mmseqs easy-search test-fasta against train-fasta with
    --min-seq-id 0.0 (all hits) and parse pident percent.
    """
    test_seqs = {up: seq for up, seq in uniprot_to_seq.items() if uniprot_to_split.get(up) == "test"}
    train_seqs = {up: seq for up, seq in uniprot_to_seq.items() if uniprot_to_split.get(up) == "train"}
    if not test_seqs or not train_seqs:
        return set(), {}

    threshold_pct = identity_threshold * 100.0
    violators: set[str] = set()
    max_pid: dict[str, float] = {}

    with tempfile.TemporaryDirectory(dir=tmp_root) as td:
        fa_test = os.path.join(td, "test.fasta")
        fa_train = os.path.join(td, "train.fasta")
        _write_fasta(test_seqs, fa_test)
        _write_fasta(train_seqs, fa_train)
        out_m8 = os.path.join(td, "test_vs_train.m8")
        mtmp = os.path.join(td, "mtmp")
        cmd = [
            mmseqs, "easy-search", fa_test, fa_train, out_m8, mtmp,
            "--min-seq-id", "0.0",
            "-c", f"{coverage}", "--cov-mode", "0",
            "-s", "7.5",
            "--max-seqs", "300",
            "--format-output", "query,target,pident",
        ]
        print(f"Running: {' '.join(cmd)}")
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr, file=sys.stderr)
            raise RuntimeError("mmseqs easy-search (test vs train, GraphPart audit) failed")
        if os.path.exists(out_m8):
            with open(out_m8) as f:
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue
                    q, t, pid = parts[0], parts[1], float(parts[2])
                    if q == t:
                        continue
                    if pid >= threshold_pct:
                        violators.add(q)
                        if pid > max_pid.get(q, 0.0):
                            max_pid[q] = pid

    return violators, max_pid


def _assign_clusters_to_splits(
    clusters: dict[str, list[str]],
    uniprot_counts: Counter,
    val_fraction: float,
    test_fraction: float,
    seed: int,
) -> tuple[dict[str, str], dict[str, int]]:
    """Seeded-shuffle cluster order, greedy-assign by record count to the
    split that is most underfilled relative to its target. Splits with a
    zero target (e.g. val_fraction=0 in an 80/20 setup) are excluded."""
    train_fraction = 1.0 - val_fraction - test_fraction
    total_records = sum(uniprot_counts.values())
    raw_targets = {
        "train": train_fraction * total_records,
        "val": val_fraction * total_records,
        "test": test_fraction * total_records,
    }
    targets = {k: v for k, v in raw_targets.items() if v > 0}
    current = {k: 0 for k in targets}

    cluster_sizes = {
        rep: sum(uniprot_counts.get(m, 0) for m in members)
        for rep, members in clusters.items()
    }

    rng = random.Random(seed)
    rep_order = list(clusters.keys())
    rng.shuffle(rep_order)

    cluster_to_split: dict[str, str] = {}
    for rep in rep_order:
        size = cluster_sizes[rep]
        # Split with the largest remaining capacity (target - current).
        pick = max(current.keys(), key=lambda s: targets[s] - current[s])
        cluster_to_split[rep] = pick
        current[pick] += size

    return cluster_to_split, current


def _max_cross_split_identity(
    mmseqs: str,
    uniprot_to_seq: dict[str, str],
    uniprot_to_split: dict[str, str],
    tmp_root: str,
) -> dict:
    """All-vs-all pairwise identity between test and train uniprots.
    Reports max observed identity for each pair direction."""
    by_split: dict[str, dict[str, str]] = {"train": {}, "val": {}, "test": {}}
    for up, seq in uniprot_to_seq.items():
        split = uniprot_to_split.get(up)
        if split is not None:
            by_split[split][up] = seq

    results: dict = {}
    for a, b in [("test", "train"), ("val", "train"), ("test", "val")]:
        if not by_split[a] or not by_split[b]:
            results[f"{a}_vs_{b}_max_identity"] = None
            continue
        with tempfile.TemporaryDirectory(dir=tmp_root) as td:
            fa_a = os.path.join(td, f"{a}.fasta")
            fa_b = os.path.join(td, f"{b}.fasta")
            _write_fasta(by_split[a], fa_a)
            _write_fasta(by_split[b], fa_b)
            out_m8 = os.path.join(td, "hits.m8")
            mtmp = os.path.join(td, "mtmp")
            cmd = [
                mmseqs, "easy-search", fa_a, fa_b, out_m8, mtmp,
                "--min-seq-id", "0.0",
                "-c", "0.8", "--cov-mode", "0",
                "-s", "7.5",
                "--max-seqs", "300",
                "--format-output", "query,target,pident",
            ]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                print(r.stderr, file=sys.stderr)
                raise RuntimeError(f"mmseqs easy-search ({a} vs {b}) failed")
            max_pid = 0.0
            worst_pair = None
            if os.path.exists(out_m8):
                with open(out_m8) as f:
                    for line in f:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 3:
                            continue
                        q, t, pid = parts[0], parts[1], float(parts[2])
                        if q == t:
                            continue
                        if pid > max_pid:
                            max_pid = pid
                            worst_pair = (q, t)
            results[f"{a}_vs_{b}_max_identity"] = max_pid
            results[f"{a}_vs_{b}_worst_pair"] = worst_pair
    return results


def build_split(
    jsonl_path: str = DEFAULT_JSONL,
    identity_threshold: float = 0.4,
    coverage: float = 0.8,
    temporal_cutoff_year: int | None = None,
    val_fraction: float = 0.0,
    test_fraction: float = 0.2,
    seed: int = 42,
    mmseqs_path: str = DEFAULT_MMSEQS,
    work_dir: str | None = None,
) -> Split:
    if temporal_cutoff_year is not None:
        raise NotImplementedError("temporal_cutoff_year is not implemented")

    if not os.path.exists(jsonl_path):
        raise FileNotFoundError(jsonl_path)
    if shutil.which(mmseqs_path) is None and not os.path.exists(mmseqs_path):
        raise RuntimeError(f"mmseqs not found at {mmseqs_path}")

    print(f"[1/6] Scanning {jsonl_path}")
    uniprot_to_seq, uniprot_counts, qtype_counts, total_used, total_scanned = _scan_source(jsonl_path)
    print(f"  scanned {total_scanned:,} records; kept {total_used:,} cross-grounded records with ref_sequence")
    print(f"  unique uniprots: {len(uniprot_to_seq):,}")

    tmp_root = work_dir or tempfile.mkdtemp(prefix="crossvalqa_splits_")
    os.makedirs(tmp_root, exist_ok=True)
    print(f"[2/6] Writing FASTA to {tmp_root}")
    fasta_path = os.path.join(tmp_root, "uniprots.fasta")
    _write_fasta(uniprot_to_seq, fasta_path)

    print(f"[3/6] Clustering at identity={identity_threshold} coverage={coverage}")
    out_prefix = os.path.join(tmp_root, "clusters")
    mmseqs_tmp = os.path.join(tmp_root, "mmseqs_tmp")
    os.makedirs(mmseqs_tmp, exist_ok=True)
    cluster_tsv = _run_mmseqs_cluster(
        mmseqs_path, fasta_path, out_prefix, mmseqs_tmp, identity_threshold, coverage
    )
    clusters = _parse_cluster_tsv(cluster_tsv)
    cluster_sizes = [len(members) for members in clusters.values()]
    print(f"  {len(clusters):,} clusters; sizes: min={min(cluster_sizes)} max={max(cluster_sizes)} mean={sum(cluster_sizes)/len(clusters):.2f}")

    print(f"[4/6] Assigning clusters to splits (val={val_fraction}, test={test_fraction}, seed={seed})")
    cluster_to_split, achieved = _assign_clusters_to_splits(
        clusters, uniprot_counts, val_fraction, test_fraction, seed
    )
    uniprot_to_split: dict[str, str] = {}
    for rep, members in clusters.items():
        split = cluster_to_split[rep]
        for m in members:
            uniprot_to_split[m] = split
    total = sum(achieved.values())
    for name, count in achieved.items():
        print(f"  {name}: {count:,} records ({count/total:.2%})")

    print(f"[5/6] GraphPart-style audit: discarding test uniprots with >= {identity_threshold:.0%} identity to any train uniprot")
    discarded, discard_max_pid = _find_violating_test_uniprots(
        mmseqs_path, uniprot_to_seq, uniprot_to_split, identity_threshold, coverage, tmp_root
    )
    print(f"  flagged {len(discarded):,} test uniprot(s) for discard (test side only)")
    if discarded:
        sample = sorted(discard_max_pid.items(), key=lambda kv: -kv[1])[:5]
        for up, pid in sample:
            print(f"    {up}: max identity to train = {pid:.1f}%")

    print("[6/6] Verifying max cross-split identity (post-discard)")
    surviving_uniprot_to_split = {
        up: s for up, s in uniprot_to_split.items() if up not in discarded
    }
    cross_identity = _max_cross_split_identity(
        mmseqs_path, uniprot_to_seq, surviving_uniprot_to_split, tmp_root
    )
    for k, v in cross_identity.items():
        print(f"  {k}: {v}")

    from collections import Counter as _Counter
    cluster_size_hist = _Counter(cluster_sizes)
    record_cluster_sizes = [
        sum(uniprot_counts.get(m, 0) for m in members)
        for members in clusters.values()
    ]

    manifest: dict = {
        "source_jsonl": jsonl_path,
        "params": {
            "identity_threshold": identity_threshold,
            "coverage": coverage,
            "val_fraction": val_fraction,
            "test_fraction": test_fraction,
            "seed": seed,
            "mmseqs_path": mmseqs_path,
            "mmseqs_version": _mmseqs_version(mmseqs_path),
            "temporal_cutoff_year": temporal_cutoff_year,
            "audit_method": "graphpart-style test-side discard at identity_threshold",
        },
        "source_stats": {
            "records_scanned": total_scanned,
            "records_used": total_used,
            "unique_uniprots": len(uniprot_to_seq),
            "question_type_counts": dict(qtype_counts),
        },
        "clusters": {
            "count": len(clusters),
            "protein_size_histogram": {
                str(k): v for k, v in sorted(cluster_size_hist.items())
            },
            "record_size_max": max(record_cluster_sizes),
            "record_size_mean": sum(record_cluster_sizes) / len(record_cluster_sizes),
        },
        "assignment_targets": {
            k: v for k, v in {
                "train": (1.0 - val_fraction - test_fraction),
                "val": val_fraction,
                "test": test_fraction,
            }.items() if v > 0
        },
        "assignment_achieved_record_counts": achieved,
        "graphpart_audit": {
            "discarded_test_uniprot_count": len(discarded),
            "discarded_test_uniprot_max_identities": discard_max_pid,
        },
        "cross_split_identity": cross_identity,
    }

    return Split(
        uniprot_to_split=uniprot_to_split,
        cluster_to_split=cluster_to_split,
        cluster_members=clusters,
        source_jsonl=jsonl_path,
        discarded_test_uniprots=discarded,
        manifest=manifest,
    )


def build_random_split(
    jsonl_path: str = DEFAULT_JSONL,
    val_fraction: float = 0.0,
    test_fraction: float = 0.2,
    seed: int = 42,
) -> Split:
    """Random uniprot-level train/test split with no homology awareness.

    Same source filter as build_split (cross-grounded + uniprot_primary +
    ref_sequence) and the same uniprot-level grouping (no uniprot in both
    splits). Each uniprot is its own singleton cluster, so paralogs and
    homologs can land on either side. Useful as a baseline for measuring
    the impact of homology-aware splitting.
    """
    if not os.path.exists(jsonl_path):
        raise FileNotFoundError(jsonl_path)

    print(f"[1/3] Scanning {jsonl_path}")
    uniprot_to_seq, uniprot_counts, qtype_counts, total_used, total_scanned = _scan_source(jsonl_path)
    print(f"  scanned {total_scanned:,} records; kept {total_used:,} cross-grounded records with ref_sequence")
    print(f"  unique uniprots: {len(uniprot_to_seq):,}")

    print(f"[2/3] Building singleton clusters and assigning (val={val_fraction}, test={test_fraction}, seed={seed})")
    clusters = {up: [up] for up in uniprot_to_seq}
    cluster_to_split, achieved = _assign_clusters_to_splits(
        clusters, uniprot_counts, val_fraction, test_fraction, seed
    )
    uniprot_to_split = {rep: cluster_to_split[rep] for rep in clusters}
    total = sum(achieved.values())
    for name, count in achieved.items():
        print(f"  {name}: {count:,} records ({count/total:.2%})")

    print("[3/3] Building manifest")
    manifest: dict = {
        "source_jsonl": jsonl_path,
        "params": {
            "split_type": "random_uniprot_level",
            "val_fraction": val_fraction,
            "test_fraction": test_fraction,
            "seed": seed,
        },
        "source_stats": {
            "records_scanned": total_scanned,
            "records_used": total_used,
            "unique_uniprots": len(uniprot_to_seq),
            "question_type_counts": dict(qtype_counts),
        },
        "assignment_targets": {
            k: v for k, v in {
                "train": (1.0 - val_fraction - test_fraction),
                "val": val_fraction,
                "test": test_fraction,
            }.items() if v > 0
        },
        "assignment_achieved_record_counts": achieved,
    }

    return Split(
        uniprot_to_split=uniprot_to_split,
        cluster_to_split=cluster_to_split,
        cluster_members=clusters,
        source_jsonl=jsonl_path,
        manifest=manifest,
    )


def _mmseqs_version(mmseqs_path: str) -> str:
    try:
        r = subprocess.run([mmseqs_path, "version"], capture_output=True, text=True)
        return r.stdout.strip()
    except Exception:
        return "unknown"


def verify_split(split: Split) -> dict:
    """Post-hoc leakage checks. Returns dict with pass/fail per rule."""
    report: dict = {"passed": True, "checks": {}}

    # No uniprot in more than one split (trivially true — map is a dict — but
    # verify by checking cluster integrity).
    seen: dict[str, str] = {}
    offenders = []
    for rep, members in split.cluster_members.items():
        rep_split = split.cluster_to_split[rep]
        for m in members:
            if m in seen and seen[m] != rep_split:
                offenders.append(m)
            seen[m] = rep_split
    report["checks"]["unique_uniprot_per_split"] = {
        "passed": len(offenders) == 0,
        "offenders": offenders[:20],
    }
    if offenders:
        report["passed"] = False

    # Max cross-split identity below threshold (homology-aware splits only).
    # The random split intentionally does no homology audit, so this check is
    # skipped when no identity_threshold / cross_split_identity is recorded.
    if "identity_threshold" in split.manifest.get("params", {}) and "cross_split_identity" in split.manifest:
        thr_pct = split.manifest["params"]["identity_threshold"] * 100.0
        cross = split.manifest["cross_split_identity"]
        observed_pct = {k: v for k, v in cross.items() if k.endswith("_max_identity")}
        bad = {k: v for k, v in observed_pct.items() if isinstance(v, (int, float)) and v >= thr_pct}
        report["checks"]["max_cross_split_identity_below_threshold"] = {
            "passed": not bad,
            "threshold_percent": thr_pct,
            "observed_percent": observed_pct,
        }
        if bad:
            report["passed"] = False

    return report


def main() -> int:
    """CLI entrypoint. Usage:
        python create_splits.py [homology|random|both] [--input ...] [--mmseqs ...]

    Default is `both` — produces homology_split/ and random_split/.
    """
    ap = argparse.ArgumentParser(description="Build the homology-aware (MMseqs2) and/or random uniprot-level splits.")
    ap.add_argument("mode", nargs="?", default="both", choices=("homology", "random", "both"))
    ap.add_argument("--input", default=DEFAULT_JSONL, help="release JSONL (default: release/output/crossvalqa_release_v3.jsonl)")
    ap.add_argument("--homology-out", default=DEFAULT_OUTPUT_DIR, help="output dir for the homology split (default: release/output/homology_split)")
    ap.add_argument("--random-out", default=DEFAULT_RANDOM_OUTPUT_DIR, help="output dir for the random split (default: release/output/random_split)")
    ap.add_argument("--mmseqs", default=DEFAULT_MMSEQS, help="MMseqs2 executable (default: mmseqs on PATH)")
    ap.add_argument("--work-dir", default=None, help="scratch dir for MMseqs2 (default: a new temp dir)")
    args = ap.parse_args()
    mode = args.mode

    all_passed = True
    if mode in ("homology", "both"):
        print("=" * 60)
        print("Building homology-aware split (MMseqs2 + GraphPart audit)")
        print("=" * 60)
        h = build_split(jsonl_path=args.input, mmseqs_path=args.mmseqs, work_dir=args.work_dir)
        h.write(args.homology_out)
        print("\nVerifying homology split integrity...")
        report = verify_split(h)
        print(json.dumps(report, indent=2))
        all_passed = all_passed and report["passed"]

    if mode in ("random", "both"):
        print("\n" + "=" * 60)
        print("Building random uniprot-level split (no homology awareness)")
        print("=" * 60)
        r = build_random_split(jsonl_path=args.input)
        r.write(args.random_out)
        print("\nVerifying random split integrity...")
        report = verify_split(r)
        print(json.dumps(report, indent=2))
        all_passed = all_passed and report["passed"]

    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
