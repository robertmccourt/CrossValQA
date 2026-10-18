"""Build _gold_rows_v2.jsonl: protein-verifier input derived from v2 + UniProt cache.

Input:  uniprot/output/m2tqa_gene_corrected_v2.jsonl (423,415 rows) and
        audit/output/uniprot_sequences.json (fetch_uniprot_sequences.py).
Output: audit/output/_gold_rows_v2.jsonl, row-aligned to the v2 file.

Example (from the repo root):
    python audit/build_gold_rows_v2.py

Field mapping:

Per row:
  protein, Hugo_Symbol  <- corrected_gene
  mutation              <- mutation
  HGVSp_Short           <- mutation (verifier falls back to `mutation` anyway)
  ref_sequence          <- uniprot_sequences.json[uniprot_primary]
  mutation_type         <- parse_hgvs(mutation)["kind"], title-cased
  qagen_answer_final    <- phase1_answer
  llm2_answer_final     <- phase2_answer
  question_type         <- question_type
  split                 <- "expansion" if row_idx >= 416188 else "v1"
  _row_id               <- row_idx (string, as the verifier expects)
  pubmed_id             <- pmid (kept for traceability)
  uniprot_primary       <- copied through

Skips no rows; missing sequences just leave ref_sequence="".
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "verifier"))
from mutation_verifier_dual import parse_hgvs  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
V2 = REPO_ROOT / "uniprot/output/m2tqa_gene_corrected_v2.jsonl"
SEQ_CACHE = REPO_ROOT / "audit/output/uniprot_sequences.json"
OUT = REPO_ROOT / "audit/output/_gold_rows_v2.jsonl"
EXPANSION_BOUNDARY = 416188  # rows >= this index in v2 come from the expansion batch


def main() -> None:
    global V2, SEQ_CACHE, OUT
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--v2", type=Path, default=V2,
                    help=f"v2 dataset JSONL (default: {V2.relative_to(REPO_ROOT)})")
    ap.add_argument("--seq_cache", type=Path, default=SEQ_CACHE,
                    help=f"sequence cache (default: {SEQ_CACHE.relative_to(REPO_ROOT)})")
    ap.add_argument("--output", type=Path, default=OUT,
                    help=f"output JSONL (default: {OUT.relative_to(REPO_ROOT)})")
    args = ap.parse_args()
    V2, SEQ_CACHE, OUT = args.v2, args.seq_cache, args.output
    OUT.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading sequence cache from {SEQ_CACHE.name}...")
    seq_cache: dict[str, str] = json.loads(SEQ_CACHE.read_text())
    print(f"Cached sequences: {len(seq_cache):,}")

    written = 0
    no_seq = 0
    parse_ok = 0
    expansion_count = 0
    print(f"Reading {V2.name} -> writing {OUT.name}")

    with V2.open() as fin, OUT.open("w") as fout:
        for i, line in enumerate(fin):
            d = json.loads(line)
            mutation = d.get("mutation") or ""
            parsed = parse_hgvs(mutation)
            mtype = parsed["kind"].title() if parsed else ""
            if parsed:
                parse_ok += 1

            acc = d.get("uniprot_primary") or ""
            ref_seq = seq_cache.get(acc, "")
            if not ref_seq:
                no_seq += 1

            corrected = d.get("corrected_gene") or d.get("gene") or ""
            split = "expansion" if i >= EXPANSION_BOUNDARY else "v1"
            if split == "expansion":
                expansion_count += 1

            row = {
                "_row_id": str(i),
                "pubmed_id": d.get("pmid"),
                "protein": corrected,
                "Hugo_Symbol": corrected,
                "mutation": mutation,
                "HGVSp_Short": mutation,
                "mutation_type": mtype,
                "ref_sequence": ref_seq,
                "qagen_answer_final": d.get("phase1_answer") or "",
                "llm2_answer_final": d.get("phase2_answer") or "",
                "question": d.get("question") or "",
                "question_type": d.get("question_type") or "",
                "split": split,
                "uniprot_primary": acc,
            }
            fout.write(json.dumps(row, ensure_ascii=False) + "\n")
            written += 1
            if (written % 50000) == 0:
                print(f"  written {written:,}", flush=True)

    print()
    print(f"Wrote {written:,} rows to {OUT}")
    print(f"  parse_hgvs ok: {parse_ok:,} ({parse_ok/written:.1%})")
    print(f"  no ref_sequence: {no_seq:,} ({no_seq/written:.1%})")
    print(f"  split=expansion: {expansion_count:,}")
    print(f"  output size: {OUT.stat().st_size:,} bytes")


if __name__ == "__main__":
    main()
