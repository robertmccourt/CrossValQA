"""DNA-side deterministic verifier for MutQA (both answers per row).

Mirror of `mutation_verifier_dual.py` but for DNA-notation rows. Uses
`dna_facts.mutation_facts_dna()` to produce a facts dict shaped like the
protein verifier's, then runs:

  DNA-specific checks (always):
    - inline c.* HGVS notation  (`_check_dna_hgvs_inline`)
    - inline p.* protein derivation  (`_check_protein_derivation`)
    - splice-region claims  (`_check_splice_region`)
    - mutation-kind language  (`_check_dna_muttype`)
    - transition / transversion  (`_check_transition`)

  Reused protein checks (only when the substitution translates to a
  clean wt → mut AA pair, i.e. exonic non-synonymous substitution):
    - position, protein_length, charge, polarity, hydrophobicity, size,
      special-residue, region (from `mutation_verifier_dual`)

Aggregate label per row: `verified` / `contradicted` / `not_applicable` /
`parse_failed`. Answers are read from `phase1_answer` (qagen) and
`phase2_answer` (llm2).

Input:  dataset JSONL (the 411,060-row working set
        release/output/m2tqa_gene_corrected_v2_verified_answerable.jsonl);
        non-DNA-notation rows are skipped. `row_id` is the 0-based line index.
Output: <out_dir>/mutation_verification_dna.jsonl and
        <out_dir>/mutation_verification_dna_summary.md.
Needs:  pysam, the GRCh38 FASTA and the UniProt-hg38 CDS tables under
        $CROSSVALQA_REF_DIR (default reference_data/ at the repo root; see audit/README.md).

Example (from the repo root):
    python audit/verifier/mutation_verifier_dna.py \
        --input release/output/m2tqa_gene_corrected_v2_verified_answerable.jsonl \
        --out_dir audit/output/dna_v3 \
        --features audit/output/uniprot_features.json
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import triggers as T
import dna_triggers as D
from dna_facts import mutation_facts_dna
from mutation_verifier_dual import (
    _is_hypothetical,
    _has_negation_before,
    _has_pre_disqualifier,
    _has_post_disqualifier,
    _dedup_checks,
    _check_charge,
    _check_polarity,
    _check_hydrophobicity,
    _check_size,
    _check_special,
    _check_region,
    aggregate_label,
    iterate_jsonl,
)
from topic_bucket import bucket_topic
from uniprot_features import CACHE_PATH as FEATURES_CACHE, FeatureStore, load_features
from feature_checks import run_uniprot_feature_checks


# ────────────────────────────────────────────── DNA-notation detection

DNA_PREFIX = re.compile(r"^[cgmn]\.", re.IGNORECASE)
BARE_NUC = re.compile(r"^\d+[ACGT]>[ACGT]$|^[ACGT]\d+[ACGT]$")


def is_dna_notation(s):
    if not s:
        return False
    return bool(DNA_PREFIX.match(s) or BARE_NUC.match(s))


# ─────────────────────────────────────────────── DNA-specific checks

def _check_dna_hgvs_inline(answer, facts):
    """Find inline c.\\d+X>Y / bare X\\d+Y / \\d+X>Y mentions and verify
    pos / ref / alt against `facts`."""
    out = []
    cds_pos = facts.get("cds_position")
    ref, alt = facts.get("ref_base"), facts.get("alt_base")
    if not cds_pos:
        return out

    seen = set()  # (pos, ref, alt) tuples we've already emitted on
    for regex, fields in [
        (D.DNA_HGVS_SUB, ("pos", "off", "ref", "alt")),
        (D.DNA_BARE_SUB, ("pos", "ref", "alt")),
        (D.DNA_PRE_SUB, ("ref", "pos", "alt")),
    ]:
        for m in regex.finditer(answer):
            if _is_hypothetical(answer, m.start()):
                continue
            groups = m.groups()
            if regex is D.DNA_HGVS_SUB:
                pos = int(groups[0]); off = groups[1]; r = groups[2]; a = groups[3]
            elif regex is D.DNA_BARE_SUB:
                pos = int(groups[0]); off = None; r = groups[1]; a = groups[2]
            else:  # DNA_PRE_SUB
                pos = int(groups[1]); off = None; r = groups[0]; a = groups[2]
            if pos != cds_pos:
                continue  # mention is about a different position; not our claim
            key = (pos, r.upper(), a.upper())
            if key in seen:
                continue
            seen.add(key)

            # check ref base
            if ref:
                status = "verified" if r.upper() == ref.upper() else "contradicted"
                out.append({"claim_type": "dna_ref_base",
                            "claimed": r.upper(), "ref": ref.upper(),
                            "status": status})
            # check alt base
            if alt:
                status = "verified" if a.upper() == alt.upper() else "contradicted"
                out.append({"claim_type": "dna_alt_base",
                            "claimed": a.upper(), "ref": alt.upper(),
                            "status": status})
    return out


def _check_protein_derivation(answer, facts):
    """When the LLM gives an inline `p.Gly331Arg` / `G331R` derivation, verify
    against the codon-translated prediction. Only fires when we have a
    predicted wt/mut/position from the codon table."""
    out = []
    pred_wt = facts.get("wt_aa")
    pred_mut = facts.get("mut_aa")
    pred_pos = facts.get("position")
    if not (pred_wt and pred_mut and pred_pos):
        return out
    # Don't mark LLM contradicted on protein derivation when our reference
    # base disagrees with GRCh38 — usually means a transcript mismatch and
    # our predicted protein change is unreliable.
    if not facts.get("ref_matches_genome", True):
        return out

    seen = set()
    for regex in (D.PROT_THREE, D.PROT_ONE):
        for m in regex.finditer(answer):
            if _is_hypothetical(answer, m.start()):
                continue
            wt_tok, pos_tok, mut_tok = m.group(1), m.group(2), m.group(3)
            wt = D.norm_aa_token(wt_tok)
            mut = D.norm_aa_token(mut_tok)
            try:
                pos = int(pos_tok)
            except ValueError:
                continue
            if pos != pred_pos:
                continue  # different protein-position claim; not ours
            key = (wt, pos, mut)
            if key in seen:
                continue
            seen.add(key)

            if wt:
                status = "verified" if wt == pred_wt else "contradicted"
                out.append({"claim_type": "predicted_wt_aa",
                            "claimed": wt, "ref": pred_wt, "status": status})
            if mut:
                status = "verified" if mut == pred_mut else "contradicted"
                out.append({"claim_type": "predicted_mut_aa",
                            "claimed": mut, "ref": pred_mut, "status": status})
    return out


def _check_splice_region(answer, facts):
    """Verify splice-region claims only when the variant is at an intronic
    offset. Exonic mutations CAN affect splicing through enhancer disruption
    or cryptic-donor creation, so an exonic answer that mentions splicing
    isn't necessarily wrong — those fire as `not_applicable` (no row emitted)
    rather than contradicted."""
    out = []
    offset = facts.get("intronic_offset")
    in_splice = facts.get("splice_region", False)

    # Generic splice language — verified only when actually intronic.
    m = D.SPLICE_GENERIC_REGEX.search(answer)
    if m and not _is_hypothetical(answer, m.start()):
        if offset is not None:
            out.append({"claim_type": "splice_region",
                        "claimed": True, "ref": True, "status": "verified"})
        # else: exonic but answer mentions splicing — out of deterministic
        # scope, don't fire.

    # Donor specifics: intronic offset > 0 (5' end of intron after the exon).
    m = D.SPLICE_DONOR_REGEX.search(answer)
    if m and not _is_hypothetical(answer, m.start()):
        if in_splice and offset is not None and offset > 0:
            out.append({"claim_type": "splice_donor",
                        "claimed": True, "ref": True, "status": "verified"})
        # else: not at a donor offset — don't fire

    m = D.SPLICE_ACCEPTOR_REGEX.search(answer)
    if m and not _is_hypothetical(answer, m.start()):
        if in_splice and offset is not None and offset < 0:
            out.append({"claim_type": "splice_acceptor",
                        "claimed": True, "ref": True, "status": "verified"})
        # else: not at an acceptor offset — don't fire

    return out


def _check_dna_muttype(answer, facts):
    """Verify mutation-kind language against the derived kind.

    For exonic substitutions we know whether it's missense / synonymous /
    nonsense from the codon translation. For other kinds (del / ins / dup /
    delins / intronic-sub), we only check that the LLM doesn't explicitly
    misclassify it (e.g. calling a deletion 'missense' is a contradiction).
    """
    out = []
    kind = facts.get("kind", "")
    is_intronic = facts.get("intronic_offset") is not None

    # Step 1: synonymous claim — verified iff our kind is synonymous
    m = D.SYNONYMOUS_REGEX.search(answer)
    if m and not _is_hypothetical(answer, m.start()) \
            and not _has_negation_before(answer, m.start()):
        status = "verified" if kind == "synonymous" else "contradicted"
        out.append({"claim_type": "dna_muttype",
                    "claimed": "synonymous", "ref": kind, "status": status})

    # Step 2: missense claim — verified for substitution kind, contradicted for
    # synonymous / nonsense / del / ins / dup / delins
    m = D.MISSENSE_REGEX.search(answer)
    if m and not _is_hypothetical(answer, m.start()) \
            and not _has_negation_before(answer, m.start()):
        status = "verified" if kind == "substitution" else "contradicted"
        out.append({"claim_type": "dna_muttype",
                    "claimed": "missense", "ref": kind, "status": status})

    # Step 3a: strict nonsense-kind claim — verified iff kind is nonsense.
    # Only marks contradicted when the actual kind is a substitution-class
    # that we can confidently distinguish (synonymous, missense substitution).
    # For indels we can't know without checking length; treat as not_applicable.
    m = D.NONSENSE_KIND_REGEX.search(answer)
    if m and not _is_hypothetical(answer, m.start()) \
            and not _has_negation_before(answer, m.start()):
        if kind == "nonsense":
            out.append({"claim_type": "dna_muttype",
                        "claimed": "nonsense_kind", "ref": kind,
                        "status": "verified"})
        elif kind in ("substitution", "synonymous"):
            out.append({"claim_type": "dna_muttype",
                        "claimed": "nonsense_kind", "ref": kind,
                        "status": "contradicted"})
        # else: indel — out of scope without a length check

    # Step 3b: premature-stop / truncation language — describes the protein
    # consequence, not the kind directly. Verified for nonsense substitutions
    # AND for indels (frameshift indels typically produce premature stops).
    # Contradicted only for clean substitutions (which can't truncate via
    # frameshift) or synonymous.
    m = D.PREMATURE_STOP_REGEX.search(answer)
    if m and not _is_hypothetical(answer, m.start()) \
            and not _has_negation_before(answer, m.start()):
        if kind == "nonsense" or kind in ("del", "ins", "dup", "delins"):
            out.append({"claim_type": "premature_stop_consequence",
                        "claimed": True, "ref": True, "status": "verified"})
        elif kind in ("substitution", "synonymous"):
            out.append({"claim_type": "premature_stop_consequence",
                        "claimed": True, "ref": False, "status": "contradicted"})

    # Step 4: frameshift — verified for ins/del that change frame, contradicted
    # for substitution / synonymous / nonsense (single-base subs never shift
    # frame). Note: del/dup of multiples of 3 are in-frame, but we don't model
    # length here, so we treat any del/ins/dup/delins frameshift claim as
    # `not_applicable` — too noisy to mark contradicted without checking length.
    m = D.FRAMESHIFT_REGEX.search(answer)
    if m and not _is_hypothetical(answer, m.start()) \
            and not _has_negation_before(answer, m.start()):
        if kind in ("substitution", "synonymous", "nonsense"):
            out.append({"claim_type": "dna_muttype",
                        "claimed": "frameshift", "ref": kind,
                        "status": "contradicted"})
        # else: indel — frameshift status depends on length (out of scope)

    return out


def _check_transition(answer, facts):
    """Verify transition / transversion claims against ref/alt nucleotide
    purine/pyrimidine class."""
    out = []
    ref, alt = facts.get("ref_base"), facts.get("alt_base")
    is_trans = D.is_transition(ref, alt)
    if is_trans is None:
        return out

    m = D.TRANSITION_REGEX.search(answer)
    if m and not _is_hypothetical(answer, m.start()) \
            and not _has_negation_before(answer, m.start()):
        status = "verified" if is_trans else "contradicted"
        out.append({"claim_type": "nucleotide_class",
                    "claimed": "transition", "ref": "transition" if is_trans
                    else "transversion", "status": status})
    m = D.TRANSVERSION_REGEX.search(answer)
    if m and not _is_hypothetical(answer, m.start()) \
            and not _has_negation_before(answer, m.start()):
        status = "verified" if not is_trans else "contradicted"
        out.append({"claim_type": "nucleotide_class",
                    "claimed": "transversion", "ref": "transition" if is_trans
                    else "transversion", "status": status})
    return out


# ─────────────────────────────────────────────── orchestration

def verify_qa_dna(answer, facts, feat_store=None):
    if not answer or not facts:
        return []
    checks = []
    checks += _check_dna_hgvs_inline(answer, facts)
    checks += _check_protein_derivation(answer, facts)
    checks += _check_splice_region(answer, facts)
    checks += _check_dna_muttype(answer, facts)
    checks += _check_transition(answer, facts)

    # Reused protein-side AA-chemistry checks — only fire when the codon
    # translation gave us a clean substitution.
    if facts.get("change", {}).get("applicable"):
        checks += _check_charge(answer, facts)
        checks += _check_polarity(answer, facts)
        checks += _check_hydrophobicity(answer, facts)
        checks += _check_size(answer, facts)
        checks += _check_special(answer, facts)
        checks += _check_region(answer, facts)

    # UniProt feature-aware checks — fire when we have a clean
    # gene/transcript resolution. Skip when ref_matches_genome=False because
    # those rows have transcript / gene-correction issues that make the
    # position index unreliable.
    if (feat_store is not None
            and facts.get("position")
            and facts.get("ref_matches_genome", True)):
        checks += run_uniprot_feature_checks(answer, facts, feat_store)
    return _dedup_checks(checks)


def _resolve_gene(row):
    return (row.get("corrected_gene") or row.get("recovered_gene")
            or row.get("gene") or "")


ANSWER_SOURCES = ("qagen", "llm2")
ANSWER_FIELDS = {"qagen": "phase1_answer", "llm2": "phase2_answer"}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--input", required=True,
                   help="JSONL with v2 rows (e.g. m2tqa_gene_corrected_v2_"
                        "verified_answerable.jsonl)")
    p.add_argument("--out_dir", required=True)
    p.add_argument("--limit", type=int, default=0,
                   help="Stop after N DNA-notation rows (0 = no limit)")
    p.add_argument("--all_rows", action="store_true",
                   help="Don't skip non-DNA rows; verify any row whose "
                        "mutation field looks DNA-like")
    p.add_argument("--features", type=Path, default=FEATURES_CACHE,
                   help="UniProt feature cache from fetch_uniprot_features.py "
                        "(default: audit/output/uniprot_features.json). If "
                        "missing, the uniprot_* checks are skipped.")
    p.add_argument("--no_features", action="store_true",
                   help="Skip the UniProt-feature checks (uniprot_* claim "
                        "families) even if the cache exists.")
    args = p.parse_args()

    src = Path(args.input)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    per_qa_path = out_dir / "mutation_verification_dna.jsonl"

    # Load UniProt features once at startup. None disables the new checks.
    try:
        if args.no_features:
            raise FileNotFoundError("--no_features given.")
        feat_store = FeatureStore(load_features(args.features))
        print(f"Loaded UniProt features for {sum(1 for _ in feat_store._raw):,} accessions")
    except FileNotFoundError as e:
        print(f"WARNING: {e}\nProceeding without UniProt feature checks.")
        feat_store = None

    bucket_totals = {
        s: defaultdict(
            lambda: Counter({"checkable_qa": 0, "total_qa": 0,
                             "verified_claims": 0, "contradicted_claims": 0,
                             "total_claims": 0})
        )
        for s in ANSWER_SOURCES
    }
    by_claim_type = {s: Counter() for s in ANSWER_SOURCES}
    by_claim_type_verified = {s: Counter() for s in ANSWER_SOURCES}

    n_dna_seen = 0
    n_resolved = 0
    n_parse_fail = 0
    n_ref_mismatch = 0

    with per_qa_path.open("w") as fout:
        for i, row in enumerate(iterate_jsonl(src)):
            mut = row.get("mutation", "")
            if not is_dna_notation(mut):
                continue
            n_dna_seen += 1
            if args.limit and n_dna_seen > args.limit:
                break

            gene = _resolve_gene(row)
            facts = mutation_facts_dna(gene, mut)
            qtype = row.get("question_type")
            bucket = bucket_topic(qtype) if qtype else "unknown"
            split = row.get("split") or "unknown"

            if facts is None:
                n_parse_fail += 1
                fout.write(json.dumps({
                    "row_id": i, "gene": gene, "mutation": mut,
                    "parse_failed": True,
                }) + "\n")
                for s in ANSWER_SOURCES:
                    bucket_totals[s][(bucket, split)]["total_qa"] += 1
                continue

            n_resolved += 1
            if not facts["ref_matches_genome"]:
                n_ref_mismatch += 1

            results = {}
            for s in ANSWER_SOURCES:
                ans = row.get(ANSWER_FIELDS[s], "") or ""
                ck = verify_qa_dna(ans, facts, feat_store)
                results[s] = {"checks": ck, "label": aggregate_label(ck)}

            out_row = {
                "row_id": i, "gene": gene, "mutation": mut,
                "question_type": qtype, "bucket": bucket, "split": split,
                "facts": {k: v for k, v in facts.items()
                          if k != "change" and not k.startswith("_")},
                "predicted_protein_change": (
                    f"{facts['wt_aa']}{facts['position']}{facts['mut_aa']}"
                    if facts["wt_aa"] and facts["position"] else None
                ),
                "ref_matches_genome": facts["ref_matches_genome"],
            }
            for s in ANSWER_SOURCES:
                out_row[f"{s}_checks"] = results[s]["checks"]
                out_row[f"{s}_label"] = results[s]["label"]
            fout.write(json.dumps(out_row, ensure_ascii=False) + "\n")

            for s in ANSWER_SOURCES:
                checks = results[s]["checks"]
                btot = bucket_totals[s][(bucket, split)]
                btot["total_qa"] += 1
                if checks:
                    btot["checkable_qa"] += 1
                    btot["total_claims"] += len(checks)
                    for c in checks:
                        by_claim_type[s][c["claim_type"]] += 1
                        if c["status"] == "verified":
                            btot["verified_claims"] += 1
                            by_claim_type_verified[s][c["claim_type"]] += 1
                        else:
                            btot["contradicted_claims"] += 1

            if n_dna_seen % 5000 == 0:
                print(f"  processed {n_dna_seen:,} DNA rows "
                      f"({n_resolved:,} resolved, {n_parse_fail:,} failed)",
                      flush=True)

    summary_path = out_dir / "mutation_verification_dna_summary.md"
    lines = []
    lines.append("# DNA-Coding Deterministic Verification — Summary\n")
    lines.append(
        f"Source: `{src}`  \nDNA rows seen: **{n_dna_seen:,}**  "
        f"(resolved: {n_resolved:,}; parse/gene failures: {n_parse_fail:,}; "
        f"ref-base mismatches with GRCh38: {n_ref_mismatch:,})\n"
        f"Answer sources: `qagen`={ANSWER_FIELDS['qagen']}, "
        f"`llm2`={ANSWER_FIELDS['llm2']}\n"
    )

    def _fmt(vals):
        total = vals["total_qa"]; checkable = vals["checkable_qa"]
        verified = vals["verified_claims"]; contradicted = vals["contradicted_claims"]
        claims = vals["total_claims"]
        cov = checkable / total if total else 0
        ver = verified / claims if claims else 0
        con = contradicted / claims if claims else 0
        return (f"{total:,}", f"{checkable:,}", f"{cov:.3%}",
                f"{claims:,}", f"{ver:.3%}", f"{con:.3%}")

    for s in ANSWER_SOURCES:
        lines.append(f"\n---\n\n# Source: `{s}` ({ANSWER_FIELDS[s]})\n")
        overall = Counter()
        for v in bucket_totals[s].values():
            overall.update(v)
        lines.append("## Overall\n")
        lines.append("| total Q&A | checkable | coverage | claims | verified% | contradicted% |")
        lines.append("|---|---|---|---|---|---|")
        lines.append("| " + " | ".join(_fmt(overall)) + " |\n")

        by_bucket = defaultdict(Counter)
        for (bucket, split), v in bucket_totals[s].items():
            by_bucket[bucket].update(v)
        lines.append("## By bucket\n")
        lines.append("| bucket | total Q&A | checkable | coverage | claims | verified% | contradicted% |")
        lines.append("|---|---|---|---|---|---|---|")
        for bucket in sorted(by_bucket.keys()):
            lines.append(f"| {bucket} | " + " | ".join(_fmt(by_bucket[bucket])) + " |")

        lines.append("\n## Per claim type\n")
        lines.append("| claim_type | total fired | verified | verified% |")
        lines.append("|---|---|---|---|")
        for ct in sorted(by_claim_type[s].keys()):
            tot = by_claim_type[s][ct]
            ver = by_claim_type_verified[s][ct]
            rate = ver / tot if tot else 0
            lines.append(f"| {ct} | {tot:,} | {ver:,} | {rate:.3%} |")

    summary_path.write_text("\n".join(lines))
    print(f"\nWrote:\n  {per_qa_path}\n  {summary_path}")
    print(f"  DNA seen: {n_dna_seen:,}  resolved: {n_resolved:,}  "
          f"parse fail: {n_parse_fail:,}  ref mismatch: {n_ref_mismatch:,}")


if __name__ == "__main__":
    main()
