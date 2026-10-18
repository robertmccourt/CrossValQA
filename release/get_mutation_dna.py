#!/usr/bin/env python3
"""
Get a DNA sequence window around a coding DNA mutation in the human genome (GRCh38).

Supports substitutions, deletions, insertions, duplications, delins,
and intronic/splice-site variants.

Standalone helper for dataset users (not a step of the release build). Needs
pysam and the GRCh38 / UniProt CDS reference files listed below.

Usage:
    python get_mutation_dna.py BRAF  "c.1799T>A"        --window 60
    python get_mutation_dna.py KRAS  "35G>A"             --window 40
    python get_mutation_dna.py BRCA1 "c.5266dupC"        --window 60
    python get_mutation_dna.py CFTR  "c.1521_1523delCTT" --window 60
    python get_mutation_dna.py BRCA2 "c.316+1G>A"        --window 40

As a library:
    from get_mutation_dna import get_mutation_dna_window
    result = get_mutation_dna_window("BRAF", "c.1799T>A", window=60)
"""

import os
import re
import sys
import argparse
import pysam
from pathlib import Path

# ── Reference-data paths ───────────────────────────────────────────────────────
# Override the root via: export CROSSVALQA_REF_DIR=/path/to/reference_root
# The directory must contain:
#   reference_genome/GRCh38.primary_assembly.genome.fa (+ .fai index)
#   uniprot/hg38_uniprot_sequences.csv
#   uniprot/uniprotGtfHg38_swiss.csv
#   uniprot/uniprotGtfHg38_trembl.csv

# Default: <repo>/reference_data (gitignored data; not shipped with the repo).
_REF_DIR = Path(os.environ.get("CROSSVALQA_REF_DIR",
                               Path(__file__).resolve().parents[1] / "reference_data"))
GENOME_FA = _REF_DIR / "reference_genome/GRCh38.primary_assembly.genome.fa"
UNIPROT_SEQ = _REF_DIR / "uniprot/hg38_uniprot_sequences.csv"
GTF_SWISS = _REF_DIR / "uniprot/uniprotGtfHg38_swiss.csv"
GTF_TREMBL = _REF_DIR / "uniprot/uniprotGtfHg38_trembl.csv"


def _check_ref_data():
    missing = [p for p in (GENOME_FA, UNIPROT_SEQ, GTF_SWISS, GTF_TREMBL) if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Reference data not found:\n"
            + "\n".join(f"  {p}" for p in missing)
            + f"\n\nSet CROSSVALQA_REF_DIR to the directory containing these files "
            f"(currently: {_REF_DIR})."
        )

# ── Helpers ────────────────────────────────────────────────────────────────────

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def _complement(seq):
    return seq.translate(_COMP)


def _revcomp(seq):
    return seq.translate(_COMP)[::-1]


# ── Lazy-loaded data ───────────────────────────────────────────────────────────

_cache = {}


def _get_genome():
    if "genome" not in _cache:
        _check_ref_data()
        _cache["genome"] = pysam.FastaFile(str(GENOME_FA))
    return _cache["genome"]


def _get_gene_map():
    """gene_name (upper) → {uid, chrom, strand}  (canonical entry)."""
    if "genes" in _cache:
        return _cache["genes"]

    genes = {}
    with open(UNIPROT_SEQ) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {h.lstrip("#"): i for i, h in enumerate(header)}
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            gene = cols[ci["geneName"]].strip()
            if not gene:
                continue
            is_primary = cols[ci["isMain"]] == "primary sequence"
            gu = gene.upper()
            if gu in genes and genes[gu]["is_primary"] and not is_primary:
                continue
            rec = {
                "uid": cols[ci["name"]],
                "chrom": cols[ci["chrom"]],
                "strand": cols[ci["strand"]],
                "is_primary": is_primary,
            }
            genes[gu] = rec
            for syn in re.split(r"[;,]", cols[ci["geneSynonyms"]]):
                su = syn.strip().upper()
                if su and not (su in genes and genes[su]["is_primary"]):
                    genes[su] = rec

    _cache["genes"] = genes
    return genes


def _get_cds_map():
    """uid → list of CDS exon dicts [{chrom, start, end, strand}]."""
    if "cds" in _cache:
        return _cache["cds"]

    cds = {}
    swiss_ids = set()
    for gtf_path, is_swiss in [(GTF_SWISS, True), (GTF_TREMBL, False)]:
        with open(gtf_path) as fh:
            for line in fh:
                cols = line.split("\t")
                if len(cols) < 9 or cols[2] != "CDS":
                    continue
                m = re.search(r'gene_id "([^"]+)"', cols[8])
                if not m:
                    continue
                uid = m.group(1)
                if not is_swiss and uid in swiss_ids:
                    continue
                if is_swiss:
                    swiss_ids.add(uid)
                cds.setdefault(uid, []).append({
                    "start": int(cols[3]),  # 1-based inclusive (GTF)
                    "end": int(cols[4]),
                })

    _cache["cds"] = cds
    return cds


# ── CDS base map (cached per transcript) ──────────────────────────────────────

_cds_base_cache = {}


def _get_cds_bases(uid, cds_exons, strand):
    """Return list of 0-based genomic positions in coding order (5'→3')."""
    if uid in _cds_base_cache:
        return _cds_base_cache[uid]

    if strand == "+":
        exons = sorted(cds_exons, key=lambda e: e["start"])
        bases = []
        for ex in exons:
            bases.extend(range(ex["start"] - 1, ex["end"]))
    else:
        exons = sorted(cds_exons, key=lambda e: e["start"], reverse=True)
        bases = []
        for ex in exons:
            bases.extend(range(ex["end"] - 1, ex["start"] - 2, -1))

    _cds_base_cache[uid] = bases
    return bases


# ── Parsing ────────────────────────────────────────────────────────────────────

def parse_dna_mutation(mutation_str):
    """Parse a DNA mutation string into a structured dict.

    Returns dict with keys:
        type:       "sub" | "del" | "ins" | "dup" | "delins"
        pos:        primary CDS position (int, 1-based)
        end_pos:    end CDS position for ranges (int or None)
        offset:     intronic offset at pos (int or None, e.g. +1, -12)
        end_offset: intronic offset at end_pos (int or None)
        ref_bases:  specified reference bases (str or None)
        alt_bases:  specified alt/inserted bases (str or None)
        raw:        original input string
    """
    raw = mutation_str.strip()
    s = raw

    # strip c. prefix
    if s.lower().startswith("c."):
        s = s[2:]

    # normalise arrows: T→C  →  T>C
    s = s.replace("\u2192", ">").replace("→", ">")

    result = {"raw": raw, "end_pos": None, "offset": None, "end_offset": None,
              "ref_bases": None, "alt_bases": None}

    # ── Substitutions ──────────────────────────────────────────────────────

    # c.1799T>A  or  c.316+1G>A
    m = re.match(r"^(\d+)([+-]\d+)?([ACGTacgt])>([ACGTacgt])$", s)
    if m:
        result.update(type="sub", pos=int(m.group(1)),
                      offset=int(m.group(2)) if m.group(2) else None,
                      ref_bases=m.group(3).upper(), alt_bases=m.group(4).upper())
        return result

    # T1799A  (bare ref-pos-alt)
    m = re.match(r"^([ACGTacgt])(\d+)([ACGTacgt])$", s)
    if m:
        result.update(type="sub", pos=int(m.group(2)),
                      ref_bases=m.group(1).upper(), alt_bases=m.group(3).upper())
        return result

    # Old style: C538T  (same as 538C>T but letter first with c. prefix)
    # Already caught by bare ref-pos-alt above if no c. prefix;
    # with c. prefix the s has "C538T" after stripping
    # This is identical to the T1799A pattern so it's already handled.

    # ── Deletions ──────────────────────────────────────────────────────────

    # range: c.241_242delCT  or  c.1355_1356del
    m = re.match(r"^(\d+)([+-]\d+)?_(\d+)([+-]\d+)?del([ACGTacgt]*)$", s, re.I)
    if m:
        result.update(type="del", pos=int(m.group(1)), end_pos=int(m.group(3)),
                      offset=int(m.group(2)) if m.group(2) else None,
                      end_offset=int(m.group(4)) if m.group(4) else None,
                      ref_bases=m.group(5).upper() or None)
        return result

    # single: c.210delA  or  c.210del  or  c.226+2delT
    m = re.match(r"^(\d+)([+-]\d+)?del([ACGTacgt]*)$", s, re.I)
    if m:
        result.update(type="del", pos=int(m.group(1)),
                      offset=int(m.group(2)) if m.group(2) else None,
                      ref_bases=m.group(3).upper() or None)
        return result

    # ── Duplications ───────────────────────────────────────────────────────

    # range: c.123_125dup  or  c.123_125dupACG
    m = re.match(r"^(\d+)([+-]\d+)?_(\d+)([+-]\d+)?dup([ACGTacgt]*)$", s, re.I)
    if m:
        result.update(type="dup", pos=int(m.group(1)), end_pos=int(m.group(3)),
                      offset=int(m.group(2)) if m.group(2) else None,
                      end_offset=int(m.group(4)) if m.group(4) else None,
                      ref_bases=m.group(5).upper() or None)
        return result

    # single: c.2527dupG  or  c.2527dup
    m = re.match(r"^(\d+)([+-]\d+)?dup([ACGTacgt]*)$", s, re.I)
    if m:
        result.update(type="dup", pos=int(m.group(1)),
                      offset=int(m.group(2)) if m.group(2) else None,
                      ref_bases=m.group(3).upper() or None)
        return result

    # ── Insertions ─────────────────────────────────────────────────────────

    # range: c.123_124insACG
    m = re.match(r"^(\d+)([+-]\d+)?_(\d+)([+-]\d+)?ins([ACGTacgt]+)$", s, re.I)
    if m:
        result.update(type="ins", pos=int(m.group(1)), end_pos=int(m.group(3)),
                      offset=int(m.group(2)) if m.group(2) else None,
                      end_offset=int(m.group(4)) if m.group(4) else None,
                      alt_bases=m.group(5).upper())
        return result

    # non-standard single-pos: c.640insC
    m = re.match(r"^(\d+)ins([ACGTacgt]+)$", s, re.I)
    if m:
        p = int(m.group(1))
        result.update(type="ins", pos=p, end_pos=p + 1,
                      alt_bases=m.group(2).upper())
        return result

    # ── Delins ─────────────────────────────────────────────────────────────

    # c.123_125delinsACG  or  c.123delinsACG
    m = re.match(r"^(\d+)(_(\d+))?delins([ACGTacgt]+)$", s, re.I)
    if m:
        result.update(type="delins", pos=int(m.group(1)),
                      end_pos=int(m.group(3)) if m.group(3) else None,
                      alt_bases=m.group(4).upper())
        return result

    raise ValueError(f"Could not parse: {raw}")


# ── Position resolution ────────────────────────────────────────────────────────

def _cds_to_genomic(cds_pos, cds_bases, strand, intronic_offset=None):
    """Map 1-based CDS position (+ optional intronic offset) to 0-based genomic pos."""
    if cds_pos < 1 or cds_pos > len(cds_bases):
        raise ValueError(
            f"CDS position {cds_pos} out of range (CDS length {len(cds_bases)} bp)"
        )
    gpos = cds_bases[cds_pos - 1]
    if intronic_offset:
        # HGVS offset is in coding direction; on - strand that's reversed genomically
        gpos += intronic_offset if strand == "+" else -intronic_offset
    return gpos


def _fetch_coding_bases(gpositions, chrom, strand, genome):
    """Fetch bases at genomic positions and return on the coding strand."""
    plus_bases = "".join(genome.fetch(chrom, p, p + 1) for p in gpositions).upper()
    return _complement(plus_bases) if strand == "-" else plus_bases


# ── Core ───────────────────────────────────────────────────────────────────────

def get_mutation_dna_window(gene, mutation_str, window=100):
    """Look up a coding-DNA mutation and return a reference-genome DNA window.

    Parameters
    ----------
    gene : str          Gene symbol, e.g. "TP53"
    mutation_str : str  DNA mutation, e.g. "c.1799T>A", "c.241_242delCT"
    window : int        Total window size in bp, centered on the mutation site.

    Returns
    -------
    dict with: gene, mutation, uniprot_id, chrom, strand, mut_type,
               cds_pos, genomic_positions, ref_matches_genome,
               window_start, window_end,
               ref_dna_coding, mut_dna_coding,
               ref_dna_genomic, mut_dna_genomic,
               warnings
    """
    mut = parse_dna_mutation(mutation_str)
    gene_upper = gene.strip().upper()
    warnings = []

    # ── gene → transcript ──
    gene_map = _get_gene_map()
    if gene_upper not in gene_map:
        raise ValueError(f"Gene '{gene}' not found in UniProt catalog")
    info = gene_map[gene_upper]
    uid, chrom, strand = info["uid"], info["chrom"], info["strand"]

    # ── CDS base map ──
    cds_map = _get_cds_map()
    if uid not in cds_map:
        raise ValueError(f"No CDS exons found for {uid} ({gene})")
    cds_bases = _get_cds_bases(uid, cds_map[uid], strand)
    genome = _get_genome()

    # ── resolve genomic positions for each mutation type ──

    mtype = mut["type"]
    pos, end_pos = mut["pos"], mut["end_pos"]
    offset, end_offset = mut["offset"], mut.get("end_offset")

    if mtype == "sub":
        # single genomic position (possibly intronic)
        gpos = _cds_to_genomic(pos, cds_bases, strand, offset)
        affected_gpos = [gpos]

        # validate ref base
        coding_base = _fetch_coding_bases([gpos], chrom, strand, genome)
        ref_ok = (mut["ref_bases"] == coding_base)
        if not ref_ok:
            warnings.append(
                f"Ref base mismatch at c.{pos}: expected '{mut['ref_bases']}', "
                f"genome has '{coding_base}'"
            )

    elif mtype == "del":
        if end_pos is not None:
            if offset or end_offset:
                # intronic range del — resolve both endpoints
                g_start = _cds_to_genomic(pos, cds_bases, strand, offset)
                g_end = _cds_to_genomic(end_pos, cds_bases, strand, end_offset)
                lo, hi = sorted([g_start, g_end])
                affected_gpos = list(range(lo, hi + 1))
            else:
                # exonic range — use CDS base map directly
                affected_gpos = [cds_bases[i] for i in range(pos - 1, end_pos)]
        else:
            gpos = _cds_to_genomic(pos, cds_bases, strand, offset)
            affected_gpos = [gpos]

        # validate ref bases if specified
        coding_bases = _fetch_coding_bases(
            [cds_bases[i] for i in range(pos - 1, (end_pos or pos))] if not offset else affected_gpos,
            chrom, strand, genome,
        )
        ref_ok = True
        if mut["ref_bases"] and coding_bases != mut["ref_bases"]:
            ref_ok = False
            warnings.append(
                f"Ref bases mismatch at c.{pos}: expected '{mut['ref_bases']}', "
                f"genome has '{coding_bases}'"
            )

    elif mtype == "dup":
        if end_pos is not None:
            affected_gpos = [cds_bases[i] for i in range(pos - 1, end_pos)]
        else:
            gpos = _cds_to_genomic(pos, cds_bases, strand, offset)
            affected_gpos = [gpos]

        coding_bases = _fetch_coding_bases(affected_gpos, chrom, strand, genome)
        ref_ok = True
        if mut["ref_bases"] and coding_bases != mut["ref_bases"]:
            ref_ok = False
            warnings.append(
                f"Ref bases mismatch at c.{pos}: expected '{mut['ref_bases']}', "
                f"genome has '{coding_bases}'"
            )

    elif mtype == "ins":
        g_start = _cds_to_genomic(pos, cds_bases, strand, offset)
        g_end = _cds_to_genomic(end_pos, cds_bases, strand, end_offset)
        affected_gpos = sorted([g_start, g_end])
        ref_ok = True  # nothing to validate for pure insertions

    elif mtype == "delins":
        if end_pos is not None:
            affected_gpos = [cds_bases[i] for i in range(pos - 1, end_pos)]
        else:
            affected_gpos = [cds_bases[pos - 1]]
        ref_ok = True

    else:
        raise ValueError(f"Unknown mutation type: {mtype}")

    # ── build genomic window ──
    affected_sorted = sorted(affected_gpos)
    center = affected_sorted[len(affected_sorted) // 2]
    half = window // 2
    win_start = max(0, center - half)
    win_end = center + (window - half)

    ref_window = genome.fetch(chrom, win_start, win_end).upper()
    mut_window = list(ref_window)

    # ── apply mutation to + strand window ──

    if mtype == "sub":
        idx = affected_gpos[0] - win_start
        mut_window[idx] = _complement(mut["alt_bases"]) if strand == "-" else mut["alt_bases"]

    elif mtype == "del":
        for gp in sorted(affected_sorted, reverse=True):
            idx = gp - win_start
            if 0 <= idx < len(mut_window):
                del mut_window[idx]

    elif mtype == "dup":
        # fetch + strand bases at affected positions (ascending genomic order)
        dup_seq = genome.fetch(chrom, affected_sorted[0],
                               affected_sorted[-1] + 1).upper()
        # insert after the highest genomic pos for + strand,
        # before the lowest for - strand
        if strand == "+":
            ins_idx = affected_sorted[-1] - win_start + 1
        else:
            ins_idx = affected_sorted[0] - win_start
        for i, base in enumerate(dup_seq):
            mut_window.insert(ins_idx + i, base)

    elif mtype == "ins":
        lo = min(affected_gpos)
        ins_idx = lo - win_start + 1
        ins_seq = _revcomp(mut["alt_bases"]) if strand == "-" else mut["alt_bases"]
        for i, base in enumerate(ins_seq):
            mut_window.insert(ins_idx + i, base)

    elif mtype == "delins":
        # delete affected bases
        for gp in sorted(affected_sorted, reverse=True):
            idx = gp - win_start
            if 0 <= idx < len(mut_window):
                del mut_window[idx]
        # insert alt bases where the first deleted base was
        lo = affected_sorted[0]
        ins_idx = lo - win_start
        ins_seq = _revcomp(mut["alt_bases"]) if strand == "-" else mut["alt_bases"]
        for i, base in enumerate(ins_seq):
            mut_window.insert(ins_idx + i, base)

    mut_window = "".join(mut_window)

    # ── coding strand (5'→3') ──
    if strand == "-":
        ref_coding = _revcomp(ref_window)
        mut_coding = _revcomp(mut_window)
    else:
        ref_coding = ref_window
        mut_coding = mut_window

    return {
        "gene": gene_upper,
        "mutation": mutation_str.strip(),
        "uniprot_id": uid,
        "chrom": chrom,
        "strand": strand,
        "mut_type": mtype,
        "cds_pos": pos,
        "cds_end_pos": end_pos,
        "intronic_offset": offset,
        "genomic_positions": affected_sorted,
        "ref_matches_genome": ref_ok,
        "window_start": win_start,
        "window_end": win_end,
        "ref_dna_coding": ref_coding,
        "mut_dna_coding": mut_coding,
        "ref_dna_genomic": ref_window,
        "mut_dna_genomic": mut_window,
        "warnings": warnings,
    }


# ── CLI formatting ─────────────────────────────────────────────────────────────

def format_result(r):
    gp = r["genomic_positions"]
    if len(gp) == 1:
        loc = f"{r['chrom']}:{gp[0] + 1}"
    else:
        loc = f"{r['chrom']}:{gp[0] + 1}-{gp[-1] + 1}"

    lines = [
        f"Gene:     {r['gene']} ({r['uniprot_id']})",
        f"Mutation: {r['mutation']}  [{r['mut_type']}]",
        f"Location: {loc} ({r['strand']} strand)",
        f"Ref check: {'OK' if r['ref_matches_genome'] else 'MISMATCH'}",
        "",
    ]

    for w in r["warnings"]:
        lines.append(f"  *** {w}")
    if r["warnings"]:
        lines.append("")

    lines.append(f"Coding strand (5'->3'):")
    lines.append(f"  ref ({len(r['ref_dna_coding'])} bp): {r['ref_dna_coding']}")
    lines.append(f"  mut ({len(r['mut_dna_coding'])} bp): {r['mut_dna_coding']}")
    lines.append("")
    lines.append(f"Genomic (+ strand), {r['chrom']}:{r['window_start'] + 1}-{r['window_end']}:")
    lines.append(f"  ref: {r['ref_dna_genomic']}")
    lines.append(f"  mut: {r['mut_dna_genomic']}")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description="DNA window around a coding-DNA mutation")
    ap.add_argument("gene", help="Gene symbol, e.g. TP53")
    ap.add_argument("mutation", help='e.g. "c.1799T>A", "c.241_242delCT", "c.5266dupC"')
    ap.add_argument("--window", type=int, default=100, help="window size in bp (default 100)")
    args = ap.parse_args()

    try:
        result = get_mutation_dna_window(args.gene, args.mutation, window=args.window)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(format_result(result))


if __name__ == "__main__":
    main()
