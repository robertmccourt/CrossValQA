"""DNA-coding mutation fact builder for the deterministic verifier.

Wraps `release/get_mutation_dna.py` (the same GRCh38 window extractor that is
shipped as a user-facing helper) to produce a facts dict shaped like
`mutation_facts()` in `mutation_verifier_dual.py`, so the existing
AA-chemistry checks (charge / polarity / hydrophobicity / size / special)
reuse for substitutions that translate to a clean predicted protein change.

For non-substitution kinds (del / ins / dup / delins) and intronic / splice
variants, only the DNA-side facts are populated — AA-chem checks are skipped
because the predicted protein change is either a frameshift (junk past the
indel) or doesn't apply.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
# Reference data (GRCh38 FASTA + UniProt-hg38 CDS tables) for the auditor.
# get_mutation_dna reads $CROSSVALQA_REF_DIR at import time; default it to
# <repo>/reference_data/ unless the caller already set it.
os.environ.setdefault("CROSSVALQA_REF_DIR", str(REPO_ROOT / "reference_data"))

sys.path.insert(0, str(REPO_ROOT / "release"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from get_mutation_dna import (  # noqa: E402
    get_mutation_dna_window,
    parse_dna_mutation,
    _get_genome,
    _get_cds_map,
    _get_cds_bases,
    _complement,
)
from aa_properties import aa_change  # noqa: E402


# Standard codon table (NCBI translation table 1).
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _get_codon(cds_bases, strand, cds_pos, chrom, genome):
    """Return (codon_str, offset_in_codon, codon_start_cds_1based) for a
    1-based exonic CDS position. None if out of CDS bounds.

    Uses cds_bases (the per-transcript map from coding-order index → genomic
    position), so the codon is always assembled from the right exons even
    when it spans an intron — the genome-window approach in
    get_mutation_dna_window cannot do this.
    """
    if cds_pos < 1 or cds_pos > len(cds_bases):
        return None
    offset_in_codon = (cds_pos - 1) % 3
    codon_start_idx = (cds_pos - 1) - offset_in_codon
    if codon_start_idx + 3 > len(cds_bases):
        return None
    gpositions = cds_bases[codon_start_idx:codon_start_idx + 3]
    plus = "".join(genome.fetch(chrom, p, p + 1) for p in gpositions).upper()
    if strand == "-":
        codon = _complement(plus)
    else:
        codon = plus
    return codon, offset_in_codon, codon_start_idx + 1


def _apply_sub_to_codon(codon, offset_in_codon, alt_base):
    return codon[:offset_in_codon] + alt_base.upper() + codon[offset_in_codon + 1:]


def mutation_facts_dna(gene, mutation_str, ref_seq=None):
    """Build a verifier-ready facts dict for a DNA-coding mutation.

    Returns None on parse failure or gene-resolution failure. On success,
    returns a dict that is a strict superset of the protein-side facts dict
    so it can be passed to the existing per-family check functions:

      Required-by-existing-checks fields:
        wt_aa, mut_aa, position, change, kind, protein_length, sanity_ok,
        wt_in_glyco_sequon (False for DNA), region_class (None for DNA-only)

      DNA-specific fields:
        notation="dna", cds_position, intronic_offset, ref_base, alt_base,
        ref_matches_genome, splice_region, region_class_dna,
        wt_codon, mut_codon, synonymous,
        uniprot_id, chrom, strand, raw_mutation, gene_used, warnings
    """
    try:
        parsed = parse_dna_mutation(mutation_str)
    except ValueError:
        return None
    try:
        result = get_mutation_dna_window(gene, mutation_str, window=60)
    except (ValueError, KeyError, IndexError):
        return None

    intronic_offset = result["intronic_offset"]
    region_class_dna = "cds"
    splice_region = False
    if intronic_offset is not None:
        if abs(intronic_offset) <= 2:
            region_class_dna = "splice_site"
            splice_region = True
        else:
            region_class_dna = "intron"

    facts = {
        # ── shape required by existing per-family check functions ──
        "wt_aa": "",
        "mut_aa": "",
        "position": None,
        "kind": result["mut_type"],
        "protein_length": len(ref_seq) if ref_seq else None,
        "sanity_ok": result["ref_matches_genome"],
        "ref_aa_at_pos": None,
        "wt_in_glyco_sequon": False,
        "region_class": None,  # AA-position region class — N/A for DNA-only
        "change": {"applicable": False},
        "mutation_type_label": result["mut_type"],
        # ── DNA-specific facts ──
        "notation": "dna",
        "cds_position": result["cds_pos"],
        "cds_end_position": result["cds_end_pos"],
        "intronic_offset": intronic_offset,
        "ref_base": parsed.get("ref_bases"),
        "alt_base": parsed.get("alt_bases"),
        "ref_matches_genome": result["ref_matches_genome"],
        "splice_region": splice_region,
        "region_class_dna": region_class_dna,
        "synonymous": False,
        "wt_codon": None,
        "mut_codon": None,
        "uniprot_id": result["uniprot_id"],
        "chrom": result["chrom"],
        "strand": result["strand"],
        "raw_mutation": mutation_str,
        "gene_used": result["gene"],
        "warnings": list(result["warnings"]),
    }

    # Predicted protein change — only for exonic substitutions with a known
    # alt base. For intronic substitutions (splice region), the AA change
    # depends on splicing outcome, not the codon, so we don't predict it.
    if (result["mut_type"] == "sub"
            and intronic_offset is None
            and facts["alt_base"]):
        cds_map = _get_cds_map()
        if result["uniprot_id"] in cds_map:
            cds_bases = _get_cds_bases(
                result["uniprot_id"],
                cds_map[result["uniprot_id"]],
                result["strand"],
            )
            genome = _get_genome()
            codon_info = _get_codon(
                cds_bases, result["strand"], result["cds_pos"],
                result["chrom"], genome,
            )
            if codon_info:
                codon, offset_in_codon, _ = codon_info
                mut_codon = _apply_sub_to_codon(
                    codon, offset_in_codon, facts["alt_base"]
                )
                wt_aa = CODON_TABLE.get(codon, "")
                mut_aa = CODON_TABLE.get(mut_codon, "")
                aa_pos = (result["cds_pos"] - 1) // 3 + 1
                facts["wt_aa"] = wt_aa
                facts["mut_aa"] = mut_aa
                facts["position"] = aa_pos
                facts["wt_codon"] = codon
                facts["mut_codon"] = mut_codon
                if wt_aa and mut_aa:
                    if wt_aa == mut_aa:
                        facts["synonymous"] = True
                        facts["kind"] = "synonymous"
                    elif mut_aa == "*":
                        facts["kind"] = "nonsense"
                    else:
                        facts["kind"] = "substitution"
                        facts["change"] = aa_change(wt_aa, mut_aa)

                # Region-class on the protein for AA-position reuse.
                if facts["protein_length"]:
                    L = facts["protein_length"]
                    third = max(1, L // 3)
                    if aa_pos <= third:
                        facts["region_class"] = "n_terminal"
                    elif aa_pos > L - third:
                        facts["region_class"] = "c_terminal"
                    else:
                        facts["region_class"] = "middle"

    return facts


if __name__ == "__main__":
    import json
    cases = [
        ("BRAF", "c.1799T>A"),       # V600E — classic missense, minus strand
        ("KRAS", "c.35G>A"),         # G12D
        ("CFTR", "c.1521_1523del"),  # F508del — in-frame deletion
        ("BRCA2", "c.316+1G>A"),     # canonical splice donor
        ("VPS33B", "c.89T>C"),       # L30P — bare-style mutation in v2 sample
    ]
    for gene, mut in cases:
        f = mutation_facts_dna(gene, mut)
        print(f"\n=== {gene} {mut} ===")
        if f is None:
            print("  PARSE/RESOLVE FAILED")
            continue
        print(json.dumps({k: v for k, v in f.items() if k != "warnings"},
                         indent=2, default=str))
        if f["warnings"]:
            print("  warnings:", f["warnings"])
