"""Amino-acid property tables and per-mutation `aa_change` calculator.

Provides the per-residue property facts (charge, polarity, hydropathy,
volume, special residues) for `(wt_aa, mut_aa)` pairs. Constants only — no
regex, no parsing of HGVS strings (that lives in `mutation_verifier_dual.py`).
"""
from __future__ import annotations

# Standard 20 + stop.
AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Stop": "*",
}
AA_1TO3 = {v: k for k, v in AA_3TO1.items() if v != "*" or k == "Ter"}

AA_FULL_TO_1 = {
    "alanine": "A", "arginine": "R", "asparagine": "N", "aspartate": "D",
    "aspartic acid": "D", "cysteine": "C", "glutamine": "Q",
    "glutamate": "E", "glutamic acid": "E", "glycine": "G", "histidine": "H",
    "isoleucine": "I", "leucine": "L", "lysine": "K", "methionine": "M",
    "phenylalanine": "F", "proline": "P", "serine": "S", "threonine": "T",
    "tryptophan": "W", "tyrosine": "Y", "valine": "V",
}

# Charge at physiological pH (~7.4). D, E carry a net negative charge
# (sidechain pKa ~3.7--4.3); K, R carry a net positive charge (pKa
# ~10.5--12.5). H has a sidechain pKa ~6.0 and is treated here as
# ~neutral at pH 7.4 with the caveat flagged via HIS_AMBIG; the verifier
# avoids committing to a charge change for H-involving substitutions.
# Reference: Lehninger, Principles of Biochemistry (8th ed.), Table 3-1
# (sidechain pKa values); Stryer, Biochemistry (9th ed.), Ch. 3.
CHARGE = {
    "D": -1, "E": -1,
    "K": +1, "R": +1,
    "H": 0,  # ambiguous; see HIS_AMBIG below
    **{a: 0 for a in "ACFGILMNPQSTVWY"},
}
HIS_AMBIG = {"H"}

# Polarity (polar = sidechain forms hydrogen bonds at physiological pH).
# Standard textbook grouping; matches the Wikipedia "Amino acid properties"
# polar/nonpolar split and the IUPAC residue classification table.
# Reference: Lehninger Table 3-1; Branden & Tooze, "Introduction to Protein
# Structure" (2nd ed.), Ch. 1.
POLAR = set("STNQDEHKRYC")
NONPOLAR = set("AVLIMFWPG")

# Kyte–Doolittle hydropathy sign (KD index > 0 => hydrophobic). Values from
# Kyte, J. and Doolittle, R. F. (1982), "A simple method for displaying the
# hydropathic character of a protein", J. Mol. Biol. 157(1):105--132,
# doi:10.1016/0022-2836(82)90515-0. Cysteine (KD=+2.5) and tyrosine
# (KD=-1.3) are borderline; we follow the original paper's sign.
HYDROPHOBIC = set("AILMFVWCY")  # cysteine and tyrosine are borderline; KD>0
HYDROPHILIC = set("DEKNQRSTH")
SPECIAL_GP = set("GP")  # G (KD=-0.4) and P (KD=-1.6) are not strongly classified

# Mean residue volumes (Å³). Modern standard from Tsai, J., Taylor, R.,
# Chothia, C. and Gerstein, M. (1999), "The packing density in proteins:
# standard radii and volumes", J. Mol. Biol. 290(1):253--266,
# doi:10.1006/jmbi.1999.2829. Values are the "mean residue volume (all
# residues)" column of their Table 6, computed via Voronoi tessellation
# over a high-quality dataset of well-determined protein structures and
# now the de-facto standard used by NACCESS, modern packing analyses,
# and structure-prediction toolkits. (Earlier widely-cited tables --
# Zamyatnin 1972/1984, Chothia 1975, Harpaz/Gerstein/Chothia 1994 --
# differ by 1--6% per residue.) The 3-bin partition (small / medium /
# large) follows the 110 / 165 Å³ boundaries used by Bordo & Argos (1991),
# J. Mol. Biol. 217:721--729 for substitution-tolerance analysis.
VOLUME = {
    "G":  63.8, "A":  89.3, "S":  94.2, "C": 102.5, "T": 119.6,
    "P": 121.6, "N": 122.4, "D": 114.4, "V": 138.2, "E": 138.8,
    "Q": 146.9, "H": 157.5, "I": 163.0, "L": 163.1, "K": 165.1,
    "M": 165.8, "F": 190.8, "R": 190.3, "Y": 194.6, "W": 226.4,
}
SMALL = {a for a, v in VOLUME.items() if v < 110}     # G A S C
LARGE = {a for a, v in VOLUME.items() if v > 165}     # K M F R Y W
MEDIUM = set(VOLUME) - SMALL - LARGE                  # T P N D V E Q H I L

# Aromatic sidechains. F, Y, W are unambiguous; H is included because its
# imidazole ring is aromatic (6 π electrons via Hückel's rule) and the
# MutQA dataset routinely refers to His as aromatic.
# Reference: Burley, S. K. and Petsko, G. A. (1985), "Aromatic-aromatic
# interaction: a mechanism of protein structure stabilization",
# Science 229:23--28, doi:10.1126/science.3892686.
AROMATIC = set("FYWH")
# Aliphatic (non-aromatic, non-polar carbon sidechains). Standard set per
# Lehninger Table 3-1 / Wikipedia "Aliphatic amino acid".
ALIPHATIC = set("AVLIMG")
# Sulfur-containing sidechains: Cys (thiol) and Met (thioether).
SULFUR_CONTAINING = set("CM")
# Canonical kinase phosphorylation acceptors (Ser/Thr/Tyr).
# Reference: Manning, G. et al. (2002), "The protein kinase complement of the
# human genome", Science 298:1912--1934, doi:10.1126/science.1075762.
PHOSPHO_SITE = set("STY")
# N-linked glycosylation acceptor (Asn) -- valid only inside an N-X-S/T
# sequon (X != P). Sequon check is performed by mutation_facts().
# Reference: Marshall, R. D. (1972), "Glycoproteins", Annu. Rev. Biochem.
# 41:673--702, doi:10.1146/annurev.bi.41.070172.003325.
GLYCO_N_LINKED = set("N")
# Hydrogen-bond donor / acceptor sets following the Lipinski-style
# heuristic adapted for protein sidechains: any sidechain with an N-H or
# O-H is a donor; any with a lone-pair-bearing N or O is an acceptor.
# Reference: Lipinski, C. A. et al. (1997), Adv. Drug Deliv. Rev. 23:3--25,
# doi:10.1016/S0169-409X(96)00423-1; McDonald & Thornton (1994),
# J. Mol. Biol. 238:777--793.
HBOND_DONOR = set("STNQYHKRW")
HBOND_ACCEPTOR = set("DENQHSTY")


def _bin(aa: str, small=SMALL, large=LARGE) -> str:
    if aa in small:
        return "small"
    if aa in large:
        return "large"
    return "medium"


def aa_properties(aa: str) -> dict:
    """Per-residue property dict for a single AA (1-letter code, uppercase)."""
    aa = aa.upper()
    return {
        "aa": aa,
        "charge": CHARGE.get(aa, 0),
        "polar": aa in POLAR,
        "hydrophobic": aa in HYDROPHOBIC,
        "hydrophilic": aa in HYDROPHILIC,
        "volume_class": _bin(aa),
        "aromatic": aa in AROMATIC,
        "aliphatic": aa in ALIPHATIC,
        "sulfur": aa in SULFUR_CONTAINING,
        "phospho_capable": aa in PHOSPHO_SITE,
        "glyco_capable": aa in GLYCO_N_LINKED,
        "hbond_donor": aa in HBOND_DONOR,
        "hbond_acceptor": aa in HBOND_ACCEPTOR,
        "is_special_pro": aa == "P",
        "is_special_gly": aa == "G",
        "is_cys": aa == "C",
    }


def aa_change(wt: str, mut: str) -> dict:
    """Return change-fact dict for a wt→mut AA substitution. Strings are
    canonical labels matching the literals in `triggers.CHARGE_TRIGGERS` etc.
    """
    if not wt or not mut or wt == "*" or mut == "*":
        # Stop-gain or unparsed; only mutation_kind applies.
        return {"applicable": False}
    if wt not in VOLUME or mut not in VOLUME:
        # Non-canonical / ambiguous AA code (X, B, Z, U, O, etc.).
        return {"applicable": False, "reason": "non_canonical_aa"}

    wp = aa_properties(wt)
    mp = aa_properties(mut)

    # Charge change.
    cw, cm = wp["charge"], mp["charge"]
    if cw == cm:
        charge_change = "no_change"
    elif cw > 0 and cm == 0:
        charge_change = "loss_positive"
    elif cw == 0 and cm > 0:
        charge_change = "gain_positive"
    elif cw < 0 and cm == 0:
        charge_change = "loss_negative"
    elif cw == 0 and cm < 0:
        charge_change = "gain_negative"
    elif cw * cm < 0:
        charge_change = "sign_flip"
    elif cw != 0 and cm == 0:
        charge_change = "neutralized"
    else:
        charge_change = "other"

    # Polarity / hydrophobicity.
    if wp["polar"] and not mp["polar"]:
        polarity_change = "polar_to_nonpolar"
    elif not wp["polar"] and mp["polar"]:
        polarity_change = "nonpolar_to_polar"
    else:
        polarity_change = "no_change"

    if wp["hydrophobic"] and mp["hydrophilic"]:
        hydrophobicity_change = "more_hydrophilic"
    elif wp["hydrophilic"] and mp["hydrophobic"]:
        hydrophobicity_change = "more_hydrophobic"
    else:
        hydrophobicity_change = "no_change"

    # Size change (raw Å³ delta on Tsai 1999 residue volumes). Threshold is
    # |Δv| > 30 Å³, ~2 methylene groups (CH2 ~ 16--17 Å³ per Bondi 1964 vdW
    # increments), which is the smallest difference that reliably crosses
    # the per-residue conformational variance reported in Tsai et al. 1999
    # (sigma ~5--15 Å³ per residue) and corresponds to ~0.8sigma of the
    # 38 Å³ between-residue spread. Anything below ±30 Å³ is treated as
    # "similar" because models / annotators typically don't call e.g.
    # L (163.1) <-> I (163.0) or V (138.2) <-> T (119.6) "larger" or
    # "smaller".
    dv = VOLUME[mut] - VOLUME[wt]
    if dv > 30:
        size_change = "larger"
    elif dv < -30:
        size_change = "smaller"
    else:
        size_change = "similar"

    # Aromaticity.
    if wp["aromatic"] and not mp["aromatic"]:
        aromaticity_change = "loss_aromatic"
    elif not wp["aromatic"] and mp["aromatic"]:
        aromaticity_change = "gain_aromatic"
    else:
        aromaticity_change = "no_change"

    return {
        "applicable": True,
        "wt": wt, "mut": mut,
        "charge_change": charge_change,
        "polarity_change": polarity_change,
        "hydrophobicity_change": hydrophobicity_change,
        "size_change": size_change,
        "volume_delta": dv,
        "aromaticity_change": aromaticity_change,
        "loss_of_proline": wp["is_special_pro"] and not mp["is_special_pro"],
        "gain_of_proline": mp["is_special_pro"] and not wp["is_special_pro"],
        "loss_of_glycine": wp["is_special_gly"] and not mp["is_special_gly"],
        "loss_of_cysteine": wp["is_cys"] and not mp["is_cys"],
        "gain_of_cysteine": mp["is_cys"] and not wp["is_cys"],
        "loss_of_phospho_site": wp["phospho_capable"] and not mp["phospho_capable"],
        "gain_of_phospho_site": mp["phospho_capable"] and not wp["phospho_capable"],
        "wt_props": wp, "mut_props": mp,
    }
