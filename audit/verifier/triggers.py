"""Trigger / disqualifier / hedge regex tables for the MutQA
deterministic verifier. Patterns derived from a 10-way scan of the 19,577-row
gold answer set; each pattern below is grounded in one or more verbatim
phrasings from `discovery/slice_*.md` (pattern-discovery notes, not
distributed; "slice NN row MM" comments below refer to them).

This module contains ONLY constants (regex strings + dicts/lists) — the
verifier logic lives in mutation_verifier_dual.py / mutation_verifier_dna.py and imports from here.

# Out of scope:
# The following claim families recur in the gold answers (see section-9 of
# every slice) but are NOT deterministically verifiable from the (wt_aa,
# mut_aa, position, mutation_type) tuple alone, so the verifier returns
# `not_applicable` for them:
#   - NMD / nonsense-mediated mRNA decay assertions
#   - Splicing claims ("cryptic splice donor", "exon skipping", "out-of-frame
#     deletion of 22 nucleotides of exon 5", etc.)
#   - Kinetic / quantitative constants (Km, Vmax, kcat, IC50, n-fold change,
#     %-of-wild-type-activity, fold-decrease in catalytic efficiency)
#   - Voltage-clamp electrophysiology numbers (mV shifts in V½ of activation /
#     inactivation, ramp-current amplitude, recovery time constants)
#   - Sub-cellular localization / trafficking ("retained in Golgi", "fails to
#     reach the plasma membrane", "mislocalizes to ER")
#   - Secretion / expression-level changes ("reduced secretion", "absent in
#     extracellular matrix", "lower surface expression")
#   - Binding-partner / interactome claims ("loss of interaction with BACH1",
#     "fails to bind ADAM22", "disrupts homodimer formation")
#   - Aggregate / oligomer phenotypes ("forms larger cytoplasmic aggregates",
#     "fails to assemble into signaling platforms")
#   - Drug-resistance / pharmacology claims ("resistant to entrectinib",
#     "1.5-fold lower affinity for mexiletine")
#   - Disease-mechanism narratives ("muscle fiber hyperexcitability", "leads
#     to NDI", "causes the patient phenotype")
"""

import re

# ---------------------------------------------------------------- AA name maps
# slice 01 row 159 (BRCA1, p.Glu1250Lys): "glutamic acid" / "lysine"
# slice 01 row 6 (GPR68, p.Phe129_Asn278del): "Phe", "Asn"
AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}

AA_FULL_TO_1 = {
    "alanine": "A", "arginine": "R", "asparagine": "N", "aspartate": "D",
    "aspartic acid": "D", "cysteine": "C", "glutamine": "Q", "glutamate": "E",
    "glutamic acid": "E", "glycine": "G", "histidine": "H", "isoleucine": "I",
    "leucine": "L", "lysine": "K", "methionine": "M", "phenylalanine": "F",
    "proline": "P", "serine": "S", "threonine": "T", "tryptophan": "W",
    "tyrosine": "Y", "valine": "V",
}

AA_1TO3 = {v: k for k, v in AA_3TO1.items()}


# --------------------------------------------------------- Mutation-type triggers
# Each value is a list of regex pattern strings; alternates were collected
# from the section-6 ("Mutation-type language") and section-1 buckets.
MUTTYPE_TRIGGERS = {
    "substitution": [
        # slice 01 row 159 (BRCA1, p.Glu1250Lys): "changes an acidic amino
        # acid (glutamic acid) to a basic one (lysine)"
        # slice 05 row 7880 (PDHA1, p.Arg263Gln): "results in a glycine
        # substitution for arginine"
        r"\bsubstitut(?:ion|ed?|es|ing)\b",
        # slice 01 row 35 (NR2E3, p.Gly88Val): "introduces a valine
        # substitution at a glycine residue"
        r"\bmissense\b",
        # slice 01 row 38 (EDA, p.Met364Thr): "replaces a nonpolar methionine
        # with a polar threonine"
        r"\breplac(?:e|es|ed|ing|ement)\s+(?:a\s+|an\s+|the\s+)?\w+\s+(?:with|by|to|for)\b",
        # slice 03 row 5874 (PAH, p.Gly344Asp): "changes the glycine at amino
        # acid position 344 to aspartic acid"
        r"\bchang(?:e|es|ed|ing)\s+(?:a\s+|an\s+|the\s+)?\w+\s+to\s+\w+\b",
        # slice 01 row 159 (BRCA1, p.Glu1250Lys): "change from an acidic amino
        # acid (glutamic acid) to a basic amino acid (lysine)"
        r"\bchange\s+from\s+\w+\s+to\s+\w+\b",
    ],
    "deletion": [
        # slice 01 row 0 (GPR68, p.Phe129_Asn278del): "deletes four of the
        # seven transmembrane helices"
        # slice 01 row 328 (SLC4A1, p.Ala400_Ala408del): "deletion of nine
        # amino acids"
        r"\bdeletion\b",
        r"\bdelet(?:e|es|ed|ing)\b",
        # slice 05 row 8224 (FOXP3, p.Glu251del): "(E251 deletion)"
        r"\bdel\b",
        # slice 01 row 8 (MBOAT7, p.Glu253_Ala259del): in-frame deletion of 7
        # residues — phrasing not always present, but covered by deletion+\b
        r"\bin-?frame\s+deletion\b",
    ],
    "insertion": [
        # slice 01 row 320 (MUSK, p.Thr14fs): related "duplication" verbiage
        # slice mentions; insertion phrasings appear less frequently but real:
        r"\binsertion\b",
        r"\binsert(?:s|ed|ing)\b",
    ],
    "duplication": [
        # slice 01 row 320 (MUSK, c.40dup, p.Thr14fs): the c.40dup HGVS form;
        # phrasing examples:
        r"\bduplicat(?:ion|ed?|es|ing)\b",
        r"\bdup\b",
    ],
    "nonsense": [
        # slice 01 row 3 (OTOF, p.Tyr1497Ter): "introduces a premature stop
        # codon"
        # slice 01 row 13 (OTOF, p.Gln829Ter): "creates a premature stop codon"
        r"\bnonsense\b",
        r"\bpremature\s+stop\s+codon\b",
        r"\bpremature\s+stop\s+signal\b",
        r"\bpremature\s+termination(?:\s+codon)?\b",
        # slice 01 row 491 (RFX5, p.Arg239Ter): "creates a stop codon that
        # terminates the protein"
        r"\bstop\s+codon\b",
        # slice 01 row 11 (EPOR, p.Trp439Ter): "truncating the last 70 amino
        # acids at the C-terminal end"
        # slice 01 row 13: "leading to truncation of the protein"
        r"\btruncat(?:e|es|ed|ing|ion)\b",
        r"\btruncated\s+protein\b",
    ],
    "frameshift": [
        # slice 01 row 131 (MERTK, p.Cys738fs): "causes a frameshift, leading
        # to premature termination"
        # slice 01 row 416 (NPR2, p.Ile364fs): "causes a reading frame shift"
        r"\bframeshift\b",
        r"\bframe[- ]shift\b",
        r"\breading\s+frame\s+shift\b",
        # slice 01 row 170 (BRCA1, p.Arg71Gly): "out-of-frame deletion"
        r"\bout[- ]of[- ]frame\b",
    ],
    "synonymous": [
        # slice 01 row 317 (POMGNT1, p.Phe212=) / row 386 (NR2E3, p.Pro365=):
        # synonymous-but-splice-affecting variants. Often called silent.
        r"\bsynonymous\b",
        r"\bsilent\s+(?:mutation|variant|substitution)\b",
    ],
}


# ---------------------------------------------- Charge / polarity / hydrophobicity
# Each entry: (regex string, claim_type, expected_change_value).
# expected_change_value is one of the literal strings produced by an
# aa_properties.aa_change(wt, mut) helper:
#   "loss_positive" | "loss_negative" | "gain_positive" | "gain_negative" |
#   "sign_flip" | "neutralized" | "introduce_charge" |
#   "polar_to_nonpolar" | "nonpolar_to_polar" | "polar_unchanged" |
#   "more_hydrophobic" | "less_hydrophobic" | "hydrophobicity_change"
CHARGE_TRIGGERS = [
    # slice 01 row 69 (SCN4A, p.Arg1448Cys): "loss of the positive charge at
    # the mutated residue"
    (r"\bloss\s+of\s+(?:a\s+|the\s+)?positive(?:ly\s+charged)?\s+(?:residue|amino\s+acid|charge)\b",
     "charge_change", "loss_positive"),
    # slice 01 row 384 (VCL, p.Arg975Gln): "alters a positively charged patch"
    (r"\b(?:remov(?:e|es|ed|ing|al)|eliminat(?:e|es|ed|ing))\s+(?:a\s+|the\s+)?positive(?:ly\s+charged)?\s+(?:residue|charge)\b",
     "charge_change", "loss_positive"),
    # slice 01 row 343 (SLC4A1, p.Gly771Asp): "introduces a negative charge in
    # a transmembrane segment"
    (r"\bintroduc(?:e|es|ed|ing|tion\s+of)\s+(?:a\s+|an\s+|the\s+)?negative(?:ly\s+charged)?\s+(?:residue|amino\s+acid|charge)\b",
     "charge_change", "gain_negative"),
    # slice 01 row 421 (ATP2A2, p.Gly769Arg): "replaces a hydrophobic amino
    # acid with a basic one"
    (r"\bintroduc(?:e|es|ed|ing|tion\s+of)\s+(?:a\s+|an\s+|the\s+)?(?:positive(?:ly\s+charged)?|basic)\s+(?:residue|amino\s+acid|charge)\b",
     "charge_change", "gain_positive"),
    # slice 01 row 215 (BRCA1, p.Trp1837Arg): "introducing a charged residue"
    (r"\bintroduc(?:e|es|ed|ing|tion\s+of)\s+(?:a\s+|an\s+|the\s+)?charged\s+(?:residue|amino\s+acid)\b",
     "charge_change", "introduce_charge"),
    # slice 01 row 159 (BRCA1, p.Glu1250Lys): "changes an acidic amino acid
    # (glutamic acid) to a basic one (lysine)"
    (r"\b(?:changes?|change\s+from)\s+(?:an?\s+)?acidic\s+(?:residue|amino\s+acid)?[^.]{0,30}\bto\s+(?:an?\s+)?basic\b",
     "charge_change", "sign_flip"),
    (r"\b(?:changes?|change\s+from)\s+(?:an?\s+)?basic\s+(?:residue|amino\s+acid)?[^.]{0,30}\bto\s+(?:an?\s+)?acidic\b",
     "charge_change", "sign_flip"),
    # generic "charge reversal" phrasing (commonly used in literature; observed
    # in section-2 of multiple slices)
    (r"\bcharge\s+reversal\b", "charge_change", "sign_flip"),
    (r"\breverses?\s+the\s+charge\b", "charge_change", "sign_flip"),
    # slice 01 row 411 (AR, p.Arg780Trp): "non-conservative change with
    # alterations in polarity or electric charge"
    (r"\bneutraliz(?:e|es|ed|ing|ation)\b", "charge_change", "neutralized"),
]

POLARITY_TRIGGERS = [
    # slice 01 row 38 (EDA, p.Met364Thr): "replaces a nonpolar methionine with
    # a polar threonine"
    (r"\b(?:replac(?:e|es|ed|ing)|chang(?:e|es|ed|ing)|substitut(?:e|es|ed|ing))\s+(?:a\s+|an\s+|the\s+)?(?:non[- ]?polar|nonpolar|hydrophobic|apolar)\s+\w+\s+(?:with|by|to|for)\s+(?:a\s+|an\s+)?polar\b",
     "polarity_change", "nonpolar_to_polar"),
    # slice 05 row 8336 (TSPAN7, p.Pro172His): "substitutes a non-polar amino
    # acid with a polar one"
    (r"\b(?:non[- ]?polar|nonpolar|apolar)\s+(?:amino\s+acid|residue)\s+with\s+(?:a\s+|an\s+)?polar\b",
     "polarity_change", "nonpolar_to_polar"),
    # slice 03 row 4627 (LIPC, p.Ala195Thr): "introduces a polar amino acid in
    # a highly conserved nonpolar region"
    (r"\bintroduc(?:e|es|ed|ing)\s+(?:a\s+|an\s+)?polar\s+(?:amino\s+acid|residue)\b",
     "polarity_change", "nonpolar_to_polar"),
    # slice 01 row 17 (SLC40A1, p.Asn144Asp): "replaces a neutral amino acid
    # with a polar residue"
    (r"\b(?:replac(?:e|es|ed|ing)|chang(?:e|es|ed|ing)|substitut(?:e|es|ed|ing))\s+(?:a\s+|an\s+|the\s+)?polar\s+\w+\s+(?:with|by|to|for)\s+(?:a\s+|an\s+)?(?:non[- ]?polar|nonpolar|hydrophobic|apolar)\b",
     "polarity_change", "polar_to_nonpolar"),
    # slice 01 row 411 (AR, p.Arg780Trp): "alterations in polarity"
    (r"\b(?:alteration|change|alter\w*)\s+(?:in|of)\s+polarity\b",
     "polarity_change", "polar_unchanged"),
]

HYDROPHOB_TRIGGERS = [
    # slice 01 row 17 (SLC40A1, p.Asn144Asp): "reducing the hydrophobicity of
    # the transmembrane domain"
    (r"\b(?:reduc(?:e|es|ed|ing|tion\s+of)|decreas(?:e|es|ed|ing))\s+(?:the\s+)?hydrophobicity\b",
     "hydrophobicity_change", "less_hydrophobic"),
    # slice 01 row 437 (ATP2A1, p.Leu67Arg): "alters the hydrophobicity of the
    # transmembrane helix M1"
    (r"\balter(?:s|ed|ing)?\s+(?:the\s+)?hydrophobicity\b",
     "hydrophobicity_change", "hydrophobicity_change"),
    # slice 01 row 489 (FLNC, p.Tyr2563Cys): "creates a smaller, less
    # hydrophobic residue"
    (r"\bless\s+hydrophobic\b", "hydrophobicity_change", "less_hydrophobic"),
    (r"\bmore\s+hydrophobic\b", "hydrophobicity_change", "more_hydrophobic"),
    # slice 01 row 450 (DDC, p.Arg447Cys): "increases hydrophobic surface
    # exposure"
    (r"\bincreas(?:e|es|ed|ing)\s+hydrophobic(?:ity)?\b",
     "hydrophobicity_change", "more_hydrophobic"),
    # slice 01 row 38 (EDA, p.Met364Thr): "disrupting hydrophobic interactions"
    (r"\bdisrupt(?:s|ed|ing)?\s+hydrophobic\s+(?:interactions?|packing|core)\b",
     "hydrophobicity_change", "hydrophobicity_change"),
]


# ---------------------------------------------------------- Size / volume / steric
SIZE_TRIGGERS = [
    # slice 01 row 1240 (PKD2, p.Leu656Trp): "introduces a bulky tryptophan
    # residue"
    # slice 01 row 1374 (PKLR, p.Gly275Arg): "bulky positively-charged arginine"
    (r"\bbulk(?:y|ier|iness)\b", "size_change", "larger"),
    # slice 01 row 836 (COL2A1, p.Gly504Ser): "bulkier serine"
    (r"\b(?:bulk(?:ier|y)|larger|bigger)\s+(?:residue|side[- ]?chain|amino\s+acid)\b",
     "size_change", "larger"),
    # slice 01 row 1350 (AK1, p.Gly40Arg): "replacing a small residue with a
    # larger one"
    (r"\b(?:replac\w+|substitut\w+|chang\w+)\s+(?:a\s+|an\s+|the\s+)?small\w*\s+\w+\s+(?:with|by|to|for)\s+(?:a\s+|an\s+)?(?:larger|bigger|bulky|bulkier)\b",
     "size_change", "larger"),
    # slice 01 row 489 (FLNC, p.Tyr2563Cys): "creates a smaller, less
    # hydrophobic residue"
    (r"\bsmaller\s+(?:residue|side[- ]?chain|amino\s+acid)\b",
     "size_change", "smaller"),
    # slice 01 row 1322 (ALB, p.Arg242His): "smaller side chain of the
    # substituted amino acid"
    (r"\bsmaller\s+side[- ]?chain\b", "size_change", "smaller"),
    # slice 01 row 8 (MBOAT7, p.Glu253_Ala259del): "steric hindrance to the
    # formation of an α-helix"
    (r"\bsteric\s+(?:clash|hindrance|occlusion)\b", "size_change", "larger"),
    # slice 01 row 190 (BRCA1, p.Met1775Arg): "creating a steric clash"
    (r"\bsteric\s+clash\b", "size_change", "larger"),
    # slice 01 row 1374 (PKLR, p.Gly275Arg): "smallest amino acid, glycine,
    # with a bulky positively-charged arginine"
    (r"\bsmallest\s+amino\s+acid\b", "size_change", "smaller"),
    # slice 01 row 1082 (PDZD7, p.Gly103Arg): "small polar glycine with a
    # large bulky basic arginine"
    (r"\blarge\s+bulky\b", "size_change", "larger"),
]


# ---------------------------------------------- Special-residue presence / loss
# Triggered only if the WT or mutant AA matches the relevant group.
SPECIAL_TRIGGERS = {
    "disulfide": {
        # slice 01 row 43 (EFEMP2, p.Arg279Cys): "near a conserved disulfide
        # bond" / "proper disulfide bond formation"
        # slice 01 row 1421 (GDF1, p.Leu248fs): "cystine-knot" / "disulfide
        # bond formation"
        # slice 01 row 1672 (DLL4, p.Cys390Arg): "critical cysteine residue
        # involved in forming a disulfide bond"
        "patterns": [
            r"\bdisulfide\s+(?:bond|bridge|bonds?)\b",
            r"\bdisulphide\s+(?:bond|bridge)\b",
            r"\bcystine[- ]?knot\b",
            r"\bS-?S\s+(?:bond|bridge)\b",
        ],
        "valid_when_aa_in": {"C"},
    },
    "proline_rigidity": {
        # slice 01 row 293 (CTSC, p.Leu196Pro): "introduction of a proline
        # residue, which disrupts the usual organization of the polypeptide
        # backbone"
        # slice 01 row 487 (FLNC, p.Pro2298Leu): "rigid proline to a leucine"
        # slice 03 row 7514 (SET, p.Trp82Gly): "reduced rigidity"
        "patterns": [
            r"\brigid\s+proline\b",
            r"\bproline\s+(?:residue|kink|backbone|rigid\w*|cis-?trans)\b",
            r"\bdisrupts?\s+(?:the\s+)?(?:usual\s+)?(?:organization\s+of\s+the\s+)?(?:polypeptide\s+)?backbone\b",
            r"\bbackbone\s+(?:rigid\w*|kink)\b",
        ],
        "valid_when_aa_in": {"P"},
    },
    "glycine_flexibility": {
        # slice 01 row 836 (COL2A1, p.Gly504Ser): "replaces a glycine residue
        # with a bulkier serine in the triple-helical domain"
        # slice 01 row 1082 (PDZD7, p.Gly103Arg): "replacing a small polar
        # glycine with a large bulky basic arginine"
        # slice 01 row 727 (APOA1, p.Lys131del): "more flexible structure"
        "patterns": [
            r"\bglycine\s+(?:flexibility|residue|backbone)\b",
            r"\bflexibility\s+of\s+(?:the\s+)?glycine\b",
            r"\b(?:loss|absence)\s+of\s+(?:backbone\s+)?flexibility\b",
        ],
        "valid_when_aa_in": {"G"},
    },
    "aromatic": {
        # slice 01 row 206 (BRCA1, p.Trp1718Cys): "invariant tryptophan in the
        # BRCT core"
        # slice 01 row 489 (FLNC, p.Tyr2563Cys): replaces tyrosine
        # slice 01 row 1870 (NTRK1, p.Arg761Trp): "bulky tryptophan residue"
        "patterns": [
            r"\baromatic\s+(?:ring|side[- ]?chain|stack(?:ing)?|residue|amino\s+acid)\b",
            r"\bπ-?stacking\b",
            r"\bpi-?stacking\b",
            r"\baromatic\s+interaction\b",
        ],
        "valid_when_aa_in": {"F", "Y", "W", "H"},
    },
    "phospho_site": {
        # slice 01 row 591 (ELP1, p.Arg696Pro): "potential threonine
        # phosphorylation site"
        # slice 01 row 595 (ELP1, p.Arg696Pro): "consensus serine/threonine
        # kinase phosphorylation site"
        # slice 01 row 1361 (SLC6A3, p.Ala559Val): "phosphorylation at
        # specific N-terminal serine residues (Ser-7 and Ser-13)"
        # slice 01 row 660 (EYS, p.Gly2945Glu): "loss of a phosphorylation
        # site"
        "patterns": [
            r"\bphosphoryl(?:at(?:ion|ed?|es))?\s+site\b",
            r"\bphosphorylation\s+(?:site|residue)\b",
            r"\b(?:serine|threonine|tyrosine)\s+phosphorylation\b",
            r"\bautophosphoryl\w+\b",
        ],
        "valid_when_aa_in": {"S", "T", "Y"},
    },
    "glycosylation": {
        # slice 01 row 278 (GYPC, p.Asn8Ser): "prevents the addition of an
        # N-linked glycan to the protein by changing the asparagine residue at
        # position 8 to serine, which cannot be glycosylated"
        "patterns": [
            r"\bN-?linked\s+(?:glycan|glycosylation)\b",
            r"\bN-?glycosylation\b",
            r"\bO-?linked\s+(?:glycan|glycosylation)\b",
            r"\bglycosylation\s+site\b",
            r"\bglycosylated?\b",
            r"\bglycan\s+(?:attachment|addition)\b",
        ],
        "valid_when_aa_in": {"N", "S", "T"},
    },
    "histidine_pH": {
        # slice 01 row 0 (GPR68, p.Phe129_Asn278del): "histidine residues
        # crucial for pH sensitivity"
        "patterns": [
            r"\bhistidine\s+(?:residues?|side[- ]?chain)\b",
            r"\bpH\s+sens(?:itivity|ing)\b",
        ],
        "valid_when_aa_in": {"H"},
    },
}


# ---------------------------------------------------------- Region / position
REGION_TRIGGERS = {
    "n_terminal": [
        # slice 01 row 294 (GALC, p.Ile250Thr): "50-54 kDa N-terminal fragment"
        # slice 01 row 1361 (SLC6A3, p.Ala559Val): "N-terminal serine residues"
        r"\bN[- ]?terminal\b",
        r"\bN[- ]?terminus\b",
        r"\bamino[- ]?terminal\b",
        r"\bamino[- ]?terminus\b",
    ],
    "c_terminal": [
        # slice 01 row 11 (EPOR, p.Trp439Ter): "C-terminal end of the protein"
        # slice 01 row 249 (EXPH5, p.Pro1929fs): "truncating the protein by 52
        # residues at the C-terminus"
        # slice 01 row 301 (CYP11A1, p.Ile279fs): "carboxy-terminal region"
        r"\bC[- ]?terminal\b",
        r"\bC[- ]?terminus\b",
        r"\bcarboxy[- ]?terminal\b",
        r"\bcarboxy[- ]?terminus\b",
    ],
    "active_site": [
        # slice 01 row 39 (SGPL1, p.Arg222Gln): "distorting the active site"
        # slice 01 row 1374 (PKLR, p.Gly275Arg): "near the putative active site"
        # slice 03 row 4426 (RGS9, p.Trp299Arg): "catalytic domain"
        r"\bactive\s+site\b",
        r"\bcatalytic\s+(?:site|residue|residues|center|centre|core|domain)\b",
        r"\benzymatic\s+site\b",
    ],
    "binding_site": [
        # slice 01 row 108 (BRCA1, p.Met1775Arg): "hydrophobic pocket
        # involved in phosphopeptide binding"
        # slice 01 row 306 (SUOX, p.Lys379Arg): "ligand-binding site"
        # slice 01 row 857 (FKTN, p.Gly300Arg): "active site cleft"
        # slice 01 row 1056 (KCNJ11, p.Arg192His): "ATP-binding pocket"
        r"\bbinding\s+(?:site|pocket|interface|cleft|domain|motif|region)\b",
        r"\bligand[- ]?binding\b",
        r"\bsubstrate[- ]?binding\b",
        r"\bATP[- ]?binding\b",
        r"\bDNA[- ]?binding\b",
        r"\bphosphopeptide\s+binding\b",
    ],
    "transmembrane": [
        # slice 01 row 0 (GPR68, p.Phe129_Asn278del): "four of the seven
        # transmembrane helices"
        # slice 01 row 2 (SLC33A1, p.Ser113Arg): "second transmembrane domain"
        # slice 01 row 17 (SLC40A1, p.Asn144Asp): "transmembrane domain"
        # slice 01 row 343 (SLC4A1, p.Gly771Asp): "transmembrane segment"
        # slice 01 row 437 (ATP2A1, p.Leu67Arg): "transmembrane helix M1"
        r"\btransmembrane\s+(?:domain|region|helix|helices|segment|segments)\b",
        r"\bTM\d+\b",
        r"\bTM\s+\d+\b",
        r"\bmembrane[- ]?spanning\s+(?:domain|region|helix)\b",
    ],
    "extracellular": [
        # slice 01 row 226 (CLDN16, p.Leu151Phe): "first extracellular loop"
        # slice 03 row 3961 (ABCA4, p.Arg1443Cys): "extracellular domain"
        r"\bextracellular\s+(?:domain|region|loop|segment)\b",
    ],
    "intracellular": [
        # slice 03 row 4071 (MPL, p.Gln186Ter): "intracellular domains needed
        # for signal transduction"
        r"\bintracellular\s+(?:domain|region|loop|segment|tail)\b",
        r"\bcytoplasmic\s+(?:domain|region|loop|tail)\b",
    ],
    "signal_peptide": [
        # slice 03 row 4556 (GNRH1, p.Leu7fs): "truncated signal peptide
        # missing its hydrophobic core"
        r"\bsignal\s+(?:peptide|sequence)\b",
        r"\bleader\s+(?:peptide|sequence)\b",
    ],
    "kinase_domain": [
        # slice 01 row 320 (MUSK, p.Thr14fs): "cytoplasmic tyrosine kinase
        # domain"
        # slice 01 row 1870 (NTRK1, p.Arg761Trp): "kinase catalytic domain"
        r"\bkinase\s+(?:domain|catalytic\s+domain)\b",
        r"\btyrosine\s+kinase\s+domain\b",
        r"\bhelicase\s+domain\b",
    ],
    "nuclear_localization": [
        # slice 01 row 323 (PARN, p.Arg237Ter): "nuclear localization sequence"
        r"\bnuclear\s+localization\s+(?:signal|sequence)\b",
        r"\bNLS\b",
    ],
    "loop": [
        # slice 01 row 109 (SCN4A, p.Gly1306Glu): "III-IV cytoplasmic loop"
        # slice 01 row 487 (FLNC, p.Pro2298Leu): "loop between beta strands"
        r"\b(?:cytoplasmic|extracellular|intracellular|interdomain)\s+loop\b",
        r"\bloop\s+between\s+\w+\s+(?:strands?|helices)\b",
    ],
    "helix": [
        # slice 01 row 39 (SGPL1, p.Arg222Gln): "helices at the dimer
        # interface"
        # slice 01 row 209 (BRCA1, p.Ala1752Pro): "linker helix"
        # slice 01 row 548 (HFE, p.Leu183Pro): "α2-helix"
        # slice 01 row 616 (SIX3, p.Arg257Pro): "helix 3 of the homeodomain"
        r"\b(?:α|alpha|a)-?helix\b",
        r"\b(?:β|beta|b)-?(?:sheet|strand|hairpin|turn)\b",
        r"\blinker\s+helix\b",
        r"\bcoiled[- ]?coil\b",
    ],
}

# Numeric position claims like "residue 253", "position 253", "codon 253",
# "amino acid 253".
#   slice 01 row 281 (GYG1, p.Thr83Met): "glucose-O-tyrosine linkage at
#     position 195"
#   slice 01 row 591 (ELP1, p.Arg696Pro): "phosphorylation site at residue 699"
#   slice 01 row 131 (MERTK, p.Cys738fs): "premature termination before codon
#     770"
#   slice 01 row 82 (SCN4A, p.Gly1306Val): "mutations at position 1306"
POSITION_NUM_PATTERN = re.compile(
    r"\b(?:residue|position|codon|amino\s+acid)\s+(\d+)\b",
    re.IGNORECASE,
)

# Protein-length statements like "1,230 amino acid length" or "509 amino
# acids".
#   slice 01 row 3 (OTOF, p.Tyr1497Ter): "instead of its full 1,230 amino
#     acid length"
#   slice 01 row 588 (APOB, p.Arg2522Ter): "truncated protein that is 2494
#     amino acids long"
#   slice 03 row 3923 (SOX9, p.Ser99fs): "509 amino acids" / "108 amino acids"
PROTEIN_LENGTH_PATTERN = re.compile(
    r"\b(\d{1,4}(?:,\d{3})?)\s*[- ]?(?:amino\s+acid|residue|aa)s?\s+(?:protein|polypeptide|long|length|in\s+length)\b",
    re.IGNORECASE,
)


# --------------------------------------- HGVS-style direct claim notation
# slice 01 row 82 (SCN4A): "G1306E"
# slice 01 row 442 (RARA): "Arg276Gln" / "Arg276Ala"
# slice 01 row 600 (VAV1, p.Tyr841Phe): "Tyr-826" "Tyr-841"
HGVS_SHORT_PATTERN = re.compile(
    r"\b(?:p\.)?([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY*])\b"
)
HGVS_LONG_PATTERN = re.compile(
    r"\b(?:p\.)?(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter)"
    r"(\d+)"
    r"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter|\*)\b",
    re.IGNORECASE,
)


# ---------------------------------------- Hedge / hypothetical / negation
# slice 01 row 0 (GPR68, p.Phe129_Asn278del): "could be unstable"
# slice 01 row 8 (MBOAT7, p.Glu253_Ala259del): "potentially impairing" / "may
#   affect"
# slice 01 row 17 (SLC40A1, p.Asn144Asp): "likely affects" / "potentially
#   affecting"
# slice 01 row 56 (HNF4A, p.Val177Ile): "potentially modifying"
# slice 07 hedging samples: "likely", "would", "potentially", "appears",
#   "expected", "predicted"
HEDGE_REGEX = re.compile(
    r"\b(?:may|might|could|would|likely|possibly|presumably|potentially|"
    r"appears?|seems?|suggests?|suggesting|predict(?:s|ed|ing|ion)?|"
    r"expect(?:s|ed|ing|ation)?|hypothetic\w*|theoretic\w*|"
    r"approximately|about|roughly|around|nearly|almost|"
    r"~|≈)\b",
    re.IGNORECASE,
)

# slice 01 row 56 (HNF4A, p.Val177Ile): "may alter the protein's transcription
#   function"
# slice 01 row 57 (HNF4A, p.Val177Ile): "may cause a conformational change"
# slice 01 row 442 (RARA, p.Arg276Trp): "while the Arg276Gln mutation caused"
#   (compares to a different variant — hypothetical for THIS row)
# slice 01 row 30 (NTRK1, p.Gly667Cys): "remains sensitive to cabozantinib"
HYPOTHETICAL_REGEX = re.compile(
    r"\b(?:would|could|might|"
    r"if\s+(?:it|the\s+\w+|this\s+\w+)\s+(?:were|was|is|had|has)|"
    # comparative / analog discussion
    r"unlike|in\s+contrast\s+to|compared\s+to|versus|vs\.?|"
    r"similar\s+to|like\s+(?:the|an?)\s+\w+|by\s+analogy|"
    # other-variant discussion
    r"in\s+(?:other|similar)\s+(?:proteins?|variants?)|"
    r"in\s+(?:a\s+)?(?:paralog|ortholog|homolog)\w*|"
    # synthesis / modification framing
    r"replac(?:e|ing|ement)\s+(?:with|by)|"
    r"substitut(?:e|ion|ed|ing)\s+(?:with|by|for)|"
    r"upon\s+(?:substitution|mutation|replacement)|"
    r"if\s+used|if\s+the\s+mutation\s+were|"
    # other-mutation framing — see slice 01 row 442 (RARA, p.Arg276Trp)
    r"the\s+(?:Arg|Lys|Glu|Asp|Cys|His|Phe|Tyr|Trp|Met|Pro|Gly|Ala|Val|Leu|Ile|Ser|Thr|Asn|Gln)\d+(?:Arg|Lys|Glu|Asp|Cys|His|Phe|Tyr|Trp|Met|Pro|Gly|Ala|Val|Leu|Ile|Ser|Thr|Asn|Gln|\*)\b)",
    re.IGNORECASE,
)

# slice 01 row 32 (SLC40A1, p.Gly80Ser): "while maintaining normal
#   localization at the cell surface"
# slice 01 row 68 (WRN, p.Lys577Met): "does not affect its exonuclease
#   activity"
# slice 01 row 1696 (KCNJ10, p.Leu218Phe): "does not destroy"
# slice 07 (none observed) — but standard negation is well attested across
#   slices.
NEGATION_REGEX = re.compile(
    r"\b(?:no|not|without|"
    r"lack(?:s|ing|ed)?|absen(?:t|ce)(?:\s+of)?|devoid\s+of|free\s+of|"
    r"does\s+not\s+(?:have|contain|possess|affect|alter|change|disrupt)|"
    r"doesn'?t\s+(?:have|contain|possess|affect|alter|change|disrupt)|"
    r"do\s+not\s+(?:have|contain|possess|affect|alter|change|disrupt)|"
    r"unaffected|unchanged|preserved|retained|retains?|maintain(?:s|ed|ing)?|"
    r"missing|neither|"
    r"normal\s+(?:function|expression|localization|activity)|"
    r"fails?\s+to|failed\s+to)\b",
    re.IGNORECASE,
)


# --------------------------------------- Disqualifier per claim family
# Pre-window disqualifier: ~30-char window BEFORE the trigger that suppresses
# the claim (e.g. "in similar proteins", "homologous to", "the paralog X
# has", "literature suggests", "in contrast to", "compared to").
DISQUALIFIER_PRE = {
    # Residue identity: comparisons with other variants / paralogs / homologs
    # — see slice 01 row 442 (RARA, p.Arg276Trp): "While the Arg276Gln
    # mutation caused..." (the model is talking about a *different* variant).
    "residue_identity": re.compile(
        r"\b(?:paralog|orthol(?:og|ogue)|homolog(?:ue|ous)|"
        r"in\s+similar\s+proteins?|in\s+the\s+wild[- ]?type|by\s+analogy|"
        r"compared\s+to|in\s+contrast\s+to|unlike|"
        r"the\s+(?:Arg|Lys|Glu|Asp|Cys|His|Phe|Tyr|Trp|Met|Pro|Gly|Ala|Val|Leu|Ile|Ser|Thr|Asn|Gln)\d+\w{3}\s+mutation)\s*$",
        re.IGNORECASE,
    ),
    # Charge change: "compared to ...", "the wild-type protein has..."
    "charge_change": re.compile(
        r"\b(?:wild[- ]?type|wt|control|normal\s+protein|reference|"
        r"in\s+contrast\s+to|compared\s+to|unlike)\s*$",
        re.IGNORECASE,
    ),
    # Polarity: same comparison disqualifiers
    "polarity_change": re.compile(
        r"\b(?:wild[- ]?type|wt|control|normal\s+protein|"
        r"in\s+contrast\s+to|compared\s+to|unlike)\s*$",
        re.IGNORECASE,
    ),
    # Hydrophobicity: same
    "hydrophobicity_change": re.compile(
        r"\b(?:wild[- ]?type|wt|control|normal\s+protein|"
        r"in\s+contrast\s+to|compared\s+to|unlike)\s*$",
        re.IGNORECASE,
    ),
    # Size: ditto
    "size_change": re.compile(
        r"\b(?:wild[- ]?type|wt|control|normal\s+protein|"
        r"in\s+contrast\s+to|compared\s+to|unlike)\s*$",
        re.IGNORECASE,
    ),
    # Special residue: "such as cysteine residues in similar proteins" type
    # framings; also negation "loss of a..." should NOT disqualify here (it's
    # a real claim).
    "special_residue": re.compile(
        r"\b(?:in\s+similar\s+proteins?|in\s+other\s+proteins?|"
        r"by\s+analogy|in\s+the\s+homolog\w*)\s*$",
        re.IGNORECASE,
    ),
    # Region: "in the homologous protein..."
    "region": re.compile(
        r"\b(?:in\s+similar\s+proteins?|by\s+analogy|"
        r"in\s+the\s+homolog\w*|in\s+contrast\s+to)\s*$",
        re.IGNORECASE,
    ),
    # Mutation type: "unlike a frameshift, this mutation..."
    "mutation_type": re.compile(
        r"\b(?:unlike\s+(?:a|an)|in\s+contrast\s+to\s+(?:a|an)|"
        r"as\s+opposed\s+to\s+(?:a|an))\s*$",
        re.IGNORECASE,
    ),
}

# Post-window disqualifiers — context AFTER the trigger that re-frames it.
# slice 03 row 3934 (ARHGEF10, p.Thr332Ile): "similar to the effect observed
#   with an N-terminal truncated mutant" — the N-terminal claim is about a
#   different construct, not THIS mutation.
DISQUALIFIER_POST = {
    "residue_identity": re.compile(
        r"^\s*(?:in\s+(?:similar|other)\s+proteins?|"
        r"of\s+(?:the\s+)?(?:paralog|ortholog|homolog)\w*)",
        re.IGNORECASE,
    ),
    "region": re.compile(
        r"^\s*(?:mutant|construct|variant\s+of\s+a\s+different\s+protein)",
        re.IGNORECASE,
    ),
    "special_residue": re.compile(
        r"^\s*(?:in\s+(?:similar|other)\s+proteins?|of\s+the\s+homolog\w*)",
        re.IGNORECASE,
    ),
}


# --------------------------------------- Pre-compilation
# Compiled regex dicts/lists so the verifier can use them without recompiling.
MUTTYPE_TRIGGER_REGEX = {
    k: re.compile("|".join(f"(?:{p})" for p in v), re.IGNORECASE)
    for k, v in MUTTYPE_TRIGGERS.items()
}

CHARGE_TRIGGER_COMPILED = [
    (re.compile(p, re.IGNORECASE), claim_type, expected)
    for (p, claim_type, expected) in CHARGE_TRIGGERS
]
POLARITY_TRIGGER_COMPILED = [
    (re.compile(p, re.IGNORECASE), claim_type, expected)
    for (p, claim_type, expected) in POLARITY_TRIGGERS
]
HYDROPHOB_TRIGGER_COMPILED = [
    (re.compile(p, re.IGNORECASE), claim_type, expected)
    for (p, claim_type, expected) in HYDROPHOB_TRIGGERS
]
SIZE_TRIGGER_COMPILED = [
    (re.compile(p, re.IGNORECASE), claim_type, expected)
    for (p, claim_type, expected) in SIZE_TRIGGERS
]

SPECIAL_TRIGGER_COMPILED = {
    name: {
        "regex": re.compile(
            "|".join(f"(?:{p})" for p in spec["patterns"]), re.IGNORECASE,
        ),
        "valid_when_aa_in": set(spec["valid_when_aa_in"]),
    }
    for name, spec in SPECIAL_TRIGGERS.items()
}

REGION_TRIGGER_REGEX = {
    k: re.compile("|".join(f"(?:{p})" for p in v), re.IGNORECASE)
    for k, v in REGION_TRIGGERS.items()
}


def _summary() -> str:
    """Total counts per category — printed when run as __main__."""
    n_muttype = sum(len(v) for v in MUTTYPE_TRIGGERS.values())
    n_charge = len(CHARGE_TRIGGERS)
    n_polar = len(POLARITY_TRIGGERS)
    n_hydro = len(HYDROPHOB_TRIGGERS)
    n_size = len(SIZE_TRIGGERS)
    n_special = sum(len(v["patterns"]) for v in SPECIAL_TRIGGERS.values())
    n_region = sum(len(v) for v in REGION_TRIGGERS.values())
    n_dq_pre = len(DISQUALIFIER_PRE)
    n_dq_post = len(DISQUALIFIER_POST)
    return (
        f"triggers.py: muttype={n_muttype} ({len(MUTTYPE_TRIGGERS)} groups), "
        f"charge={n_charge}, polarity={n_polar}, hydrophob={n_hydro}, "
        f"size={n_size}, special={n_special} ({len(SPECIAL_TRIGGERS)} groups), "
        f"region={n_region} ({len(REGION_TRIGGERS)} groups), "
        f"disqualifier_pre={n_dq_pre}, disqualifier_post={n_dq_post}"
    )


if __name__ == "__main__":
    print(_summary())
