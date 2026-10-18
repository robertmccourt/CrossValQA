"""DNA-specific trigger patterns for the deterministic verifier.

Companion to `triggers.py` (which targets protein-level claims). This module
adds patterns that only make sense against DNA-coding mutation facts:

  - Inline c.* HGVS notation in answer text (`c.1799T>A`, bare `1799T>A`,
    `T1799A`, `c.5266dupC`, `c.241_242delCT`)
  - Inline p.* derivations from a DNA mutation (`p.Gly331Arg`, `p.G331R`,
    `(p.Val600Glu)`) — verified against the codon-translated prediction
  - Splice-region language (donor / acceptor / canonical splice site)
  - Mutation-kind aliases that overlap with DNA notation
    (`silent variant` / `synonymous` / `stop-gain` / `nonsense`)
  - Nucleotide-context claims (`G to A transition`, `transversion`)

Trigger families that already work for both protein and DNA — `region_class`
(N-/C-terminal), special-residue, and the AA-chemistry families (charge /
polarity / hydrophobicity / size) — live in `triggers.py` and are reused
unchanged when the DNA mutation is an exonic substitution that translates to
a known wt→mut AA pair.
"""
from __future__ import annotations

import re

# ─────────────────────────────────────────────────── inline DNA notation

# Match `c.1799T>A`, `c.5266dupC`, `c.241_242delCT`, `c.123_124insACG`, etc.
# Captures (cds_pos, [end_pos], [intronic_offset], kind_token, [bases])
DNA_HGVS_SUB = re.compile(
    r"\bc\.(\d+)([+-]\d+)?([ACGT])>([ACGT])\b",
    re.IGNORECASE,
)
# Bare nucleotide style without c. prefix: `1799T>A`
DNA_BARE_SUB = re.compile(
    r"(?<![A-Za-z0-9.])(\d+)([ACGT])>([ACGT])(?![A-Za-z0-9])",
)
# Pre-style: `T1799A`
DNA_PRE_SUB = re.compile(
    r"(?<![A-Za-z0-9])([ACGT])(\d+)([ACGT])(?![A-Za-z0-9])",
)
DNA_HGVS_DEL = re.compile(
    r"\bc\.(\d+)(?:_(\d+))?del([ACGT]*)\b",
    re.IGNORECASE,
)
DNA_HGVS_DUP = re.compile(
    r"\bc\.(\d+)(?:_(\d+))?dup([ACGT]*)\b",
    re.IGNORECASE,
)
DNA_HGVS_INS = re.compile(
    r"\bc\.(\d+)_(\d+)ins([ACGT]+)\b",
    re.IGNORECASE,
)
DNA_HGVS_DELINS = re.compile(
    r"\bc\.(\d+)(?:_(\d+))?delins([ACGT]+)\b",
    re.IGNORECASE,
)


# ─────────────────────────────────────────────────── inline protein derivation

# Match `p.Gly331Arg`, `p.G331R`, `Gly331Arg`, `G331R`, with optional p. prefix.
# Captures (wt_aa_token, position, mut_aa_token).
PROT_THREE = re.compile(
    r"(?:p\.)?\b(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|"
    r"Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter)(\d+)"
    r"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|"
    r"Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter|\*)\b",
    re.IGNORECASE,
)
PROT_ONE = re.compile(
    r"(?:p\.)?\b([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY*])\b"
)


# ─────────────────────────────────────────────────── splice-region language

# Splice-site claims — only fire as `verified` when intronic_offset ∈ {±1, ±2}
# AND the corresponding canonical splice phrasing appears.
SPLICE_DONOR_PATTERNS = [
    # +1 / +2 of intron from preceding exon = canonical 5' splice donor (GT)
    r"\bsplice\s+donor\b",
    r"\b5['′][- ]?splice\s+site\b",
    r"\bdonor\s+(?:splice\s+)?site\b",
    r"\b(?:GT|GU)\s+(?:donor|splice)\b",
    r"\bcanonical\s+(?:5['′][- ]?)?splice\s+donor\b",
]
SPLICE_ACCEPTOR_PATTERNS = [
    # -1 / -2 of intron preceding exon = canonical 3' splice acceptor (AG)
    r"\bsplice\s+acceptor\b",
    r"\b3['′][- ]?splice\s+site\b",
    r"\bacceptor\s+(?:splice\s+)?site\b",
    r"\bAG\s+(?:acceptor|splice)\b",
    r"\bcanonical\s+(?:3['′][- ]?)?splice\s+acceptor\b",
]
# Generic splice-related language — fires `verified` when ANY splice region
# (cds intronic offset != None), regardless of donor/acceptor.
SPLICE_GENERIC_PATTERNS = [
    r"\bsplice[- ]?(?:site|region|junction|variant|effect|alter\w*)\b",
    r"\bsplicing\s+(?:defect|alter\w*|abolish\w*|disrupt\w*)\b",
    r"\bcryptic\s+splice\b",
    r"\bexon\s+skipping\b",
    r"\bintron\s+retention\b",
    r"\baberrant\s+splicing\b",
    r"\bcanonical\s+splice\b",
]


SPLICE_DONOR_REGEX = re.compile(
    "|".join(f"(?:{p})" for p in SPLICE_DONOR_PATTERNS), re.IGNORECASE,
)
SPLICE_ACCEPTOR_REGEX = re.compile(
    "|".join(f"(?:{p})" for p in SPLICE_ACCEPTOR_PATTERNS), re.IGNORECASE,
)
SPLICE_GENERIC_REGEX = re.compile(
    "|".join(f"(?:{p})" for p in SPLICE_GENERIC_PATTERNS), re.IGNORECASE,
)


# ─────────────────────────────────────────────────── mutation-kind aliases

# These overlap with MUTTYPE_TRIGGERS in triggers.py but are listed here as a
# DNA-specific extension. The protein verifier already covers most of these;
# we add `silent` / `synonymous` since those only apply to DNA-level analysis
# of a coding substitution.
SYNONYMOUS_REGEX = re.compile(
    r"\b(?:synonymous|silent\s+(?:mutation|variant|substitution|change)|"
    r"silent\s+nucleotide\s+change|does\s+not\s+(?:change|alter)\s+the\s+amino\s+acid|"
    r"no\s+amino\s+acid\s+change|amino\s+acid\s+(?:is\s+)?(?:un|not\s+)?changed)\b",
    re.IGNORECASE,
)
MISSENSE_REGEX = re.compile(
    r"\b(?:missense|amino\s+acid\s+(?:change|substitution|replacement)|"
    r"non[- ]?synonymous)\b",
    re.IGNORECASE,
)
# Strict "nonsense" kind classification — only fires when the LLM is
# committing to nonsense as the *mutation kind*, not just describing a
# downstream consequence.
NONSENSE_KIND_REGEX = re.compile(
    r"\b(?:nonsense\s+(?:mutation|variant|substitution)|stop[- ]?gain(?:ed)?|"
    r"nonsense\s+codon)\b",
    re.IGNORECASE,
)
# Premature-stop / truncation language — describes the protein-level
# *consequence*. Verified for nonsense substitutions AND frameshift indels;
# only contradicted when the kind is a clean substitution (which can't
# truncate via frameshift) or synonymous.
PREMATURE_STOP_REGEX = re.compile(
    r"\b(?:premature\s+(?:stop|termination)(?:\s+codon)?|"
    r"truncation|truncated\s+protein|"
    r"introduces?\s+(?:a\s+)?stop\s+codon|"
    r"creates?\s+(?:a\s+)?stop\s+codon)\b",
    re.IGNORECASE,
)
# Backward-compat alias retained so any external imports keep working; aimed
# at strict-kind matching only.
NONSENSE_REGEX = NONSENSE_KIND_REGEX
FRAMESHIFT_REGEX = re.compile(
    r"\b(?:frameshift|frame[- ]shift|reading[- ]frame\s+shift|"
    r"out[- ]of[- ]frame)\b",
    re.IGNORECASE,
)
INFRAME_DEL_REGEX = re.compile(
    r"\bin[- ]?frame\s+(?:deletion|del)\b",
    re.IGNORECASE,
)


# ─────────────────────────────────────────────────── nucleotide context

# `G to A transition` / `C>T transition` / `transversion`
TRANSITION_REGEX = re.compile(
    r"\b(?:[ACGT]\s*(?:to|>)\s*[ACGT]\s+)?transition\b",
    re.IGNORECASE,
)
TRANSVERSION_REGEX = re.compile(
    r"\b(?:[ACGT]\s*(?:to|>)\s*[ACGT]\s+)?transversion\b",
    re.IGNORECASE,
)
PURINES = set("AG")
PYRIMIDINES = set("CT")


def is_transition(ref, alt):
    """Transition = purine↔purine or pyrimidine↔pyrimidine."""
    if not (ref and alt):
        return None
    r, a = ref.upper(), alt.upper()
    return (r in PURINES) == (a in PURINES)


# ─────────────────────────────────────────────────── 3-letter AA name maps

AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}


def norm_aa_token(token):
    if not token:
        return ""
    t = token.strip()
    if t == "*" or t.lower() == "ter":
        return "*"
    if len(t) == 1:
        return t.upper()
    return AA_3TO1.get(t.title(), "")


def _summary():
    return (
        "dna_triggers.py: "
        f"dna_hgvs_patterns=6, prot_inline_patterns=2, "
        f"splice_donor={len(SPLICE_DONOR_PATTERNS)}, "
        f"splice_acceptor={len(SPLICE_ACCEPTOR_PATTERNS)}, "
        f"splice_generic={len(SPLICE_GENERIC_PATTERNS)}, "
        "kind_regexes=5, nucleotide_context=2"
    )


if __name__ == "__main__":
    print(_summary())
