#!/usr/bin/env python3
"""
Step 2: Run Gene Recovery (Parallel)

Uses pre-extracted cited sentences from Step 1 to run gene recovery on all records.
No article lookup needed - uses the 'cited_sentences' field directly.

Input:
- m2tqa_with_cited.jsonl (473,243 records with cited_sentences field)

Output:
- m2tqa_gene_recovered.jsonl (same records with recovery fields added)

New fields added:
- recovered_gene: str or null
- recovery_tier: str
- recovery_confidence: "high", "medium", or "unresolved"
- recovery_score: int (0-350)
"""

import argparse
import json
import re
from collections import defaultdict
from typing import List, Dict, Optional, Set, Tuple
from dataclasses import dataclass, field
from multiprocessing import Pool, cpu_count
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

# ============================================================================
# CONFIGURATION
# ============================================================================

INPUT_PATH = str(config.CLAIM_EXTRACTOR_DIR / "phase_4" / "output" / "m2tqa_with_cited.jsonl")
OUTPUT_PATH = str(config.CLAIM_EXTRACTOR_DIR / "phase_4" / "output" / "m2tqa_gene_recovered.jsonl")
NUM_WORKERS = 64


# ============================================================================
# CONSTANTS - NON-GENE WORDS
# ============================================================================

NON_GENES = {
    # Common English words
    'THE', 'AND', 'FOR', 'NOT', 'WAS', 'ARE', 'HAS', 'HAD', 'BUT', 'THIS',
    'THAT', 'WITH', 'FROM', 'HAVE', 'BEEN', 'WERE', 'WILL', 'CAN', 'MAY',
    'HOW', 'WHAT', 'DOES', 'WHICH', 'WHERE', 'WHEN', 'WHO', 'WHY', 'ALL',
    'ONE', 'TWO', 'NEW', 'OLD', 'USE', 'USED', 'TYPE', 'BASED', 'KNOWN',
    'NOVEL', 'THESE', 'THOSE', 'SUCH', 'OTHER', 'EACH', 'BOTH', 'SAME',
    'INTO', 'ONLY', 'ALSO', 'MOST', 'SOME', 'MANY', 'MORE', 'LESS', 'THAN',
    'THEN', 'WHEN', 'WHILE', 'SINCE', 'AFTER', 'BEFORE', 'DURING', 'UNTIL',
    'GENE', 'GENES',
    # Prepositions
    'IN', 'ON', 'AT', 'TO', 'OF', 'BY', 'AS', 'OR', 'IF', 'SO',
    # Molecular biology terms
    'DNA', 'RNA', 'PCR', 'SNP', 'YES', 'NGS', 'WES', 'VAF', 'LOH', 'CNV',
    'HGVS', 'PMID', 'EXON', 'INTRON', 'CODON', 'ALLELE', 'ALLELES',
    # Lab techniques / imaging
    'IHC', 'FISH', 'CT', 'MRI', 'PET', 'ELISA', 'FACS', 'CRISPR',
    # Organizations
    'WHO', 'MSKCC', 'TCGA', 'COSMIC', 'ACMG', 'AMP', 'CAP', 'NCCN',
    # Figure/table references
    'FIG', 'FIGS', 'TABLE', 'TABLES', 'TAB', 'REF', 'REFS', 'SUPP', 'SUPPL',
    # Protein domains
    'DBD', 'HRDC', 'RQC', 'NLS', 'NES', 'TAD', 'SH2', 'SH3', 'PH', 'PDZ',
    # Amino acid codes (3-letter) and related
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'AMINO', 'ACID', 'ACIDS', 'SERINE', 'ALANINE', 'GLYCINE', 'VALINE',
    'LEUCINE', 'ISOLEUCINE', 'PROLINE', 'PHENYLALANINE', 'TYROSINE',
    'TRYPTOPHAN', 'METHIONINE', 'CYSTEINE', 'ASPARAGINE', 'GLUTAMINE',
    'ASPARTATE', 'GLUTAMATE', 'LYSINE', 'ARGININE', 'HISTIDINE', 'THREONINE',
    # Mutation terms
    'MUTATION', 'MUTATIONS', 'VARIANT', 'VARIANTS', 'DELETION', 'DELETIONS',
    'INSERTION', 'INSERTIONS', 'FRAMESHIFT', 'MISSENSE', 'NONSENSE', 'SPLICE',
    'TRUNCATING', 'TRUNCATED', 'PATHOGENIC', 'BENIGN', 'CHANGE', 'CHANGES',
    'GENETIC', 'GENOMIC', 'PROTEIN', 'SEQUENCE', 'DOMAIN', 'RESIDUE',
    'DISRUPTIVE', 'HETEROZYGOUS', 'HOMOZYGOUS', 'HEMIZYGOUS', 'COMPOUND',
    'NULL', 'LOSS', 'GAIN', 'FUNCTION', 'FUNCTIONAL', 'NONFUNCTIONAL',
    # Clinical descriptors
    'SEVERE', 'MILD', 'MODERATE', 'CLASSIC', 'CLASSICAL', 'ATYPICAL',
    'TYPICAL', 'EARLY', 'LATE', 'ONSET', 'ADULT', 'CHILDHOOD', 'INFANTILE',
    'CONGENITAL', 'ACQUIRED', 'PROGRESSIVE', 'STABLE', 'RECURRENT',
    # Ethnicity/nationality terms
    'RUSSIAN', 'CHINESE', 'JAPANESE', 'KOREAN', 'INDIAN', 'EUROPEAN',
    'AFRICAN', 'AMERICAN', 'ASIAN', 'CAUCASIAN', 'HISPANIC', 'TAIWANESE',
    'FOUNDER', 'ETHNIC', 'POPULATION',
    # Medical terms
    'PATIENT', 'PATIENTS', 'PROBAND', 'FAMILY', 'FAMILIES', 'CASE', 'CASES',
    'DISEASE', 'DISORDER', 'SYNDROME', 'TUMOR', 'TUMOUR', 'CANCER', 'NORMAL',
    'AFFECTED', 'UNAFFECTED', 'CARRIER', 'CARRIERS', 'HEALTHY', 'CONTROL',
    # Common false positives from output
    'CONTEXT', 'UNKNOWN', 'PRESENT', 'ABSENT', 'POSITIVE', 'NEGATIVE',
    'RESULTS', 'ANALYSIS', 'STUDY', 'STUDIES', 'REPORT', 'REPORTS',
    'NOVEL', 'KNOWN', 'IDENTIFIED', 'DETECTED', 'FOUND', 'DESCRIBED',
    'ASSOCIATED', 'RELATED', 'LINKED', 'CAUSED', 'CAUSING',
    # ACMG evidence codes (pathogenic)
    'PVS1', 'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    # ACMG evidence codes (benign)
    'BA1', 'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7',
    # Disease abbreviations (not genes)
    'PKU', 'CF', 'DMD', 'SMA', 'HD', 'ALS', 'MS', 'AD', 'PD',
    'MODY', 'LQTS', 'HCM', 'DCM', 'ARVC', 'CPVT', 'BrS',
    'MSD', 'MPS', 'NCL', 'CMT', 'HSAN', 'HSP', 'SCA', 'DRPLA',
    'FTD', 'PSP', 'CBD', 'MSA', 'DLB', 'VCP', 'ARC',
    # Cofactors and metabolites
    'BH4', 'FAD', 'NAD', 'NADH', 'NADP', 'NADPH', 'ATP', 'ADP', 'GTP', 'GDP',
    'SAM', 'THF', 'PLP', 'TPP', 'COA', 'ACOA',
    # Intron/exon notation
    'IVS1', 'IVS2', 'IVS3', 'IVS4', 'IVS5', 'IVS6', 'IVS7', 'IVS8', 'IVS9',
    'IVS10', 'IVS11', 'IVS12', 'IVS13', 'IVS14', 'IVS15', 'IVS16', 'IVS17',
    'IVS18', 'IVS19', 'IVS20', 'IVS21', 'IVS22', 'IVS23', 'IVS24', 'IVS25',
    # Family/pedigree IDs
    'OP1', 'OP2', 'OP3', 'OP4', 'OP5', 'OP6', 'OP7', 'OP8', 'OP9', 'OP10',
    # Misc false positives
    'ID', 'MSD1', 'PH1', 'C6B', 'C6A',
}


# ============================================================================
# PATTERNS
# ============================================================================

GENE_TOKEN = r'(?P<gene>(?=[A-Z0-9-]*[A-Z])[A-Z][A-Z0-9-]{1,14})'

GENE_PATTERNS = [
    re.compile(rf'\b{GENE_TOKEN}\s*[:(,]\s*(?:p\.|c\.)'),
    re.compile(rf'\b{GENE_TOKEN}\s+(?:p\.|c\.)'),
    re.compile(rf'\b{GENE_TOKEN}\s+[Ee][Xx][Oo][Nn]\s+\d+\s*[:(, ]?\s*(?:p\.|c\.)'),
    re.compile(rf'\b{GENE_TOKEN}\s+[Ii][Nn][Tt][Rr][Oo][Nn]\s+\d+\s*[:(, ]?\s*(?:p\.|c\.)'),
    re.compile(rf'(?:p\.|c\.)[^\s,;){{}}]{{1,30}}\s+[Ii][Nn]\s+(?:[Tt][Hh][Ee]\s+)?{GENE_TOKEN}\b'),
    re.compile(rf'(?:p\.|c\.)[^\s,;){{}}]{{1,30}}\s+[Oo][Ff]\s+(?:[Tt][Hh][Ee]\s+)?{GENE_TOKEN}\b'),
    re.compile(rf'\b{GENE_TOKEN}\s*[:(, ]\s*[A-Z]\d+[A-Z*X]\b'),
    re.compile(rf'\b{GENE_TOKEN}\s*\(\s*[A-Z]\d+[A-Z*X]\s*\)'),
    re.compile(rf'[Mm][Uu][Tt][Aa][Tt][Ii][Oo][Nn]\s+[Ii][Nn]\s+(?:[Tt][Hh][Ee]\s+)?{GENE_TOKEN}\b'),
    re.compile(rf'[Vv][Aa][Rr][Ii][Aa][Nn][Tt]\s+[Ii][Nn]\s+(?:[Tt][Hh][Ee]\s+)?{GENE_TOKEN}\b'),
    re.compile(rf'[Ii][Nn]\s+(?:[Tt][Hh][Ee]\s+)?{GENE_TOKEN}\s+[Gg][Ee][Nn][Ee]\b'),
    re.compile(rf'\b{GENE_TOKEN}\s+[Mm][Uu][Tt][Aa][Tt][Ii][Oo][Nn]s?\b'),
    re.compile(rf'\b{GENE_TOKEN}\s+[Vv][Aa][Rr][Ii][Aa][Nn][Tt]s?\b'),
]

GENE_MENTION_PATTERN = re.compile(rf'\b{GENE_TOKEN}\b')

AA_PATTERN = re.compile(
    r'^(?:p\.)?'
    r'(?:[ARNDCQEGHILKMFPSTWYV]|Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)'
    r'\d+'
    r'(?:[ARNDCQEGHILKMFPSTWYV*X]|Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter|fs|del|dup|ins)?$',
    re.IGNORECASE
)

PEDIGREE_PATTERN = re.compile(r'^[A-Z]{1,2}-[IVX]+\d+$|^[IVX]+-\d+$')
CODING_PATTERN = re.compile(r'^(?:c\.)?[ACGT]\d+[ACGT>]?(?:DEL|DUP|INS)?$', re.IGNORECASE)
IVS_PATTERN = re.compile(r'^IVS\d+', re.IGNORECASE)
TRAILING_HYPHEN = re.compile(r'-$')
CH_COMPLEMENT = re.compile(r'^CH\d+$', re.IGNORECASE)

GENE_CONTEXT_PHRASES = [
    r'{gene}\s+gene\b',
    r'screening\s+of\s+{gene}\b',
    r'{gene}\s+mutation',
    r'{gene}\s+variant',
    r'mutation\s+in\s+{gene}\b',
    r'variant\s+in\s+{gene}\b',
]


# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class GeneScore:
    """Tracks scoring details for a gene candidate."""
    gene: str
    total_score: int = 0
    cited_sentence_comention: bool = False
    question_comention: bool = False
    rationale_comention: bool = False
    answer_comention: bool = False
    gene_context_phrase: bool = False
    multi_source: bool = False
    sources_found: Set[str] = field(default_factory=set)

    def compute_total(self) -> int:
        """Compute total score based on components."""
        score = 0
        if self.cited_sentence_comention:
            score += 100
        if self.question_comention:
            score += 90
        if self.rationale_comention:
            score += 80
        if self.answer_comention:
            score += 40
        if self.gene_context_phrase:
            score += 30
        if len(self.sources_found) >= 2:
            score += 10
            self.multi_source = True
        self.total_score = score
        return score


@dataclass
class RecoveryResult:
    """Result of gene recovery for a record."""
    gene: Optional[str]
    tier: str
    confidence: str
    score: int = 0


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def is_valid_gene_symbol(s: str) -> bool:
    """Check if string could be a valid gene symbol."""
    if s != s.upper():
        return False
    if s in NON_GENES:
        return False
    if not (2 <= len(s) <= 15):
        return False
    if AA_PATTERN.match(s):
        return False
    if CODING_PATTERN.match(s):
        return False
    if PEDIGREE_PATTERN.match(s):
        return False
    if IVS_PATTERN.match(s):
        return False
    if TRAILING_HYPHEN.search(s):
        return False
    if CH_COMPLEMENT.match(s):
        return False
    if not any(c.isalpha() for c in s):
        return False
    return True


# ============================================================================
# MUTATION ALIAS GENERATION
# ============================================================================

def generate_mutation_aliases(mutation: str) -> Set[str]:
    """Generate variations of a mutation string for flexible matching."""
    aliases = set()
    mutation = mutation.strip()

    if not mutation:
        return aliases

    aliases.add(mutation)

    if not mutation.startswith(('p.', 'c.', 'g.', 'r.')):
        aliases.add(f'p.{mutation}')
        aliases.add(f'c.{mutation}')

    if mutation.startswith(('p.', 'c.', 'g.', 'r.')):
        stripped = mutation[2:]
        aliases.add(stripped)

    stop_variants = []
    for alias in list(aliases):
        if alias.endswith('*'):
            stop_variants.append(alias[:-1] + 'X')
            stop_variants.append(alias[:-1] + 'Ter')
        if alias.endswith('X') and not alias.endswith('Ter'):
            stop_variants.append(alias[:-1] + '*')
            stop_variants.append(alias[:-1] + 'Ter')
        if alias.endswith('Ter'):
            stop_variants.append(alias[:-3] + '*')
            stop_variants.append(alias[:-3] + 'X')

        mid_stop_match = re.search(r'([A-Z])(\d+)([*X]|Ter)$', alias)
        if mid_stop_match:
            base = alias[:mid_stop_match.start()] + mid_stop_match.group(1) + mid_stop_match.group(2)
            stop_variants.append(base + '*')
            stop_variants.append(base + 'X')
            stop_variants.append(base + 'Ter')

    aliases.update(stop_variants)
    return aliases


# ============================================================================
# GROUNDED BLOCK BUILDING (uses pre-extracted cited_sentences)
# ============================================================================

def build_grounded_blocks(record: Dict) -> Dict[str, List[str]]:
    """
    Build grounded text blocks from record.
    Uses pre-extracted 'cited_sentences' field (no article lookup).
    """
    blocks = {
        'cited_sentences': record.get('cited_sentences', []),
        'rationale': [],
        'question': [],
        'answers': [],
    }

    for field in ['phase1_rationale', 'phase2_rationale']:
        text = record.get(field, '')
        if text and text.strip():
            blocks['rationale'].append(text.strip())

    question = record.get('question', '')
    if question and question.strip():
        blocks['question'].append(question.strip())

    for field in ['phase1_answer', 'phase2_answer']:
        text = record.get(field, '')
        if text and text.strip():
            blocks['answers'].append(text.strip())

    return blocks


# ============================================================================
# CANDIDATE GENERATION
# ============================================================================

def extract_genes_with_patterns(text: str) -> List[str]:
    """Extract gene symbols using high-precision patterns."""
    candidates = []
    seen = set()

    for pattern in GENE_PATTERNS:
        for match in pattern.finditer(text):
            try:
                gene = match.group('gene').upper()
            except (IndexError, AttributeError):
                continue
            if gene in seen:
                continue
            if not is_valid_gene_symbol(gene):
                continue
            seen.add(gene)
            candidates.append(gene)

    return candidates


def find_all_gene_tokens(text: str) -> List[str]:
    """Find ALL HGNC-like tokens, filtered by is_valid_gene_symbol."""
    genes = []
    seen = set()

    for match in GENE_MENTION_PATTERN.finditer(text):
        try:
            gene = match.group('gene')
        except (IndexError, AttributeError):
            continue
        if gene != gene.upper():
            continue
        if gene in seen:
            continue
        if not is_valid_gene_symbol(gene):
            continue
        seen.add(gene)
        genes.append(gene)

    return genes


def generate_candidates(blocks: Dict[str, List[str]]) -> List[str]:
    """Generate gene candidates from grounded blocks."""
    all_text = ' '.join(
        ' '.join(texts) for texts in blocks.values()
    )

    candidates = extract_genes_with_patterns(all_text)
    if candidates:
        return candidates

    return find_all_gene_tokens(all_text)


# ============================================================================
# SCORING SYSTEM
# ============================================================================

def check_comention(text: str, gene: str, mutation_aliases: Set[str]) -> bool:
    """Check if gene and any mutation alias co-occur in the same text."""
    gene_pattern = re.compile(r'\b' + re.escape(gene) + r'\b', re.IGNORECASE)
    if not gene_pattern.search(text):
        return False

    text_normalized = re.sub(r'\s+', '', text)

    for alias in mutation_aliases:
        if not alias:
            continue
        try:
            alias_pattern = re.compile(re.escape(alias), re.IGNORECASE)
            if alias_pattern.search(text):
                return True
            if alias_pattern.search(text_normalized):
                return True
        except re.error:
            continue

    return False


def check_gene_context_phrase(text: str, gene: str) -> bool:
    """Check for gene-context phrases."""
    for phrase_template in GENE_CONTEXT_PHRASES:
        pattern_str = phrase_template.format(gene=re.escape(gene))
        try:
            if re.search(pattern_str, text, re.IGNORECASE):
                return True
        except re.error:
            continue
    return False


def score_candidate(
    gene: str,
    blocks: Dict[str, List[str]],
    mutation_aliases: Set[str]
) -> GeneScore:
    """Score a gene candidate based on grounded evidence."""
    score = GeneScore(gene=gene)

    for sent in blocks.get('cited_sentences', []):
        if check_comention(sent, gene, mutation_aliases):
            score.cited_sentence_comention = True
            score.sources_found.add('cited_sentence')
            break

    for q in blocks.get('question', []):
        if check_comention(q, gene, mutation_aliases):
            score.question_comention = True
            score.sources_found.add('question')
            break

    for r in blocks.get('rationale', []):
        if check_comention(r, gene, mutation_aliases):
            score.rationale_comention = True
            score.sources_found.add('rationale')
            break

    for a in blocks.get('answers', []):
        if check_comention(a, gene, mutation_aliases):
            score.answer_comention = True
            score.sources_found.add('answer')
            break

    all_text = ' '.join(' '.join(texts) for texts in blocks.values())
    if check_gene_context_phrase(all_text, gene):
        score.gene_context_phrase = True

    score.compute_total()
    return score


# ============================================================================
# DECISION LOGIC
# ============================================================================

def apply_decision_thresholds(
    scored_candidates: Dict[str, GeneScore]
) -> RecoveryResult:
    """Apply decision thresholds and ambiguity guard."""
    if not scored_candidates:
        return RecoveryResult(
            gene=None,
            tier='unresolved_no_candidates',
            confidence='unresolved',
        )

    sorted_genes = sorted(
        scored_candidates.items(),
        key=lambda x: x[1].total_score,
        reverse=True
    )

    top_gene, top_score = sorted_genes[0]

    if len(sorted_genes) >= 2:
        second_gene, second_score = sorted_genes[1]
        score_diff = top_score.total_score - second_score.total_score

        if score_diff < 20:
            if not top_score.cited_sentence_comention:
                return RecoveryResult(
                    gene=None,
                    tier='unresolved_close_scores',
                    confidence='unresolved',
                    score=top_score.total_score,
                )

    if len(top_gene) <= 4 and top_score.cited_sentence_comention:
        has_question = top_score.question_comention
        has_gene_context = top_score.gene_context_phrase
        has_multi_source = len(top_score.sources_found) >= 2
        has_rationale = top_score.rationale_comention

        if not (has_question or has_gene_context or has_multi_source or has_rationale):
            return RecoveryResult(
                gene=None,
                tier='unresolved_short_cited_only',
                confidence='unresolved',
                score=top_score.total_score,
            )

    if top_score.total_score >= 90:
        return RecoveryResult(
            gene=top_gene,
            tier='resolved_high_comention',
            confidence='high',
            score=top_score.total_score,
        )
    elif top_score.total_score >= 80 and len(top_score.sources_found) >= 2:
        return RecoveryResult(
            gene=top_gene,
            tier='resolved_medium_comention',
            confidence='medium',
            score=top_score.total_score,
        )
    else:
        return RecoveryResult(
            gene=None,
            tier='unresolved_low_score',
            confidence='unresolved',
            score=top_score.total_score,
        )


# ============================================================================
# TIER 1 FALLBACKS
# ============================================================================

def tier1_field_search(record: Dict) -> Optional[str]:
    """Tier 1: Search data fields for gene mentions using patterns."""
    fields_text = ' '.join([
        record.get('question', ''),
        record.get('phase1_answer', ''),
        record.get('phase1_rationale', ''),
        record.get('phase2_answer', ''),
        record.get('phase2_rationale', ''),
    ])

    genes = extract_genes_with_patterns(fields_text)
    return genes[0] if genes else None


# ============================================================================
# MAIN RECOVERY FUNCTION
# ============================================================================

def recover_gene(record: Dict, pmid_gene_counts: Dict[str, Dict[str, int]] = None) -> RecoveryResult:
    """
    Main gene recovery function.

    Tier order:
    - Tier 0: Grounded scoring
    - Tier 1: Field search (backup)
    - Tier 1b: Single-gene PMID fallback (if pmid_gene_counts provided)
    """
    pmid = str(record.get('pmid', ''))
    mutation = record.get('mutation', '')

    mutation_aliases = generate_mutation_aliases(mutation)
    blocks = build_grounded_blocks(record)
    candidates = generate_candidates(blocks)

    if not candidates:
        tier1_gene = tier1_field_search(record)
        if tier1_gene:
            return RecoveryResult(
                gene=tier1_gene,
                tier='resolved_medium_pattern',
                confidence='medium',
            )

        # Tier 1b: Single-gene PMID fallback
        if pmid_gene_counts and pmid in pmid_gene_counts:
            gene_counts = pmid_gene_counts[pmid]
            if len(gene_counts) == 1:
                single_gene = list(gene_counts.keys())[0]
                return RecoveryResult(
                    gene=single_gene,
                    tier='resolved_high_pmid_single',
                    confidence='high',
                )

        return RecoveryResult(
            gene=None,
            tier='unresolved_no_candidates',
            confidence='unresolved',
        )

    scored = {}
    for gene in candidates:
        scored[gene] = score_candidate(gene, blocks, mutation_aliases)

    result = apply_decision_thresholds(scored)

    if result.gene is None:
        tier1_gene = tier1_field_search(record)
        if tier1_gene:
            return RecoveryResult(
                gene=tier1_gene,
                tier='resolved_medium_pattern',
                confidence='medium',
                score=result.score,
            )

        # Tier 1b: Single-gene PMID fallback
        if pmid_gene_counts and pmid in pmid_gene_counts:
            gene_counts = pmid_gene_counts[pmid]
            if len(gene_counts) == 1:
                single_gene = list(gene_counts.keys())[0]
                return RecoveryResult(
                    gene=single_gene,
                    tier='resolved_high_pmid_single',
                    confidence='high',
                    score=result.score,
                )

    return result


# Global variable for worker processes (set via initializer)
_pmid_gene_counts = None


def _init_worker(pmid_gene_counts):
    """Initialize worker process with shared data."""
    global _pmid_gene_counts
    _pmid_gene_counts = pmid_gene_counts


def process_record(record: Dict) -> Dict:
    """Process a single record for gene recovery."""
    result = recover_gene(record, _pmid_gene_counts)

    record['recovered_gene'] = result.gene
    record['recovery_tier'] = result.tier
    record['recovery_confidence'] = result.confidence
    record['recovery_score'] = result.score

    return record


def process_batch(records: List[Dict]) -> List[Dict]:
    """Process a batch of records."""
    return [process_record(r) for r in records]


# ============================================================================
# MAIN
# ============================================================================

def main():
    global INPUT_PATH, OUTPUT_PATH, NUM_WORKERS
    parser = argparse.ArgumentParser(description="Heuristic gene recovery over cited sentences (no API calls)")
    parser.add_argument("--input", default=INPUT_PATH, help="m2tqa_with_cited.jsonl (extract_cited_sentences.py output)")
    parser.add_argument("--output", default=OUTPUT_PATH)
    parser.add_argument("--workers", type=int, default=NUM_WORKERS)
    args = parser.parse_args()
    INPUT_PATH, OUTPUT_PATH, NUM_WORKERS = args.input, args.output, args.workers

    print("=" * 70)
    print("STEP 2: RUN GENE RECOVERY (PARALLEL)")
    print("=" * 70)

    print(f"\nUsing {NUM_WORKERS} workers")
    print(f"Input: {INPUT_PATH}")
    print(f"Output: {OUTPUT_PATH}")

    # First pass: Build PMID gene counts for Tier 1b fallback
    print("\nPass 1: Building PMID gene counts for Tier 1b fallback...")
    pmid_gene_counts = defaultdict(lambda: defaultdict(int))
    record_count = 0

    with open(INPUT_PATH) as f:
        for line in f:
            record = json.loads(line)
            pmid = str(record.get('pmid', ''))
            record_count += 1

            fields_text = ' '.join([
                record.get('question', ''),
                record.get('phase1_answer', ''),
                record.get('phase1_rationale', ''),
                record.get('phase2_answer', ''),
                record.get('phase2_rationale', ''),
            ])

            genes = extract_genes_with_patterns(fields_text)
            for gene in genes:
                pmid_gene_counts[pmid][gene] += 1

            if record_count % 100000 == 0:
                print(f"  Scanned {record_count}...")

    print(f"  Scanned {record_count} records")
    print(f"  Found genes for {len(pmid_gene_counts)} PMIDs")

    # Convert to regular dicts for serialization
    pmid_gene_counts = {k: dict(v) for k, v in pmid_gene_counts.items()}

    # Second pass: Process all records
    print("\nPass 2: Processing records...")

    # Load all records into memory for batching
    print("  Loading records...")
    records = []
    with open(INPUT_PATH) as f:
        for line in f:
            records.append(json.loads(line))

    print(f"  Loaded {len(records)} records")

    # Process in batches using multiprocessing
    batch_size = 1000
    batches = [records[i:i + batch_size] for i in range(0, len(records), batch_size)]
    print(f"  Created {len(batches)} batches of ~{batch_size} records")

    print(f"  Processing with {NUM_WORKERS} workers...")

    # Collect all processed records
    all_results = []

    with Pool(NUM_WORKERS, initializer=_init_worker, initargs=(pmid_gene_counts,)) as pool:
        processed = 0

        for batch_results in pool.imap(process_batch, batches):
            all_results.extend(batch_results)

            processed += len(batch_results)
            if processed % 50000 == 0:
                print(f"    Processed {processed}...")

    # =========================================================================
    # Pass 3: Propagate genes within (PMID, mutation) groups
    # =========================================================================
    print("\nPass 3: Propagating genes within (PMID, mutation) groups...")

    # Group records by (PMID, mutation)
    groups = defaultdict(list)
    for i, record in enumerate(all_results):
        key = (str(record['pmid']), record['mutation'])
        groups[key].append(i)

    propagated_count = 0
    ambiguous_groups = 0

    for key, indices in groups.items():
        if len(indices) <= 1:
            continue  # Single record, nothing to propagate

        # Find high-confidence seeds (score >= 90 or resolved with high confidence)
        seeds = []
        unresolved_indices = []

        for idx in indices:
            record = all_results[idx]
            score = record.get('recovery_score', 0)
            gene = record.get('recovered_gene')
            confidence = record.get('recovery_confidence')

            if gene and (score >= 90 or confidence == 'high'):
                seeds.append(gene)
            elif gene is None or confidence == 'unresolved':
                unresolved_indices.append(idx)

        if not seeds or not unresolved_indices:
            continue  # No seeds or nothing to propagate to

        # Check if all seeds agree
        unique_seed_genes = set(seeds)

        if len(unique_seed_genes) == 1:
            # All seeds agree - propagate to unresolved records
            propagate_gene = list(unique_seed_genes)[0]

            for idx in unresolved_indices:
                all_results[idx]['recovered_gene'] = propagate_gene
                all_results[idx]['recovery_tier'] = 'resolved_high_sibling'
                all_results[idx]['recovery_confidence'] = 'high'
                # Keep original score for debugging
                propagated_count += 1
        else:
            # Seeds disagree - mark as ambiguous, don't propagate
            ambiguous_groups += 1

    print(f"  Propagated gene to {propagated_count} records from sibling seeds")
    print(f"  Ambiguous groups (conflicting seeds): {ambiguous_groups}")

    # =========================================================================
    # Write final results and compute stats
    # =========================================================================
    print("\nWriting final results...")

    tier_counts = defaultdict(int)
    confidence_counts = defaultdict(int)

    with open(OUTPUT_PATH, 'w') as f_out:
        for record in all_results:
            tier_counts[record['recovery_tier']] += 1
            confidence_counts[record['recovery_confidence']] += 1
            f_out.write(json.dumps(record) + '\n')

    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"Records processed: {len(records)}")

    resolved = sum(c for t, c in tier_counts.items() if t.startswith('resolved'))
    unresolved = sum(c for t, c in tier_counts.items() if t.startswith('unresolved'))

    print(f"\nResolved: {resolved} ({resolved / len(records) * 100:.1f}%)")
    print(f"Unresolved: {unresolved} ({unresolved / len(records) * 100:.1f}%)")

    print("\nTier breakdown:")
    for tier, count in sorted(tier_counts.items(), key=lambda x: -x[1]):
        print(f"  {tier:35} {count:7} ({count / len(records) * 100:.1f}%)")

    print("\nConfidence breakdown:")
    for conf, count in sorted(confidence_counts.items(), key=lambda x: -x[1]):
        print(f"  {conf:15} {count:7} ({count / len(records) * 100:.1f}%)")

    print(f"\nOutput written to: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
