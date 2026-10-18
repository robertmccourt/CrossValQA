#!/usr/bin/env python3
"""
Tiered Gene Anchoring Pipeline v3

Tier 1: Full-text NER (GNormPlus) with proximity matching
Tier 2: Abstract-only - mutations in abstract matched to PubTator gene offsets
Tier 3: Most-mentioned gene fallback

v3 Improvements:
- Directory-based GNormPlus output parsing (handles parallel batch processing)
- Batch integrity checking
- Proper PubTator offset handling (Title + newline + Body concatenation)
- Species-agnostic anchoring with ncbi_taxon_id for downstream filtering
- Memory-efficient article streaming with ijson

Anchoring step 5 (last). Inputs:
- phase_4/output/m2tqa_gene_recovered.jsonl (run_gene_recovery.py)
- anchoring/output/pubtator_annotations.jsonl (fetch_pubtator_bioc.py)
- anchoring/output/ner_output/batch_*.out (fetch_pubtator_batches.py, or GNormPlus;
  PubTator text format, see load_gnormplus_annotations)
- the source article corpus {pmid: {"article": ...}} (config.SOURCE_ARTICLES_FILE)

Outputs:
- anchoring/output/m2tqa_tiered_anchored_v3.jsonl: every input record, in order, with
  anchor_gene_id, anchor_gene_name, anchor_tier, anchor_confidence, ncbi_taxon_id and,
  when set, anchor_distance and vep_input_string appended
  (input of uniprot/scripts/master_grounding_pipeline.py)
- anchoring/output/tiered_stats_v3.json

Each record's anchor depends only on its own (pmid, mutation) and that PMID's
annotations/article, so running on a subset of records reproduces the same
per-record output; only the stats file depends on which records are processed.
"""

import argparse
import json
import re
import sys
import ijson
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional, Generator
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import config  # noqa: E402

# Paths (overridable on the command line, see main())
ANCHORING_OUTPUT_DIR = config.PHASE4_DIR / "anchoring" / "output"
PUBTATOR_FILE = ANCHORING_OUTPUT_DIR / "pubtator_annotations.jsonl"
QA_FILE = config.PHASE4_GENE_RECOVERED_FILE
ARTICLES_FILE = config.SOURCE_ARTICLES_FILE
GNORMPLUS_DIR = ANCHORING_OUTPUT_DIR / "ner_output"
OUTPUT_FILE = ANCHORING_OUTPUT_DIR / "m2tqa_tiered_anchored_v3.jsonl"
STATS_FILE = ANCHORING_OUTPUT_DIR / "tiered_stats_v3.json"

# Expected batch count (0-30 = 31 batches)
EXPECTED_BATCHES = 31

# Maximum character distance for proximity matching
MAX_PROXIMITY_DISTANCE = 2000

# False positive blacklist - common GNormPlus/text mining errors
GENE_BLACKLIST = {
    # Prediction tools mistaken for genes
    "SIFT", "sift", "PolyPhen", "POLYPHEN", "polyphen", "PolyPhen-2",
    "CADD", "REVEL", "GERP", "MutationTaster", "PROVEAN",
    # Clinical terms
    "VUS", "VOUS", "LP", "P", "LB", "B",  # ACMG classifications
    "IMPACT", "impact",
    # Common false positives
    "DNA", "RNA", "PCR", "SNP", "CNV", "LOH",
    "WT", "MUT", "DEL", "INS", "DUP",
}

# Amino acid mappings for mutation pattern matching
AA_3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': '*', 'Stop': '*'
}

AA_1TO3 = {v: k for k, v in AA_3TO1.items() if v != '*'}
AA_1TO3['*'] = 'Ter'


@dataclass
class AnchorResult:
    gene_id: Optional[str]
    gene_name: Optional[str]
    tier: str
    confidence: str
    distance: Optional[int] = None
    vep_input_string: Optional[str] = None
    ncbi_taxon_id: Optional[str] = None


def check_batch_integrity(gnormplus_dir: Path) -> list:
    """
    Check if all expected batch output files are present.
    Returns list of missing batch numbers.
    """
    missing = []
    for i in range(EXPECTED_BATCHES):
        # Check common output extensions
        batch_patterns = [
            f"batch_{i:04d}.out",
            f"batch_{i:04d}.txt.out",
            f"batch_{i:04d}.PubTator",
        ]
        found = False
        for pattern in batch_patterns:
            if (gnormplus_dir / pattern).exists():
                found = True
                break
        if not found:
            missing.append(i)
    return missing


def load_gnormplus_annotations(gnormplus_dir: Path) -> dict:
    """
    Load GNormPlus annotations from all output files in directory.
    Handles parallel batch processing output.

    GNormPlus output format (tab-separated):
    PMID    Start   End   Mention   Type   ConceptID

    Returns dict: {pmid: {gene_id: {"name": str, "offsets": [(start, end)], "taxon_id": str}}}
    """
    # Use nested defaultdicts for clean initialization
    fulltext_ner = defaultdict(lambda: defaultdict(lambda: {"name": "", "offsets": [], "taxon_id": "unknown"}))

    if not gnormplus_dir.exists():
        print(f"GNormPlus output directory not found: {gnormplus_dir}")
        return {}

    # Find all output files (various extensions)
    out_files = list(gnormplus_dir.glob("*.out")) + \
                list(gnormplus_dir.glob("*.PubTator")) + \
                list(gnormplus_dir.glob("*.txt"))

    # Filter out input files if any ended up here
    out_files = [f for f in out_files if "batch_" in f.name or "output" in f.name.lower()]

    if not out_files:
        print("No GNormPlus output files found")
        return {}

    print(f"Found {len(out_files)} GNormPlus output files")

    # Check batch integrity
    missing = check_batch_integrity(gnormplus_dir)
    if missing:
        print(f"WARNING: Missing {len(missing)} batch files: {missing[:10]}{'...' if len(missing) > 10 else ''}")
    else:
        print("All expected batch files present")

    total_genes = 0
    total_pmids = set()

    for output_file in tqdm(out_files, desc="Loading GNormPlus"):
        with open(output_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')

                # Format: PMID \t Start \t End \t Mention \t Type \t ConceptID
                if len(parts) < 6:
                    continue

                pmid = parts[0]
                entity_type = parts[4]

                # Only process Gene entities
                if entity_type != "Gene":
                    continue

                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue

                name = parts[3]
                concept_id = parts[5]

                # Skip blacklisted terms
                if name in GENE_BLACKLIST or concept_id in GENE_BLACKLIST:
                    continue

                # Split Taxon from Gene if present (format: "9606:672")
                taxon_id = "unknown"
                gene_id = concept_id
                if ":" in concept_id:
                    taxon_id, gene_id = concept_id.split(":", 1)

                # Skip if gene_id is empty or blacklisted
                if not gene_id or gene_id in GENE_BLACKLIST:
                    continue

                fulltext_ner[pmid][gene_id]["offsets"].append((start, end))
                fulltext_ner[pmid][gene_id]["name"] = name
                fulltext_ner[pmid][gene_id]["taxon_id"] = taxon_id

                total_genes += 1
                total_pmids.add(pmid)

    print(f"Loaded {total_genes:,} gene annotations across {len(total_pmids):,} PMIDs")

    # Convert defaultdict to regular dict for cleaner handling
    return {pmid: dict(genes) for pmid, genes in fulltext_ner.items()}


def load_pubtator() -> dict:
    """Load PubTator annotations indexed by PMID."""
    annotations = {}
    with open(PUBTATOR_FILE) as f:
        for line in f:
            r = json.loads(line)
            pmid = r.get("pmid")
            if pmid:
                annotations[str(pmid)] = r
    return annotations


def stream_articles(target_pmids: set) -> Generator[tuple, None, None]:
    """
    Memory-efficient article streaming using ijson.
    Yields (pmid, article_text) tuples for target PMIDs only.
    """
    print("Streaming articles with ijson (memory-efficient)...")
    with open(ARTICLES_FILE, 'rb') as f:
        parser = ijson.kvitems(f, '')
        for pmid, record in parser:
            if pmid in target_pmids:
                if isinstance(record, dict) and "article" in record:
                    yield pmid, record["article"]


def build_article_index(target_pmids: set) -> dict:
    """
    Build article index for target PMIDs only.
    More memory efficient than loading all 4.5GB.
    """
    print(f"Loading articles for {len(target_pmids):,} target PMIDs...")
    articles = {}
    count = 0
    for pmid, article_text in stream_articles(target_pmids):
        articles[pmid] = article_text
        count += 1
        if count % 10000 == 0:
            print(f"  Loaded {count:,} articles...")
    print(f"Loaded {len(articles):,} articles")
    return articles


def get_pubtator_text(article_text: str) -> str:
    """
    Reconstruct the text as GNormPlus sees it.
    GNormPlus treats input as: Title + newline + Body
    The offsets are based on this concatenated view.

    Since our input format is:
    PMID|t|Title
    PMID|a|Body

    We need to simulate: Title\nBody
    """
    if not article_text:
        return ""

    # prepare_ner_input.py extracted title as first line, body as rest
    # The article_text from JSON is the raw full text
    lines = article_text.split('\n', 1)
    title = lines[0].strip() if lines else ""
    body = lines[1] if len(lines) > 1 else ""

    # GNormPlus concatenates with single newline
    return title + "\n" + body


def normalize_mutation_text(text: str) -> str:
    """Normalize text for mutation matching."""
    # Convert Greek delta to 'del'
    text = text.replace('Δ', 'del').replace('δ', 'del')
    # Normalize arrows
    text = text.replace('→', '>').replace('➔', '>').replace('->', '>')
    return text


def find_mutation_in_text(text: str, mutation: str) -> Optional[int]:
    """
    Find character position of mutation in text.
    Handles various mutation formats including indels and Greek symbols.

    Note: For GNormPlus alignment, the text should be in PubTator format
    (Title + newline + Body) so offsets match.
    """
    if not text or not mutation:
        return None

    # Normalize both text and mutation
    text_norm = normalize_mutation_text(text)
    text_lower = text_norm.lower()
    mutation_norm = normalize_mutation_text(mutation)
    mutation_lower = mutation_norm.lower()

    # 1. Exact match
    pos = text_lower.find(mutation_lower)
    if pos >= 0:
        return pos

    # 2. Without prefix (p.R248Q -> R248Q, c.743G>A -> 743G>A)
    if mutation.startswith(('p.', 'c.', 'g.', 'r.', 'm.')):
        short = mutation[2:]
        short_norm = normalize_mutation_text(short)
        pos = text_norm.find(short_norm)
        if pos >= 0:
            return pos
        pos = text_lower.find(short_norm.lower())
        if pos >= 0:
            return pos

    # 3. Three-letter amino acid codes: p.Arg248Gln
    match = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', mutation)
    if match:
        aa1, pos_num, aa2 = match.groups()

        # Try exact 3-letter format
        pattern_3 = f"{aa1}{pos_num}{aa2}"
        idx = text.find(pattern_3)
        if idx >= 0:
            return idx

        # Try with arrow: Arg248->Gln or Arg648>Gln
        pattern_arrow = rf"{aa1}\s*{pos_num}\s*[-=]?>\s*{aa2}"
        m = re.search(pattern_arrow, text, re.IGNORECASE)
        if m:
            return m.start()

        # Convert to single letter: R248Q
        aa1_single = AA_3TO1.get(aa1, aa1[0])
        aa2_single = AA_3TO1.get(aa2, aa2[0])
        pattern_1 = f"{aa1_single}{pos_num}{aa2_single}"
        idx = text.find(pattern_1)
        if idx >= 0:
            return idx

    # 4. Single-letter amino acid: p.R248Q
    match = re.search(r'p\.([A-Z])(\d+)([A-Z*])', mutation)
    if match:
        aa1, pos_num, aa2 = match.groups()

        # Try single letter format
        pattern_1 = f"{aa1}{pos_num}{aa2}"
        idx = text.find(pattern_1)
        if idx >= 0:
            return idx

        # Try with arrow (escape aa2 in case it's * for stop codon)
        aa2_escaped = re.escape(aa2)
        pattern_arrow = rf"{aa1}\s*{pos_num}\s*[-=]?>\s*{aa2_escaped}"
        m = re.search(pattern_arrow, text, re.IGNORECASE)
        if m:
            return m.start()

        # Convert to 3-letter
        aa1_3 = AA_1TO3.get(aa1, aa1)
        aa2_3 = AA_1TO3.get(aa2, aa2)
        if len(aa1_3) == 3 and len(aa2_3) == 3:
            pattern_3 = f"{aa1_3}{pos_num}{aa2_3}"
            idx = text.find(pattern_3)
            if idx >= 0:
                return idx

    # 5. Indel patterns: p.Gly12_Gly13del, c.1234_1236del
    if 'del' in mutation_lower or 'ins' in mutation_lower or 'dup' in mutation_lower:
        # Extract position numbers
        positions = re.findall(r'(\d+)', mutation)
        if positions:
            # Search for any position number near del/ins/dup
            for pos_num in positions:
                pattern = rf"{pos_num}\s*(?:_\d+)?\s*(?:del|ins|dup)"
                m = re.search(pattern, text_lower)
                if m:
                    return m.start()

    # 6. Frameshift: p.Arg248fs, p.R248fs
    if 'fs' in mutation_lower:
        match = re.search(r'([A-Z][a-z]{0,2})(\d+)fs', mutation, re.IGNORECASE)
        if match:
            aa, pos_num = match.groups()
            pattern = rf"{aa}.*{pos_num}.*fs"
            m = re.search(pattern, text_lower)
            if m:
                return m.start()

    # 7. rsID
    if mutation.lower().startswith('rs'):
        idx = text_lower.find(mutation_lower)
        if idx >= 0:
            return idx

    # 8. Genomic coordinates: chr17:7577120:C:T
    if ':' in mutation and any(c in mutation for c in ['chr', 'g.']):
        # Try to find any part of the coordinate
        parts = mutation.replace('chr', '').split(':')
        for part in parts:
            if part.isdigit() and len(part) > 4:  # Likely a position
                idx = text.find(part)
                if idx >= 0:
                    return idx

    return None


def get_gene_offsets(pubtator_record: dict) -> dict:
    """Get all gene mentions with their offsets from PubTator record."""
    genes = {}
    # Get paper-level species info (PubTator provides species at paper level)
    paper_species = pubtator_record.get("species", [])
    default_taxon = "9606" if "9606" in paper_species else (paper_species[0] if paper_species else None)

    for g in pubtator_record.get("genes", []):
        gene_id = g.get("ncbi_gene_id", "")
        gene_name = g.get("name", "") or g.get("text", "")

        if not gene_id:
            continue

        # Skip blacklisted
        if gene_name in GENE_BLACKLIST or gene_id in GENE_BLACKLIST:
            continue

        if gene_id not in genes:
            genes[gene_id] = {
                "name": gene_name,
                "taxon_id": default_taxon,
                "offsets": []
            }

        for offset in g.get("offsets", []):
            if len(offset) >= 2:
                start = offset[0]
                end = start + offset[1]
                genes[gene_id]["offsets"].append((start, end))

    return genes


def find_closest_gene(mutation_pos: int, gene_offsets: dict, max_dist: int = MAX_PROXIMITY_DISTANCE) -> tuple:
    """
    Finds the biologically relevant gene by calculating the shortest
    distance between a mutation mention and a Gene ID mention.

    Returns (gene_id, gene_name, distance, taxon_id) or (None, None, None, None) if no gene within max_dist.
    """
    closest_gene, closest_name, closest_taxon, min_distance = None, None, None, float('inf')

    for gene_id, info in gene_offsets.items():
        for start, end in info["offsets"]:
            # Distance is 0 if mutation is inside the gene mention
            if mutation_pos < start:
                dist = start - mutation_pos
            elif mutation_pos > end:
                dist = mutation_pos - end
            else:
                dist = 0

            if dist < min_distance:
                min_distance = dist
                closest_gene = gene_id
                closest_name = info["name"]
                closest_taxon = info.get("taxon_id")

    # Only return if within a reasonable 'window' (e.g., same section/paragraph)
    if min_distance <= max_dist:
        return closest_gene, closest_name, min_distance, closest_taxon
    return None, None, None, None


def get_most_mentioned_gene(pubtator_record: dict) -> tuple:
    """Get the most frequently mentioned gene. Returns (gene_id, gene_name, taxon_id)."""
    gene_counts = defaultdict(int)
    gene_names = {}

    # Get paper-level species info
    paper_species = pubtator_record.get("species", [])
    default_taxon = "9606" if "9606" in paper_species else (paper_species[0] if paper_species else None)

    for g in pubtator_record.get("genes", []):
        gene_id = g.get("ncbi_gene_id", "")
        gene_name = g.get("name", "") or g.get("text", "")

        if gene_id and gene_name not in GENE_BLACKLIST and gene_id not in GENE_BLACKLIST:
            gene_counts[gene_id] += 1
            gene_names[gene_id] = gene_name

    if not gene_counts:
        return None, None, None

    most_common = max(gene_counts.items(), key=lambda x: x[1])[0]
    return most_common, gene_names.get(most_common, ""), default_taxon


def extract_abstract(full_text: str) -> str:
    """
    Extract abstract portion from full text.
    PubTator offsets are based on title + abstract.
    """
    # Look for common abstract markers
    abstract_markers = ['abstract', 'summary', 'background']
    text_lower = full_text.lower()

    for marker in abstract_markers:
        pos = text_lower.find(marker)
        if pos != -1 and pos < 500:
            # Find end of abstract (next section)
            section_markers = ['\nintroduction', '\nmethods', '\nmaterials', '\nresults']
            end_pos = len(full_text)
            for end_marker in section_markers:
                epos = text_lower.find(end_marker, pos)
                if epos != -1:
                    end_pos = min(end_pos, epos)
            return full_text[pos:end_pos]

    # Fallback: first 3000 chars (typical abstract length)
    return full_text[:3000]


def format_vep_input(mutation: str, gene_id: Optional[str] = None) -> Optional[str]:
    """
    Format mutation for VEP input.

    Returns string in format suitable for VEP:
    - rsID: rs12345
    - HGVS: ENST00000269305.8:c.743G>A
    - Genomic: chr:pos:ref:alt

    Returns None if mutation format is not VEP-compatible.
    """
    if not mutation:
        return None

    # rsID - ready for VEP
    if mutation.lower().startswith('rs'):
        return mutation

    # c. notation - needs transcript context
    if mutation.startswith('c.'):
        # Would need gene_id -> canonical transcript mapping
        # For now, return the c. notation with a placeholder
        return f"NEEDS_TRANSCRIPT:{mutation}"

    # Genomic coordinates: already formatted
    if ':' in mutation and ('chr' in mutation.lower() or mutation.startswith('g.')):
        # Try to parse chr:pos:ref:alt format
        return mutation

    # p. notation - protein level, need to map back to DNA
    if mutation.startswith('p.'):
        return f"PROTEIN_LEVEL:{mutation}"

    # Unknown format
    return None


def anchor_mutation(
    pmid: str,
    mutation: str,
    pubtator: dict,
    articles: dict,
    gnormplus: dict
) -> AnchorResult:
    """
    Anchor a mutation to a gene using tiered approach.
    Species-agnostic: anchors to closest gene regardless of taxon,
    includes ncbi_taxon_id in output for downstream filtering.

    Tier 1: Full-text NER proximity (GNormPlus)
    Tier 2: Abstract proximity (PubTator offsets)
    Tier 3: Most-mentioned gene fallback
    """
    pt = pubtator.get(pmid)

    # No PubTator data
    if not pt:
        return AnchorResult(None, None, "no_pubtator", "none")

    gene_offsets_pubtator = get_gene_offsets(pt)

    # No genes annotated
    if not gene_offsets_pubtator:
        return AnchorResult(None, None, "no_genes", "none")

    # Single gene - direct assignment (highest confidence)
    if len(gene_offsets_pubtator) == 1:
        gene_id = list(gene_offsets_pubtator.keys())[0]
        gene_info = gene_offsets_pubtator[gene_id]
        gene_name = gene_info["name"]
        taxon_id = gene_info.get("taxon_id")
        vep_input = format_vep_input(mutation, gene_id)
        return AnchorResult(
            gene_id,
            gene_name,
            "single_gene",
            "high",
            vep_input_string=vep_input,
            ncbi_taxon_id=taxon_id
        )

    # Multi-gene: Try tiered proximity matching
    article_text = articles.get(pmid, "")
    vep_input = format_vep_input(mutation)

    # === TIER 1: Full-text NER (GNormPlus) ===
    if gnormplus and pmid in gnormplus:
        gnorm_genes = gnormplus[pmid]

        if article_text:
            # Use PubTator-formatted text for offset alignment
            pubtator_text = get_pubtator_text(article_text)
            mutation_pos = find_mutation_in_text(pubtator_text, mutation)

            if mutation_pos is not None:
                gene_id, gene_name, dist, taxon_id = find_closest_gene(mutation_pos, gnorm_genes)

                if gene_id:
                    # Confidence based on distance
                    if dist <= 200:
                        confidence = "high"
                    elif dist <= 500:
                        confidence = "medium"
                    else:
                        confidence = "low"

                    return AnchorResult(
                        gene_id, gene_name, "tier1_fulltext",
                        confidence, dist, vep_input, taxon_id
                    )

    # === TIER 2: Abstract proximity (PubTator offsets) ===
    if article_text:
        abstract = extract_abstract(article_text)
        mutation_pos = find_mutation_in_text(abstract, mutation)

        if mutation_pos is not None:
            gene_id, gene_name, dist, taxon_id = find_closest_gene(
                mutation_pos, gene_offsets_pubtator
            )

            if gene_id:
                # Confidence based on distance
                if dist <= 200:
                    confidence = "high"
                elif dist <= 500:
                    confidence = "medium"
                else:
                    confidence = "low"

                return AnchorResult(
                    gene_id, gene_name, "tier2_abstract",
                    confidence, dist, vep_input, taxon_id
                )

    # === TIER 3: Most-mentioned fallback ===
    gene_id, gene_name, taxon_id = get_most_mentioned_gene(pt)
    if gene_id:
        return AnchorResult(
            gene_id, gene_name, "tier3_most_mentioned",
            "low", vep_input_string=vep_input, ncbi_taxon_id=taxon_id
        )

    return AnchorResult(None, None, "anchor_failed", "none")


def main():
    global PUBTATOR_FILE, QA_FILE, ARTICLES_FILE, GNORMPLUS_DIR, OUTPUT_FILE, STATS_FILE
    parser = argparse.ArgumentParser(description="Tiered gene anchoring (v3) of Phase 4 QA records")
    parser.add_argument("--input", type=Path, default=QA_FILE,
                        help="m2tqa_gene_recovered.jsonl (run_gene_recovery.py output)")
    parser.add_argument("--pubtator", type=Path, default=PUBTATOR_FILE,
                        help="pubtator_annotations.jsonl (fetch_pubtator_bioc.py)")
    parser.add_argument("--ner-dir", type=Path, default=GNORMPLUS_DIR,
                        help="Directory of Tier 1 NER batch_*.out files (fetch_pubtator_batches.py or GNormPlus)")
    parser.add_argument("--articles", type=Path, default=ARTICLES_FILE,
                        help='Source article corpus JSON {pmid: {"article": ...}}')
    parser.add_argument("--output", type=Path, default=OUTPUT_FILE)
    parser.add_argument("--stats-output", type=Path, default=STATS_FILE)
    args = parser.parse_args()
    QA_FILE, PUBTATOR_FILE, GNORMPLUS_DIR = args.input, args.pubtator, args.ner_dir
    ARTICLES_FILE, OUTPUT_FILE, STATS_FILE = args.articles, args.output, args.stats_output

    print("=" * 60)
    print("Tiered Gene Anchoring Pipeline v3")
    print("=" * 60)
    print()

    # Load PubTator annotations
    print("Loading PubTator annotations...")
    pubtator = load_pubtator()
    print(f"Loaded {len(pubtator):,} PMID annotations")

    # Load GNormPlus annotations (directory-based)
    print()
    gnormplus = load_gnormplus_annotations(GNORMPLUS_DIR)

    # Get target PMIDs from QA file
    print("\nScanning QA file for target PMIDs...")
    target_pmids = set()
    with open(QA_FILE) as f:
        for line in f:
            qa = json.loads(line)
            pmid = str(qa.get("pmid", ""))
            if pmid:
                target_pmids.add(pmid)
    print(f"Target PMIDs: {len(target_pmids):,}")

    # Load articles for target PMIDs only (memory efficient)
    articles = build_article_index(target_pmids)

    # Stats
    stats = defaultdict(int)

    print("\nProcessing QA pairs...")
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    with open(QA_FILE) as f_in, open(OUTPUT_FILE, "w") as f_out:
        for line in tqdm(f_in, desc="Anchoring"):
            qa = json.loads(line)
            stats["total"] += 1

            pmid = str(qa.get("pmid", ""))
            mutation = qa.get("mutation", "")

            result = anchor_mutation(pmid, mutation, pubtator, articles, gnormplus)

            # Update stats
            stats[f"tier_{result.tier}"] += 1
            stats[f"conf_{result.confidence}"] += 1

            # Add results to record
            qa["anchor_gene_id"] = result.gene_id
            qa["anchor_gene_name"] = result.gene_name
            qa["anchor_tier"] = result.tier
            qa["anchor_confidence"] = result.confidence
            qa["ncbi_taxon_id"] = result.ncbi_taxon_id
            if result.distance is not None:
                qa["anchor_distance"] = result.distance
            if result.vep_input_string:
                qa["vep_input_string"] = result.vep_input_string

            # Track taxon stats
            if result.ncbi_taxon_id:
                stats[f"taxon_{result.ncbi_taxon_id}"] += 1
            else:
                stats["taxon_unknown"] += 1

            f_out.write(json.dumps(qa) + "\n")

    # Save stats
    with open(STATS_FILE, "w") as f:
        json.dump(dict(stats), f, indent=2)

    # Print summary
    total = stats["total"]
    print()
    print("=" * 60)
    print("ANCHORING RESULTS")
    print("=" * 60)
    print(f"Total QA pairs: {total:,}")
    print()

    print("By Tier:")
    tier_order = ["single_gene", "tier1_fulltext", "tier2_abstract",
                  "tier3_most_mentioned", "no_pubtator",
                  "no_genes", "anchor_failed"]
    for tier in tier_order:
        key = f"tier_{tier}"
        if key in stats:
            count = stats[key]
            print(f"  {tier}: {count:,} ({count/total*100:.1f}%)")

    print()
    print("By Confidence:")
    for conf in ["high", "medium", "low", "none"]:
        key = f"conf_{conf}"
        if key in stats:
            count = stats[key]
            print(f"  {conf}: {count:,} ({count/total*100:.1f}%)")

    print()
    print("By Species (Top 10):")
    taxon_stats = {k: v for k, v in stats.items() if k.startswith("taxon_")}
    sorted_taxa = sorted(taxon_stats.items(), key=lambda x: x[1], reverse=True)[:10]
    for taxon_key, count in sorted_taxa:
        taxon = taxon_key.replace("taxon_", "")
        label = "Human" if taxon == "9606" else taxon
        print(f"  {label}: {count:,} ({count/total*100:.1f}%)")

    # Calculate anchored total
    anchored = sum(v for k, v in stats.items()
                   if k.startswith("tier_") and k not in
                   ["tier_no_pubtator", "tier_no_genes", "tier_anchor_failed"])
    print()
    print(f"Successfully anchored: {anchored:,} ({anchored/total*100:.1f}%)")
    print(f"Output: {OUTPUT_FILE}")
    print(f"Stats: {STATS_FILE}")


if __name__ == "__main__":
    main()
