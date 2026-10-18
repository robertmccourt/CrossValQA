#!/usr/bin/env python3
"""
Fetch full-text articles from PMC for pipeline expansion.

For each target gene, searches PMC for open-access articles discussing
mutations/variants, then fetches full text via E-utilities.

Input:  a targets JSONL ({"gene": ...} per line; default expansion/output/expansion_targets.jsonl)
Output: an articles JSON {pmid: {article, gene, pmc_id}} (default expansion/output/expansion_articles.json)
        + <output>_stats.json
NCBI E-utilities (network). Set NCBI_API_KEY (via config.get_ncbi_api_key) for higher
rate limits; NCBI_EMAIL (env or repo-root .env) is required as the NCBI contact email.

Usage:
    python fetch_pmc_articles.py --limit 10  # pilot with 10 genes
    python fetch_pmc_articles.py                                            # Batch 1
    python fetch_pmc_articles.py --targets output/batch2_targets.jsonl --output output/batch2_articles.json
    python fetch_pmc_articles.py --targets output/batch2_failed_genes.jsonl \
        --output output/batch2_round2_articles.json --papers-per-gene 5 --exclude-from PREVIOUS.json
"""

import json
import os
import time
import argparse
import requests
import xml.etree.ElementTree as ET
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import PIPELINE_ROOT, get_ncbi_api_key

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
DEFAULT_TARGETS = OUTPUT_DIR / "expansion_targets.jsonl"
DEFAULT_ARTICLES = OUTPUT_DIR / "expansion_articles.json"

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ELINK_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"

REQUESTS_PER_SECOND = 8
_NCBI_EMAIL = ""


def get_ncbi_email() -> str:
    """Get the NCBI E-utilities contact email from NCBI_EMAIL (env or repo-root .env)."""
    global _NCBI_EMAIL
    if _NCBI_EMAIL:
        return _NCBI_EMAIL
    email = os.environ.get("NCBI_EMAIL", "")
    if not email:
        from dotenv import load_dotenv
        load_dotenv(PIPELINE_ROOT / ".env")
        email = os.environ.get("NCBI_EMAIL", "")
    if not email:
        raise ValueError(
            "NCBI_EMAIL not set. NCBI E-utilities requires a contact email; "
            "add NCBI_EMAIL=you@example.org to .env or export it."
        )
    _NCBI_EMAIL = email
    return email


def load_targets(path: Path, limit: int = 0) -> list[dict]:
    """Load expansion targets."""
    targets = []
    with open(path) as f:
        for line in f:
            targets.append(json.loads(line))
    if limit > 0:
        targets = targets[:limit]
    return targets


TIERED_QUERIES = [
    # Tier 1: HGVS notation — strongest signal that the article discusses specific variants
    '"{gene}"[tiab] AND (c\\.*[tiab] OR p\\.*[tiab] OR missense[tiab] OR nonsense[tiab] OR frameshift[tiab])',
    # Tier 2: Clinical/pathogenic framing
    '"{gene}"[tiab] AND (pathogenic[tiab] OR "clinical significance"[tiab] OR "disease-causing"[tiab])',
    # Tier 3: Broad fallback (original query)
    '"{gene}" AND (mutation OR variant OR pathogenic)',
]


def _esearch(query: str, api_key: str, retmax: int) -> list[str]:
    params = {
        "db": "pmc",
        "term": query,
        "retmax": retmax,
        "retmode": "json",
        "sort": "relevance",
        "tool": "crossvalqa",
    }
    params["email"] = get_ncbi_email()
    if api_key:
        params["api_key"] = api_key
    try:
        resp = requests.get(ESEARCH_URL, params=params, timeout=15)
        if resp.status_code != 200:
            return []
        return resp.json().get("esearchresult", {}).get("idlist", [])
    except Exception:
        return []


def search_pmc_for_gene(gene: str, api_key: str, retmax: int = 1, exclude_pmc_ids: set | None = None) -> list[str]:
    """Search PMC using tiered queries. Returns IDs from the first tier with results."""
    for tpl in TIERED_QUERIES:
        query = tpl.format(gene=gene)
        ids = _esearch(query, api_key, retmax * 3 if exclude_pmc_ids else retmax)
        if exclude_pmc_ids:
            ids = [i for i in ids if f"PMC{i}" not in exclude_pmc_ids]
        if ids:
            return ids[:retmax]
    return []


def pmc_id_to_pmid(pmc_id: str, api_key: str) -> str:
    """Convert a PMC ID to a PMID via elink."""
    params = {
        "dbfrom": "pmc",
        "db": "pubmed",
        "id": pmc_id,
        "retmode": "json",
        "tool": "crossvalqa",
    }
    params["email"] = get_ncbi_email()
    if api_key:
        params["api_key"] = api_key

    try:
        resp = requests.get(ELINK_URL, params=params, timeout=15)
        if resp.status_code != 200:
            return ""
        data = resp.json()
        linksets = data.get("linksets", [])
        if linksets:
            links = linksets[0].get("linksetdbs", [])
            if links:
                link_ids = links[0].get("links", [])
                if link_ids:
                    return link_ids[0]
        return ""
    except Exception:
        return ""


def extract_text_from_jats_xml(xml_str: str) -> str | None:
    """Extract readable text from JATS/PMC XML."""
    try:
        root = ET.fromstring(xml_str)
    except ET.ParseError:
        return None

    body = root.find('.//body')
    if body is None:
        return None

    def get_text(elem) -> str:
        parts = []
        if elem.text:
            parts.append(elem.text)
        for child in elem:
            if child.tag in ('table-wrap', 'fig', 'disp-formula', 'inline-formula',
                             'supplementary-material'):
                if child.tail:
                    parts.append(child.tail)
                continue
            # Keep xref text (citations, references)
            parts.append(get_text(child))
            if child.tail:
                parts.append(child.tail)
        return ''.join(parts)

    title_elem = root.find('.//article-title')
    title = get_text(title_elem) if title_elem is not None else ""

    abstract_elem = root.find('.//abstract')
    abstract = get_text(abstract_elem) if abstract_elem is not None else ""

    body_text = get_text(body)

    parts = []
    if title:
        parts.append(title)
    if abstract:
        parts.append(abstract)
    if body_text:
        parts.append(body_text)

    full_text = "\n\n".join(parts)
    return full_text if len(full_text) > 500 else None


def fetch_pmc_fulltext(pmc_id: str, api_key: str) -> str | None:
    """Fetch full text from PMC via efetch. Returns plain text or None."""
    params = {
        "db": "pmc",
        "id": pmc_id,
        "rettype": "xml",
        "tool": "crossvalqa",
    }
    params["email"] = get_ncbi_email()
    if api_key:
        params["api_key"] = api_key

    try:
        resp = requests.get(EFETCH_URL, params=params, timeout=30)
        if resp.status_code != 200:
            return None
        return extract_text_from_jats_xml(resp.text)
    except Exception as e:
        print(f"  Error fetching PMC{pmc_id}: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Fetch PMC full-text articles for expansion")
    parser.add_argument("--limit", type=int, default=0, help="Limit genes (0=all)")
    parser.add_argument("--resume", action="store_true", help="Skip already-fetched genes")
    parser.add_argument("--targets", type=str, default=str(DEFAULT_TARGETS),
                        help="Path to targets JSONL file")
    parser.add_argument("--output", type=str, default=str(DEFAULT_ARTICLES),
                        help="Path to output articles JSON file")
    parser.add_argument("--papers-per-gene", type=int, default=1,
                        help="Number of papers to fetch per gene (default: 1)")
    parser.add_argument("--exclude-from", type=str, default=None,
                        help="Existing articles JSON to exclude (avoid re-fetching same papers)")
    args = parser.parse_args()

    targets_file = Path(args.targets)
    articles_file = Path(args.output)

    get_ncbi_email()  # fail fast if NCBI_EMAIL is missing
    api_key = get_ncbi_api_key()
    if api_key:
        print("Using NCBI API key")
    else:
        print("Warning: No NCBI API key. Rate limits will be lower.")

    targets = load_targets(targets_file, limit=args.limit)
    print(f"Target genes: {len(targets)}")
    print(f"Papers per gene: {args.papers_per_gene}")

    # Load existing articles if resuming
    articles = {}
    if args.resume and articles_file.exists():
        with open(articles_file) as f:
            articles = json.load(f)
        print(f"Resuming: {len(articles)} already fetched")

    # Build set of PMC IDs to exclude (from previous fetches)
    exclude_pmc_ids = set()
    if args.exclude_from:
        exclude_path = Path(args.exclude_from)
        if exclude_path.exists():
            with open(exclude_path) as f:
                prev_data = json.load(f)
            # Support both formats: list of PMC IDs or articles dict
            if isinstance(prev_data, list):
                exclude_pmc_ids.update(prev_data)
            else:
                for info in prev_data.values():
                    if "pmc_id" in info:
                        exclude_pmc_ids.add(info["pmc_id"])
            print(f"Excluding {len(exclude_pmc_ids)} previously fetched PMC IDs")

    # Also exclude PMC IDs already in current output
    for info in articles.values():
        if "pmc_id" in info:
            exclude_pmc_ids.add(info["pmc_id"])

    fetched_genes = set()
    if articles:
        for pmid, info in articles.items():
            if "gene" in info:
                fetched_genes.add(info["gene"])

    searched = 0
    pmc_found = 0
    fulltext_found = 0
    no_results = 0
    no_fulltext = 0

    for idx, target in enumerate(targets):
        gene = target["gene"]

        if gene in fetched_genes:
            continue

        # Step 1: Search PMC for this gene (multiple results)
        pmc_ids = search_pmc_for_gene(gene, api_key, retmax=args.papers_per_gene,
                                       exclude_pmc_ids=exclude_pmc_ids)
        searched += 1
        time.sleep(1.0 / REQUESTS_PER_SECOND)

        if not pmc_ids:
            no_results += 1
            continue

        pmc_found += 1
        gene_got_article = False

        for pmc_id in pmc_ids:
            # Step 2: Get PMID
            pmid = pmc_id_to_pmid(pmc_id, api_key)
            time.sleep(1.0 / REQUESTS_PER_SECOND)

            if not pmid:
                pmid = f"PMC{pmc_id}"

            # Skip if this PMID is already in output
            if pmid in articles:
                continue

            # Step 3: Fetch full text
            text = fetch_pmc_fulltext(pmc_id, api_key)
            time.sleep(1.0 / REQUESTS_PER_SECOND)

            if text:
                articles[pmid] = {
                    "article": text,
                    "gene": gene,
                    "pmc_id": f"PMC{pmc_id}",
                }
                exclude_pmc_ids.add(f"PMC{pmc_id}")
                fulltext_found += 1
                gene_got_article = True
            else:
                no_fulltext += 1

        if gene_got_article:
            fetched_genes.add(gene)

        # Progress + periodic save
        if (searched) % 25 == 0:
            print(f"  [{searched}/{len(targets)}] Searched: {searched}, PMC found: {pmc_found}, "
                  f"Full text: {fulltext_found}, No results: {no_results}, No full text: {no_fulltext}")

            articles_file.parent.mkdir(parents=True, exist_ok=True)
            with open(articles_file, 'w') as f:
                json.dump(articles, f)

    # Final save
    articles_file.parent.mkdir(parents=True, exist_ok=True)
    with open(articles_file, 'w') as f:
        json.dump(articles, f)

    print(f"\n=== RESULTS ===")
    print(f"Genes searched: {searched}")
    print(f"PMC articles found: {pmc_found} ({100 * pmc_found / max(searched, 1):.1f}%)")
    print(f"Full text fetched: {fulltext_found} ({100 * fulltext_found / max(searched, 1):.1f}%)")
    print(f"No PMC results: {no_results}")
    print(f"No full text (publisher restricted): {no_fulltext}")
    print(f"\nOutput: {articles_file}")
    print(f"Total articles: {len(articles)}")

    # Save stats
    stats = {
        "genes_searched": searched,
        "pmc_found": pmc_found,
        "fulltext_fetched": fulltext_found,
        "no_results": no_results,
        "no_fulltext": no_fulltext,
        "total_articles": len(articles),
    }
    stats_file = articles_file.parent / (articles_file.stem + "_stats.json")
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)


if __name__ == "__main__":
    main()
