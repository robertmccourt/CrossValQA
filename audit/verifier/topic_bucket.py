"""Question-type bucketer for the MutQA deterministic auditor. The dataset's `question_type` field is already a single
canonical level ('mechanism', 'function', 'pathogenicity', ...), so this is
near-identity. The wrapper exists so verifier code can call `bucket_topic(qt)`
without caring whether the source was clean.
"""
from __future__ import annotations

_CANONICAL = {
    "mechanism": "mechanism",
    "function": "function",
    "pathogenicity": "pathogenicity",
    "structure": "structure",
    "regulation": "regulation",
    "interaction": "interaction",
    "expression": "expression",
    "localization": "localization",
}


def bucket_topic(qt: str | None) -> str:
    if not qt:
        return "unknown"
    s = qt.strip().lower()
    return _CANONICAL.get(s, s)
