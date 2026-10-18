#!/usr/bin/env python3
"""
Tests for get_mutation_dna.py.

The parser tests are pure logic. The TestGenome* integration tests need the
GRCh38 + UniProt CDS reference files (see get_mutation_dna.py /
CROSSVALQA_REF_DIR) and are skipped when those files are missing.

Run:
    python -m pytest release/test_get_mutation_dna.py
    python test_get_mutation_dna.py                          # all tests
    python test_get_mutation_dna.py TestParseSubstitutions   # one class
    python test_get_mutation_dna.py -v                       # verbose
"""

import unittest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from get_mutation_dna import (
    parse_dna_mutation,
    get_mutation_dna_window,
    _complement,
    _revcomp,
    GENOME_FA, UNIPROT_SEQ, GTF_SWISS, GTF_TREMBL,
)

HAVE_REF_DATA = all(p.exists() for p in (GENOME_FA, UNIPROT_SEQ, GTF_SWISS, GTF_TREMBL))
requires_ref_data = unittest.skipUnless(
    HAVE_REF_DATA, "GRCh38/UniProt reference data not found (set CROSSVALQA_REF_DIR)")


# ═══════════════════════════════════════════════════════════════════════════════
#  Helper tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestHelpers(unittest.TestCase):

    def test_complement(self):
        self.assertEqual(_complement("ACGT"), "TGCA")

    def test_complement_lowercase(self):
        self.assertEqual(_complement("acgt"), "tgca")

    def test_revcomp_palindrome(self):
        self.assertEqual(_revcomp("ACGT"), "ACGT")

    def test_revcomp(self):
        self.assertEqual(_revcomp("AACG"), "CGTT")

    def test_complement_roundtrip(self):
        seq = "ATCGATCGATCG"
        self.assertEqual(_complement(_complement(seq)), seq)

    def test_revcomp_roundtrip(self):
        seq = "ATCGATCG"
        self.assertEqual(_revcomp(_revcomp(seq)), seq)


# ═══════════════════════════════════════════════════════════════════════════════
#  Parser tests — pure logic, no genome needed
# ═══════════════════════════════════════════════════════════════════════════════

class TestParseSubstitutions(unittest.TestCase):

    def test_standard_c_dot(self):
        r = parse_dna_mutation("c.1799T>A")
        self.assertEqual(r["type"], "sub")
        self.assertEqual(r["pos"], 1799)
        self.assertEqual(r["ref_bases"], "T")
        self.assertEqual(r["alt_bases"], "A")
        self.assertIsNone(r["offset"])

    def test_no_c_dot_prefix(self):
        r = parse_dna_mutation("1799T>A")
        self.assertEqual(r["type"], "sub")
        self.assertEqual(r["pos"], 1799)

    def test_bare_ref_pos_alt(self):
        r = parse_dna_mutation("T1799A")
        self.assertEqual(r["type"], "sub")
        self.assertEqual(r["pos"], 1799)
        self.assertEqual(r["ref_bases"], "T")
        self.assertEqual(r["alt_bases"], "A")

    def test_old_style_with_c_prefix(self):
        r = parse_dna_mutation("c.C538T")
        self.assertEqual(r["pos"], 538)
        self.assertEqual(r["ref_bases"], "C")
        self.assertEqual(r["alt_bases"], "T")

    def test_lowercase_input(self):
        r = parse_dna_mutation("c.100a>g")
        self.assertEqual(r["ref_bases"], "A")
        self.assertEqual(r["alt_bases"], "G")

    def test_intronic_positive_offset(self):
        r = parse_dna_mutation("c.316+1G>A")
        self.assertEqual(r["pos"], 316)
        self.assertEqual(r["offset"], 1)
        self.assertEqual(r["ref_bases"], "G")

    def test_intronic_negative_offset(self):
        r = parse_dna_mutation("c.754-12C>G")
        self.assertEqual(r["pos"], 754)
        self.assertEqual(r["offset"], -12)

    def test_arrow_notation(self):
        r = parse_dna_mutation("c.1178T\u2192C")
        self.assertEqual(r["type"], "sub")
        self.assertEqual(r["alt_bases"], "C")

    def test_single_digit_position(self):
        r = parse_dna_mutation("c.1A>T")
        self.assertEqual(r["pos"], 1)

    def test_large_position(self):
        r = parse_dna_mutation("c.99999G>C")
        self.assertEqual(r["pos"], 99999)

    def test_all_base_pairs(self):
        for ref in "ACGT":
            for alt in "ACGT":
                if ref == alt:
                    continue
                r = parse_dna_mutation(f"c.100{ref}>{alt}")
                self.assertEqual(r["ref_bases"], ref)
                self.assertEqual(r["alt_bases"], alt)


class TestParseDeletions(unittest.TestCase):

    def test_single_with_base(self):
        r = parse_dna_mutation("c.210delA")
        self.assertEqual(r["type"], "del")
        self.assertEqual(r["pos"], 210)
        self.assertIsNone(r["end_pos"])
        self.assertEqual(r["ref_bases"], "A")

    def test_single_no_base(self):
        r = parse_dna_mutation("c.210del")
        self.assertEqual(r["type"], "del")
        self.assertIsNone(r["ref_bases"])

    def test_range_with_bases(self):
        r = parse_dna_mutation("c.241_242delCT")
        self.assertEqual(r["type"], "del")
        self.assertEqual(r["pos"], 241)
        self.assertEqual(r["end_pos"], 242)
        self.assertEqual(r["ref_bases"], "CT")

    def test_range_no_bases(self):
        r = parse_dna_mutation("c.1355_1356del")
        self.assertEqual(r["end_pos"], 1356)
        self.assertIsNone(r["ref_bases"])

    def test_three_base_deletion(self):
        r = parse_dna_mutation("c.1521_1523delCTT")
        self.assertEqual(r["ref_bases"], "CTT")

    def test_intronic_single_del(self):
        r = parse_dna_mutation("c.226+2delT")
        self.assertEqual(r["offset"], 2)
        self.assertEqual(r["ref_bases"], "T")


class TestParseDuplications(unittest.TestCase):

    def test_single_with_base(self):
        r = parse_dna_mutation("c.2527dupG")
        self.assertEqual(r["type"], "dup")
        self.assertEqual(r["pos"], 2527)
        self.assertEqual(r["ref_bases"], "G")

    def test_single_no_base(self):
        r = parse_dna_mutation("c.2527dup")
        self.assertIsNone(r["ref_bases"])

    def test_range_with_bases(self):
        r = parse_dna_mutation("c.123_125dupACG")
        self.assertEqual(r["end_pos"], 125)
        self.assertEqual(r["ref_bases"], "ACG")

    def test_range_no_bases(self):
        r = parse_dna_mutation("c.123_125dup")
        self.assertIsNone(r["ref_bases"])

    def test_intronic_dup(self):
        r = parse_dna_mutation("c.1012+2dupT")
        self.assertEqual(r["offset"], 2)


class TestParseInsertions(unittest.TestCase):

    def test_standard_range(self):
        r = parse_dna_mutation("c.123_124insACG")
        self.assertEqual(r["type"], "ins")
        self.assertEqual(r["pos"], 123)
        self.assertEqual(r["end_pos"], 124)
        self.assertEqual(r["alt_bases"], "ACG")

    def test_single_base_insert(self):
        r = parse_dna_mutation("c.100_101insT")
        self.assertEqual(r["alt_bases"], "T")

    def test_nonstandard_single_pos(self):
        r = parse_dna_mutation("c.640insC")
        self.assertEqual(r["type"], "ins")
        self.assertEqual(r["pos"], 640)
        self.assertEqual(r["end_pos"], 641)

    def test_multi_base_insert(self):
        r = parse_dna_mutation("c.512_513insAACG")
        self.assertEqual(r["alt_bases"], "AACG")


class TestParseDelins(unittest.TestCase):

    def test_single_pos(self):
        r = parse_dna_mutation("c.123delinsACG")
        self.assertEqual(r["type"], "delins")
        self.assertEqual(r["pos"], 123)
        self.assertIsNone(r["end_pos"])
        self.assertEqual(r["alt_bases"], "ACG")

    def test_range(self):
        r = parse_dna_mutation("c.123_125delinsACGT")
        self.assertEqual(r["end_pos"], 125)
        self.assertEqual(r["alt_bases"], "ACGT")


class TestParseErrors(unittest.TestCase):

    def test_empty_string(self):
        with self.assertRaises(ValueError):
            parse_dna_mutation("")

    def test_garbage(self):
        with self.assertRaises(ValueError):
            parse_dna_mutation("not_a_mutation")

    def test_protein_notation(self):
        with self.assertRaises(ValueError):
            parse_dna_mutation("p.R248W")

    def test_just_c_dot(self):
        with self.assertRaises(ValueError):
            parse_dna_mutation("c.")

    def test_bare_intronic_pos_no_change(self):
        """c.1357+1 — no base change specified, should fail."""
        with self.assertRaises(ValueError):
            parse_dna_mutation("c.1357+1")

    def test_raw_field_preserved(self):
        r = parse_dna_mutation("  c.100A>G  ")
        self.assertEqual(r["raw"], "c.100A>G")


# ═══════════════════════════════════════════════════════════════════════════════
#  Integration tests — require genome + annotation files on disk
# ═══════════════════════════════════════════════════════════════════════════════

@requires_ref_data
class TestGenomeSubstitutions(unittest.TestCase):

    def test_braf_v600e(self):
        """BRAF c.1799T>A (V600E), chr7, minus strand."""
        r = get_mutation_dna_window("BRAF", "c.1799T>A", window=20)
        self.assertEqual(r["gene"], "BRAF")
        self.assertEqual(r["mut_type"], "sub")
        self.assertEqual(r["strand"], "-")
        self.assertTrue(r["ref_matches_genome"])
        self.assertEqual(len(r["ref_dna_coding"]), 20)
        self.assertEqual(len(r["mut_dna_coding"]), 20)
        # exactly one base differs
        diffs = sum(a != b for a, b in zip(r["ref_dna_coding"], r["mut_dna_coding"]))
        self.assertEqual(diffs, 1)

    def test_kras_g12d_bare(self):
        """KRAS 35G>A (G12D), minus strand, bare format."""
        r = get_mutation_dna_window("KRAS", "35G>A", window=20)
        self.assertTrue(r["ref_matches_genome"])
        self.assertEqual(r["strand"], "-")

    def test_bare_formats_equivalent(self):
        """35G>A and G35A should produce identical output."""
        r1 = get_mutation_dna_window("KRAS", "35G>A", window=20)
        r2 = get_mutation_dna_window("KRAS", "G35A", window=20)
        self.assertEqual(r1["ref_dna_coding"], r2["ref_dna_coding"])
        self.assertEqual(r1["mut_dna_coding"], r2["mut_dna_coding"])
        self.assertEqual(r1["genomic_positions"], r2["genomic_positions"])

    def test_plus_strand_coding_equals_genomic(self):
        """BRCA2 is on + strand — coding should equal genomic."""
        r = get_mutation_dna_window("BRCA2", "c.7397T>G", window=20)
        self.assertEqual(r["strand"], "+")
        self.assertEqual(r["ref_dna_coding"], r["ref_dna_genomic"])
        self.assertEqual(r["mut_dna_coding"], r["mut_dna_genomic"])

    def test_minus_strand_coding_is_revcomp(self):
        """For - strand, coding = revcomp(genomic)."""
        r = get_mutation_dna_window("TP53", "c.818G>A", window=30)
        self.assertEqual(r["strand"], "-")
        self.assertEqual(r["ref_dna_coding"], _revcomp(r["ref_dna_genomic"]))
        self.assertEqual(r["mut_dna_coding"], _revcomp(r["mut_dna_genomic"]))


@requires_ref_data
class TestGenomeMismatch(unittest.TestCase):

    def test_correct_ref(self):
        r = get_mutation_dna_window("BRAF", "c.1799T>A", window=20)
        self.assertTrue(r["ref_matches_genome"])
        self.assertEqual(len(r["warnings"]), 0)

    def test_wrong_ref_warns(self):
        """c.1799C>A — BRAF pos 1799 is T not C."""
        r = get_mutation_dna_window("BRAF", "c.1799C>A", window=20)
        self.assertFalse(r["ref_matches_genome"])
        self.assertEqual(len(r["warnings"]), 1)
        self.assertIn("mismatch", r["warnings"][0].lower())

    def test_mismatch_still_returns_window(self):
        r = get_mutation_dna_window("BRAF", "c.1799C>A", window=20)
        self.assertEqual(len(r["ref_dna_coding"]), 20)


@requires_ref_data
class TestGenomeIntronic(unittest.TestCase):

    def test_splice_donor(self):
        """TP53 c.672+1G>A — splice donor."""
        r = get_mutation_dna_window("TP53", "c.672+1G>A", window=20)
        self.assertTrue(r["ref_matches_genome"])
        self.assertEqual(r["intronic_offset"], 1)

    def test_intronic_pos_differs_from_exonic(self):
        """c.672+1 and c.672 should map to different genomic positions."""
        r_intr = get_mutation_dna_window("TP53", "c.672+1G>A", window=20)
        r_exon = get_mutation_dna_window("TP53", "c.672G>T", window=20)
        self.assertNotEqual(r_intr["genomic_positions"], r_exon["genomic_positions"])

    def test_intronic_negative_offset(self):
        r = get_mutation_dna_window("TP53", "c.673-2A>T", window=20)
        self.assertEqual(r["intronic_offset"], -2)


@requires_ref_data
class TestGenomeDeletions(unittest.TestCase):

    def test_single_del_shrinks(self):
        """c.5946delT — window shrinks by 1."""
        r = get_mutation_dna_window("BRCA2", "c.5946delT", window=40)
        self.assertEqual(r["mut_type"], "del")
        self.assertTrue(r["ref_matches_genome"])
        self.assertEqual(len(r["ref_dna_coding"]), 40)
        self.assertEqual(len(r["mut_dna_coding"]), 39)

    def test_range_del_cftr(self):
        """c.1521_1523delCTT — CFTR delta-F508, 3 bases deleted."""
        r = get_mutation_dna_window("CFTR", "c.1521_1523delCTT", window=40)
        self.assertTrue(r["ref_matches_genome"])
        self.assertEqual(len(r["mut_dna_coding"]), 37)

    def test_del_no_ref_bases(self):
        """c.1521_1523del — no ref bases specified, should still work."""
        r = get_mutation_dna_window("CFTR", "c.1521_1523del", window=40)
        self.assertEqual(len(r["mut_dna_coding"]), 37)

    def test_del_wrong_ref_warns(self):
        """c.5946delG — BRCA2 pos 5946 is T not G."""
        r = get_mutation_dna_window("BRCA2", "c.5946delG", window=20)
        self.assertFalse(r["ref_matches_genome"])

    def test_del_minus_strand(self):
        r = get_mutation_dna_window("TP53", "c.818del", window=30)
        self.assertEqual(r["strand"], "-")
        self.assertEqual(len(r["mut_dna_coding"]), 29)
        self.assertEqual(r["ref_dna_coding"], _revcomp(r["ref_dna_genomic"]))

    def test_del_is_subsequence(self):
        """Mutated sequence should be ref with the deleted base(s) removed."""
        r = get_mutation_dna_window("BRCA2", "c.5946delT", window=40)
        ref, mut = r["ref_dna_coding"], r["mut_dna_coding"]
        # Find which index was removed
        found = any(ref[:i] + ref[i+1:] == mut for i in range(len(ref)))
        self.assertTrue(found, "Mut is not ref with one base removed")


@requires_ref_data
class TestGenomeDuplications(unittest.TestCase):

    def test_single_dup_grows(self):
        """c.5266dupC — BRCA1, minus strand, grows by 1."""
        r = get_mutation_dna_window("BRCA1", "c.5266dupC", window=40)
        self.assertEqual(r["mut_type"], "dup")
        self.assertTrue(r["ref_matches_genome"])
        self.assertEqual(len(r["mut_dna_coding"]), 41)

    def test_dup_with_and_without_base_same_result(self):
        """c.5266dupC and c.5266dup should produce identical mutated windows."""
        r1 = get_mutation_dna_window("BRCA1", "c.5266dupC", window=40)
        r2 = get_mutation_dna_window("BRCA1", "c.5266dup", window=40)
        self.assertEqual(r1["mut_dna_coding"], r2["mut_dna_coding"])

    def test_dup_plus_strand(self):
        r = get_mutation_dna_window("BRCA2", "c.5946dupT", window=30)
        self.assertEqual(r["strand"], "+")
        self.assertEqual(len(r["mut_dna_coding"]), 31)
        self.assertTrue(r["ref_matches_genome"])


@requires_ref_data
class TestGenomeInsertions(unittest.TestCase):

    def test_single_base_insert(self):
        """c.5946_5947insA — BRCA2, + strand, grows by 1."""
        r = get_mutation_dna_window("BRCA2", "c.5946_5947insA", window=30)
        self.assertEqual(r["mut_type"], "ins")
        self.assertEqual(len(r["mut_dna_coding"]), 31)

    def test_multi_base_insert(self):
        """Insert 3 bases — grows by 3."""
        r = get_mutation_dna_window("BRCA2", "c.5946_5947insACG", window=30)
        self.assertEqual(len(r["mut_dna_coding"]), 33)

    def test_nonstandard_single_pos(self):
        r = get_mutation_dna_window("BRCA2", "c.5946insA", window=30)
        self.assertEqual(r["mut_type"], "ins")
        self.assertEqual(len(r["mut_dna_coding"]), 31)

    def test_insert_minus_strand(self):
        r = get_mutation_dna_window("TP53", "c.817_818insA", window=20)
        self.assertEqual(r["strand"], "-")
        self.assertEqual(len(r["mut_dna_coding"]), 21)
        self.assertEqual(len(r["mut_dna_genomic"]), 21)


@requires_ref_data
class TestGenomeDelins(unittest.TestCase):

    def test_del1_ins3(self):
        """c.100delinsACG — delete 1, insert 3 → net +2."""
        r = get_mutation_dna_window("BRCA2", "c.100delinsACG", window=30)
        self.assertEqual(r["mut_type"], "delins")
        self.assertEqual(len(r["mut_dna_coding"]), 32)

    def test_del3_ins1(self):
        """c.100_102delinsT — delete 3, insert 1 → net -2."""
        r = get_mutation_dna_window("BRCA2", "c.100_102delinsT", window=30)
        self.assertEqual(len(r["mut_dna_coding"]), 28)


@requires_ref_data
class TestGenomeWindowSizes(unittest.TestCase):

    def test_small_window(self):
        r = get_mutation_dna_window("BRAF", "c.1799T>A", window=10)
        self.assertEqual(len(r["ref_dna_coding"]), 10)

    def test_large_window(self):
        r = get_mutation_dna_window("BRAF", "c.1799T>A", window=500)
        self.assertEqual(len(r["ref_dna_coding"]), 500)

    def test_odd_window(self):
        r = get_mutation_dna_window("BRAF", "c.1799T>A", window=33)
        self.assertEqual(len(r["ref_dna_coding"]), 33)

    def test_window_1_returns_single_base(self):
        r = get_mutation_dna_window("BRAF", "c.1799T>A", window=1)
        self.assertEqual(len(r["ref_dna_coding"]), 1)
        self.assertEqual(r["ref_dna_coding"], "T")
        self.assertEqual(r["mut_dna_coding"], "A")

    def test_nested_windows(self):
        """Smaller window should be a substring of the larger window."""
        r_small = get_mutation_dna_window("BRAF", "c.1799T>A", window=20)
        r_large = get_mutation_dna_window("BRAF", "c.1799T>A", window=60)
        self.assertIn(r_small["ref_dna_genomic"], r_large["ref_dna_genomic"])
        self.assertIn(r_small["mut_dna_genomic"], r_large["mut_dna_genomic"])


@requires_ref_data
class TestGenomeErrorHandling(unittest.TestCase):

    def test_unknown_gene(self):
        with self.assertRaises(ValueError):
            get_mutation_dna_window("FAKEGENE123", "c.100A>G")

    def test_position_too_large(self):
        with self.assertRaises(ValueError):
            get_mutation_dna_window("BRAF", "c.999999A>G")

    def test_position_zero(self):
        with self.assertRaises(ValueError):
            get_mutation_dna_window("BRAF", "c.0A>G")

    def test_unparseable_mutation(self):
        with self.assertRaises(ValueError):
            get_mutation_dna_window("BRAF", "not_a_mutation")


@requires_ref_data
class TestGenomeCaseInsensitivity(unittest.TestCase):

    def test_lowercase_gene(self):
        r = get_mutation_dna_window("braf", "c.1799T>A", window=20)
        self.assertEqual(r["gene"], "BRAF")
        self.assertTrue(r["ref_matches_genome"])

    def test_mixed_case_gene(self):
        r = get_mutation_dna_window("Braf", "c.1799T>A", window=20)
        self.assertEqual(r["gene"], "BRAF")

    def test_gene_with_spaces(self):
        r = get_mutation_dna_window("  BRAF  ", "c.1799T>A", window=20)
        self.assertEqual(r["gene"], "BRAF")


@requires_ref_data
class TestGenomeConsistency(unittest.TestCase):
    """Cross-checks between fields that must hold for any result."""

    def test_revcomp_relationship_multiple_genes(self):
        """For all - strand results, coding = revcomp(genomic)."""
        cases = [("BRAF", "c.1799T>A"), ("TP53", "c.818G>A"), ("KRAS", "35G>A")]
        for gene, mut in cases:
            r = get_mutation_dna_window(gene, mut, window=30)
            if r["strand"] == "-":
                self.assertEqual(r["ref_dna_coding"], _revcomp(r["ref_dna_genomic"]),
                                 f"{gene} ref coding/genomic mismatch")
                self.assertEqual(r["mut_dna_coding"], _revcomp(r["mut_dna_genomic"]),
                                 f"{gene} mut coding/genomic mismatch")

    def test_sub_exactly_one_diff(self):
        """For substitutions, ref and mut differ at exactly 1 position."""
        cases = [("BRAF", "c.1799T>A"), ("TP53", "c.818G>A"), ("BRCA2", "c.7397T>G")]
        for gene, mut in cases:
            r = get_mutation_dna_window(gene, mut, window=30)
            diffs = sum(a != b for a, b in zip(r["ref_dna_coding"], r["mut_dna_coding"]))
            self.assertEqual(diffs, 1, f"{gene} {mut}: expected 1 diff, got {diffs}")

    def test_genomic_positions_sorted(self):
        cases = [("BRAF", "c.1799T>A"), ("CFTR", "c.1521_1523delCTT"),
                 ("BRCA1", "c.5266dupC")]
        for gene, mut in cases:
            r = get_mutation_dna_window(gene, mut, window=20)
            self.assertEqual(r["genomic_positions"], sorted(r["genomic_positions"]))


@requires_ref_data
class TestGenomeReturnFields(unittest.TestCase):

    def test_all_keys_present(self):
        r = get_mutation_dna_window("BRAF", "c.1799T>A", window=20)
        expected = {
            "gene", "mutation", "uniprot_id", "chrom", "strand",
            "mut_type", "cds_pos", "cds_end_pos", "intronic_offset",
            "genomic_positions", "ref_matches_genome",
            "window_start", "window_end",
            "ref_dna_coding", "mut_dna_coding",
            "ref_dna_genomic", "mut_dna_genomic",
            "warnings",
        }
        self.assertEqual(set(r.keys()), expected)

    def test_field_types(self):
        r = get_mutation_dna_window("BRAF", "c.1799T>A", window=20)
        self.assertIsInstance(r["gene"], str)
        self.assertIn(r["strand"], ("+", "-"))
        self.assertIsInstance(r["cds_pos"], int)
        self.assertIsInstance(r["genomic_positions"], list)
        self.assertIsInstance(r["ref_matches_genome"], bool)
        self.assertIsInstance(r["warnings"], list)

    def test_dna_is_uppercase_acgt(self):
        r = get_mutation_dna_window("BRAF", "c.1799T>A", window=50)
        for field in ["ref_dna_coding", "mut_dna_coding",
                      "ref_dna_genomic", "mut_dna_genomic"]:
            self.assertTrue(all(c in "ACGT" for c in r[field]),
                            f"{field} contains non-ACGT characters")


if __name__ == "__main__":
    unittest.main()
