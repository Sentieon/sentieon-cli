"""Golden-value tests for ``sentieon_cli._indel2cnv_utils``.

Expected values were captured from ``Levenshtein.ratio`` (v0.27.3) and
``scipy.signal.find_peaks`` (v1.17.1).  Keeping the deps as test-time
imports would defeat the point of this module, so the reference
outputs are frozen as literals here.
"""

import os
import sys

sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

import pytest

from sentieon_cli._indel2cnv_utils import find_peaks, ratio


LONG_A = (
    "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACC"
    "CTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATTTGAAGGAGTCTAGCAGCCGCAGT"
    "AAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATG"
    "ACCCTCTCGTCATAAAACCTTTCTACTATGTGTTCCGCAAGAATCAACAACTACAATGGCGC"
    "GTCGTGAATAACGCGACGGCTGAGACGAACGGCGCGTGAATGAAGCGCTTAAACAGCTCAGG"
    "AGCCAGTCCCCTACGTCGCATATCCTGGCCACTGGAGGTGAAGCGAATGGTATCGATACGTA"
    "GGAGGTGTGCCTTCGTAGGCTGTTTCTCAGGACGCCCAACTATTCTTTCCAATCCTACATCT"
    "GTTTCTTGCGTCGTAGCGGGACCCTCCATTGTTACTTATTAGGTTCTCGTTATGTCTCATAA"
    "TCTC"
)
LONG_B = (
    "AAGCCCAATAAACCACTCTGCCTGGGGGAATAGGGTTATAGGCAACGACATGTGCGGGCTAC"
    "CCTTGCGACAGTGACGCTTTCGCCAGTTGCCTAAACCTATTTAAAGGACTCTAGTAGCCGCA"
    "ATAAGGCCCAATACCTCGTCGGTGTTACCAGACCAAACAAGACGTGCTCTTCAATGTTTAAA"
    "TGACCCTCTCCTCATAAAACCTTTCTACTATGTGTTCCGCAAGAATCAACAACTGACCAGGG"
    "CGCGTCGTGAATAACGCGACGGCTGAGACGAACGGCGCGTGAATGATAGCGCTTAAACAGCT"
    "CAGGAGCCAGTCACCTACGTCGCATATCCTGGCCACTGGAGGTGAAGCGAATGGTATCGATA"
    "CGTAGGAGGTGTGCCGTCGTACGCTGTTGTCTCAGGACGCCCAACTATTCTTTCCAATCCTA"
    "CATCTGTTTCTTGCGTCGTAGCGGGACCCTCCATTGTTACTTATTAGGTTCTCGTTCAGTCT"
    "CATAATCTC"
)
LONG_C = (
    "TCACAGAAGTGAGATTATGTCTCGTTTGGCAGTCTTGATGCTCGGGGGACACTTCTTTAAGC"
    "TCGGTGTGGTGGGCACGACCCTGGACGCGCGACGAAGCTAAGTTTGCAGTAATTAACCGACA"
    "TCTTTGTGAACCGACCCACATTTGACGGTACGCTACCGCAACGGTATGTGTTAATGGAACAG"
    "ACTTGCTTATGTGGACGTTGTATAGGGATATTACGTTACGCGTTAACCGATACATACTGGTT"
    "TCTCTCCAGTGGAGGTCTTGGTTGCCTCTAGTTTCTACGATATACTCATGGTAGTGTAACGC"
    "ATAATCGAAGAGGGTCCTCCCATCTCCTGTGATGCATGGTGTGCTTACTGGGATGAATGCGC"
    "CGCAAGTAGCAGGTCCCGGCGTGGATACCTGATAGATGGTGACTAGCATGTACAAGTAACCT"
    "TGTCTATTGAGCTTCGAGGATGCATACAAGCCCACCCGCAGCCGCAACAGCGACGACTAATT"
    "GATC"
)
RATIO_AB = 0.9572139303482587
RATIO_AC = 0.636

TANDEM_UNIT = "ACGTACGAT"
TANDEM = TANDEM_UNIT * 6
TANDEM_VAR = (
    TANDEM[: len(TANDEM) // 2] + "T" + TANDEM[len(TANDEM) // 2 + 1:]
)
TANDEM_RATIO = 0.9814814814814815


class TestRatio:
    """``ratio`` must match ``Levenshtein.ratio`` on the captured cases."""

    @pytest.mark.parametrize(
        "a,b,expected",
        [
            ("ACGTACGT", "ACGTACGT", 1.0),
            ("", "", 1.0),
            ("ACGT", "", 0.0),
            ("", "ACGT", 0.0),
            ("AAAA", "CCCC", 0.0),
            ("ACGT", "ACGGT", 0.8888888888888888),
            ("ACGGT", "ACGT", 0.8888888888888888),
            ("ACGT", "AGGT", 0.75),
            ("ACGTACGT" * 4, "ACGTACGT" * 4, 1.0),
            ("ACGTACGT" * 4, "CGTACGTA" * 4, 0.96875),
            ("AACCGGTT", "AACCGGTG", 0.875),
            (TANDEM, TANDEM_VAR, TANDEM_RATIO),
        ],
    )
    def test_golden_values(self, a, b, expected):
        assert ratio(a, b) == pytest.approx(expected, abs=1e-12)

    def test_long_similar_dna(self):
        assert ratio(LONG_A, LONG_B) == pytest.approx(RATIO_AB, abs=1e-12)

    def test_long_unrelated_dna(self):
        assert ratio(LONG_A, LONG_C) == pytest.approx(RATIO_AC, abs=1e-12)

    def test_symmetry(self):
        assert ratio(LONG_A, LONG_B) == ratio(LONG_B, LONG_A)

    def test_threshold_pass_returns_exact_value(self):
        # Real ratio is above threshold — must return the exact value.
        assert ratio(LONG_A, LONG_B, threshold=0.9) == pytest.approx(
            RATIO_AB, abs=1e-12
        )

    def test_threshold_fail_returns_below_threshold(self):
        # Real ratio is below threshold — returned value must be < threshold.
        result = ratio("AAAA", "CCCC", threshold=0.5)
        assert result < 0.5

    def test_threshold_zero_disables_abort(self):
        # Whatever value threshold=0 returns must equal the no-threshold call.
        assert ratio(LONG_A, LONG_C, threshold=0.0) == ratio(LONG_A, LONG_C)

    def test_threshold_very_asymmetric_lengths(self):
        # Length-based upper bound rejects cheaply.
        result = ratio("A" * 10, "A" * 10000, threshold=0.9)
        assert result < 0.9


class TestFindPeaks:
    """``find_peaks`` must match ``scipy.signal.find_peaks(x, distance=d,
    height=h)`` for the call pattern used by indel2cnv.py."""

    @pytest.mark.parametrize(
        "values,distance,height,expected",
        [
            ([0.0, 0.5, 1.0, 0.5, 0.0], 1, 0.9, [2]),
            (
                [0.0, 0.95, 0.0, 0.0, 0.0, 0.97, 0.0],
                2, 0.9, [1, 5],
            ),
            ([0.0, 0.95, 0.0, 0.97, 0.0], 3, 0.9, [3]),
            ([0.0, 0.5, 0.0, 0.6, 0.0], 1, 0.9, []),
            ([0.0, 0.95, 0.95, 0.0], 1, 0.9, [1]),
            ([1.0, 0.9, 0.8], 1, 0.9, []),
            ([0.8, 0.9, 1.0], 1, 0.9, []),
            (
                [0.3, 0.5, 0.6, 0.95, 0.7, 0.6, 0.5, 0.96, 0.4, 0.3],
                2, 0.9, [3, 7],
            ),
        ],
    )
    def test_golden_values(self, values, distance, height, expected):
        assert find_peaks(
            values, distance=distance, height=height
        ) == expected

    def test_empty_input(self):
        assert find_peaks([], distance=1, height=0.9) == []

    def test_two_element_input(self):
        # Neither element has both neighbors; never a peak.
        assert find_peaks([0.0, 1.0], distance=1, height=0.5) == []

    def test_left_edge_plateau_not_a_peak(self):
        # scipy ignores plateaus without a strictly lower left neighbor.
        assert find_peaks([0.95, 0.95, 0.0], distance=1, height=0.9) == []

    def test_right_edge_plateau_not_a_peak(self):
        # scipy ignores plateaus without a strictly lower right neighbor.
        assert find_peaks([0.0, 0.95, 0.95], distance=1, height=0.9) == []

    def test_wider_plateau_midpoint(self):
        # Plateau of width 3 at indices 2..4; midpoint = 3.
        assert find_peaks(
            [0.0, 0.3, 0.95, 0.95, 0.95, 0.3, 0.0],
            distance=1, height=0.9,
        ) == [3]
