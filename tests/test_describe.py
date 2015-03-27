"""
Tests for the mutalyzer.describe module.
"""


from __future__ import unicode_literals

#import logging; logging.basicConfig()
import os

import mutalyzer
from mutalyzer import describe

from utils import MutalyzerTest


class TestDescribe(MutalyzerTest):
    """
    Test the mutalyzer.describe module.
    """

    def _single_variant(self, sample, expected):
        """
        General single variant test.
        """
        reference = "ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT"

        result = describe.describe_dna(reference, sample)
        assert result[0].type == expected[0]
        assert result[0].start == expected[1]
        assert result[0].end == expected[2]
        assert result[0].sample_start == expected[3]
        assert result[0].sample_end == expected[4]
        assert result[0].deleted[0].sequence == expected[5]
        assert result[0].inserted[0].sequence == expected[6]
        assert unicode(result[0]) == expected[7]

    def test1(self):
        """
        Test 1.
        """
        result = describe.describe_dna(
            'ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA',
            'ATGATTTGATCAGATACATGTGATACCGGTAGTTAGGACAA')
        assert unicode(result) == '[5_6insTT;17del;26A>C;35dup]'

    def test2(self):
        """
        Test 2.
        """
        result = describe.describe_dna(
            'TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTACTGTG',
            'TAAGCACCAGGAGTCCATGAAGAAGCTGGATCCTCCCATGGAATCCCCTACTCTACTGTG')
        assert unicode(result) == '[26A>C;30C>A;35G>C]'

    def test3(self):
        """
        Test 3.
        """
        result = describe.describe_dna(
            'TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTA',
            'TAAGCACCAGGAGTCCATGAAGAAGCCATGTCCTGCCATGGAATCCCCTACTCTA')
        assert unicode(result) == '[26_29inv;30C>G]'

    def test4(self):
        """
        Test 4.
        """
        result = describe.describe_dna(
            'TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTA',
            'TAAGCACCAGGAGTCCATGAAGAAGCCATGTCCTGCCATGAATCCCCTACTCTA')
        assert unicode(result) == '[26_29inv;30C>G;41del]'

    def test5(self):
        """
        Test 5.
        """
        self._single_variant("ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("none", 0, 0, 0, 0, "", "", "="))

    def test6(self):
        """
        Test 6.
        """
        self._single_variant("ACGTCGGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("subst", 7, 7, 7, 7, "A", "G", "7A>G"))

    def test7(self):
        """
        Test 7.
        """
        self._single_variant("ACGTCGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("del", 7, 7, 6, 7, "A", "", "7del"))

    def test8(self):
        """
        Test 8.
        """
        self._single_variant("ACGTCGTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("del", 7, 8, 6, 7, "AT", "", "7_8del"))

    def test9(self):
        """
        Test 9.
        """
        self._single_variant("ACGTCGCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("ins", 6, 7, 7, 7, "", "C", "6_7insC"))

    def test10(self):
        """
        Test 10.
        """
        self._single_variant("ACGTCGCCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("ins", 6, 7, 7, 8, "", "CC", "6_7insCC"))

    def test11(self):
        """
        Test 11.
        """
        self._single_variant("ACGTCGAATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("dup", 7, 7, 8, 8, "", "A", "7dup"))

    def test12(self):
        """
        Test 12.
        """
        self._single_variant("ACGTCGAGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("dup", 6, 7, 8, 9, "", "GA", "6_7dup"))

    def test13(self):
        """
        Test 13.
        """
        self._single_variant("ACGTCGACGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("dup", 5, 7, 8, 10, "", "CGA", "5_7dup"))

    def test14(self):
        """
        Test 14.
        """
        self._single_variant("ACGTCGCGAATCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("inv", 7, 11, 7, 11, "ATTCG", "CGAAT", "7_11inv"))

    def test15(self):
        """
        Test 15.
        """
        self._single_variant("ACGTCGCCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("delins", 7, 7, 7, 8, "A", "CC", "7delinsCC"))

    def test16(self):
        """
        Test 16.
        """
        self._single_variant("ACGTCGATTCGCTAGCTTCGTTTTGATAGATAGAGATATAGAGAT",
            ("delins", 21, 23, 21, 24, "GGG", "TTTT", "21_23delinsTTTT"))

    def test17(self):
        """
        Test 17.
        """
        self._single_variant("ACGTCTCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("inv", 6, 7, 6, 7, "GA", "TC", "6_7inv"))

    def test18(self):
        """
        Test 18.
        """
        self._single_variant("ACGTCGTCTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("delins", 7, 8, 7, 8, "AT", "TC", "7_8delinsTC"))
