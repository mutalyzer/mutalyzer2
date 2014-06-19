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
    Test the mytalyzer.describe module.
    """

    def _single_variant(self, sample, expected):
        """
        General single variant test.
        """
        reference = "ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT"

        result = describe.describe(reference, sample)
        assert result[0].type == expected[0]
        assert result[0].start == expected[1]
        assert result[0].end == expected[2]
        assert result[0].sample_start == expected[3]
        assert result[0].sample_end == expected[4]

    def test1(self):
        """
        Test 1.
        """
        result = describe.describe(
            'ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA',
            'ATGATTTGATCAGATACATGTGATACCGGTAGTTAGGACAA')
        description = describe.allele_description(result)
        assert description == '[5_6insTT;17del;26A>C;35dup]'

    def test2(self):
        """
        Test 2.
        """
        result = describe.describe(
            'TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTACTGTG',
            'TAAGCACCAGGAGTCCATGAAGAAGCTGGATCCTCCCATGGAATCCCCTACTCTACTGTG')
        description = describe.allele_description(result)
        assert description == '[26A>C;30C>A;35G>C]'

    def test3(self):
        """
        Test 3.
        """
        result = describe.describe(
            'TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTA',
            'TAAGCACCAGGAGTCCATGAAGAAGCCATGTCCTGCCATGGAATCCCCTACTCTA')
        description = describe.allele_description(result)
        assert description == '[26_29inv;30C>G]'

    def test4(self):
        """
        Test 4.
        """
        result = describe.describe(
            'TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTA',
            'TAAGCACCAGGAGTCCATGAAGAAGCCATGTCCTGCCATGAATCCCCTACTCTA')
        description = describe.allele_description(result)
        assert description == '[26_29inv;30C>G;41del]'

    def test5(self):
        """
        Test 5.
        """
        self._single_variant("ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("none", 0, 0, 0, 0))

    def test6(self):
        """
        Test 6.
        """
        self._single_variant("ACGTCGGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("subst", 7, 7, 7, 7))

    def test7(self):
        """
        Test 7.
        """
        self._single_variant("ACGTCGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("del", 7, 7, 6, 7))

    def test8(self):
        """
        Test 8.
        """
        self._single_variant("ACGTCGTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("del", 7, 8, 6, 7))

    def test9(self):
        """
        Test 9.
        """
        self._single_variant("ACGTCGCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("ins", 6, 7, 7, 7))

    def test10(self):
        """
        Test 10.
        """
        self._single_variant("ACGTCGCCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("ins", 6, 7, 7, 8))

    def test11(self):
        """
        Test 11.
        """
        self._single_variant("ACGTCGAATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("dup", 7, 7, 8, 8))

    def test12(self):
        """
        Test 12.
        """
        self._single_variant("ACGTCGAGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("dup", 6, 7, 8, 9))

    def test13(self):
        """
        Test 13.
        """
        self._single_variant("ACGTCGACGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("dup", 5, 7, 8, 10))

    def test14(self):
        """
        Test 14.
        """
        self._single_variant("ACGTCGCGAATCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("inv", 7, 11, 7, 11))

    def test15(self):
        """
        Test 15.
        """
        self._single_variant("ACGTCGCCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            ("delins", 7, 7, 7, 8))

    def test16(self):
        """
        Test 16.
        """
        self._single_variant("ACGTCGATTCGCTAGCTTCGTTTTGATAGATAGAGATATAGAGAT",
            ("delins", 21, 23, 21, 24))
