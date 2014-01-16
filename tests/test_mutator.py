"""
Tests for the mutalyzer.mutator module.
"""


#import logging; logging.basicConfig()
import re
import os
import random
from nose.tools import *
from Bio.Seq import Seq

import mutalyzer
from mutalyzer.output import Output
from mutalyzer import mutator

from utils import MutalyzerTest


def _seq(length):
    """
    Return random DNA sequence of given length.
    """
    sequence = ''
    for i in range(length):
        sequence += random.choice('ACGT')
    return Seq(sequence)


class TestMutator(MutalyzerTest):
    """
    Test the mutator module.
    """
    def setup(self):
        super(TestMutator, self).setup()
        self.output = Output(__file__)

    def _mutator(self, sequence):
        """
        Create a Mutator instance for a given sequence.
        """
        return mutator.Mutator(sequence, self.output)

    def test_shift_no_change(self):
        """
        No change, no shifts.
        """
        l = 10
        m = self._mutator(_seq(l))
        # Numbering is 1-based
        for i in range(1, l + 1):
            assert_equal(m.shift(i), i)

    def test_shift_del_example(self):
        """
        Example of g.2del.
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(2, 2)
        assert_equal(m.shift(1), 1)
        assert_equal(m.shift(2), 2)
        assert_equal(m.shift(3), 2)

    def test_shift_del(self):
        """
        Starting from the deleted position (not included), shift -1.
        """
        l = 10
        for d in range(1, l + 1):
            m = self._mutator(_seq(l))
            m.deletion(d, d)
            for p in range(1, d + 1):
                assert_equal(m.shift(p), p)
            for p in range(d + 1, l + 1):
                assert_equal(m.shift(p), p - 1)

    def test_shift_del2(self):
        """
        Starting from the deleted positions (not included), shift -2.
        """
        l = 10
        for d in range(1, l):
            m = self._mutator(_seq(l))
            m.deletion(d, d + 1)
            for p in range(1, d + 2):
                assert_equal(m.shift(p), p)
            for p in range(d + 2, l + 1):
                assert_equal(m.shift(p), p - 2)

    def test_shift_ins_example(self):
        """
        Example of g.2_3insA.
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'A')
        assert_equal(m.shift(1), 1)
        assert_equal(m.shift(2), 2)
        assert_equal(m.shift(3), 4)

    def test_shift_ins(self):
        """
        Starting from the interbase insertion position, shift +1.
        """
        l = 10
        for i in range(0, l + 1):
            m = self._mutator(_seq(l))
            m.insertion(i, 'T')
            for p in range(1, i + 1):
                assert_equal(m.shift(p), p)
            for p in range(i + 1, l + 1):
                assert_equal(m.shift(p), p + 1)

    def test_shift_ins2(self):
        """
        Starting from the interbase insertion position, shift +2.
        """
        l = 10
        for i in range(0, l + 1):
            m = self._mutator(_seq(l))
            m.insertion(i, 'TT')
            for p in range(1, i + 1):
                assert_equal(m.shift(p), p)
            for p in range(i + 1, l + 1):
                assert_equal(m.shift(p), p + 2)

    def test_shift_sites_no_change(self):
        """
        No change, no shifts.

        @note: Splice sites come in pairs (acceptor and donor site) and the
        numbers are the first, respectively last, position in the exon.

        So in this example we have:     ---======----======-----===---
                                           |    |    |    |     | |
                                           4    9   14    19   25 27
        """
        l = 30
        sites = [4, 9, 14, 19, 25, 27]
        m = self._mutator(_seq(l))
        assert_equal(m.shift_sites(sites), sites)

    def test_shift_sites_acc_del_before(self):
        """
        Deletion in intron directly before exon.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.deletion(13, 13)   # g.13del
        assert_equal(m.shift_sites(sites), [4, 9, 13, 16, 24, 26])

    def test_shift_sites_acc_del_after(self):
        """
        Deletion at first exon position.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.deletion(14, 14)   # g.14del
        assert_equal(m.shift_sites(sites), [4, 9, 14, 16, 24, 26])

    def test_shift_sites_don_del_before(self):
        """
        Deletion at last exon position.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.deletion(17, 17)   # g.17del
        assert_equal(m.shift_sites(sites), [4, 9, 14, 16, 24, 26])

    def test_shift_sites_don_del_after(self):
        """
        Deletion in intron directly after exon.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.deletion(18, 18)   # g.18del
        assert_equal(m.shift_sites(sites), [4, 9, 14, 17, 24, 26])

    def test_shift_sites_acc_del2_before(self):
        """
        Deletion of 2 in intron directly before exon.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.deletion(12, 13)   # g.12_13del
        assert_equal(m.shift_sites(sites), [4, 9, 12, 15, 23, 25])

    def test_shift_sites_acc_del2_on(self):
        """
        Deletion of 2 in intron/exon.

        @note: This hits a splice site, so we don't really support it.
        """
        return

        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.deletion(13, 14)   # g.13_14del
        assert_equal(m.shift_sites(sites), [4, 9, 13, 15, 23, 25])

    def test_shift_sites_acc_del2_after(self):
        """
        Deletion of 2 at first exon position.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.deletion(14, 15)   # g.14_15del
        assert_equal(m.shift_sites(sites), [4, 9, 14, 15, 23, 25])

    def test_shift_sites_don_del2_before(self):
        """
        Deletion of 2 at last exon positions.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.deletion(16, 17)   # g.16_17del
        assert_equal(m.shift_sites(sites), [4, 9, 14, 15, 23, 25])

    def test_shift_sites_don_del2_on(self):
        """
        Deletion of 2 in exon/intron.

        @note: This hits a splice site, so we don't really support it.
        """
        return

        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.deletion(17, 18)   # g.17_18del
        assert_equal(m.shift_sites(sites), [4, 9, 14, 16, 23, 25])

    def test_shift_sites_don_del2_after(self):
        """
        Deletion of 2 in intron directly after exon.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.deletion(18, 19)   # g.18_19del
        assert_equal(m.shift_sites(sites), [4, 9, 14, 17, 23, 25])

    def test_shift_sites_acc_ins_before(self):
        """
        Insertion 1 position before intron/exon boundary.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(12, 'A')   # g.12_13insA
        assert_equal(m.shift_sites(sites), [4, 9, 15, 18, 26, 28])

    def test_shift_sites_acc_ins_on(self):
        """
        Insertion in intron/exon boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(13, 'A')   # g.13_14insA
        assert_equal(m.shift_sites(sites), [4, 9, 14, 18, 26, 28])

    def test_shift_sites_first_acc_ins_on(self):
        """
        Insertion in first intron/exon boundary not be included.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(3, 'A')   # g.3_4insA
        assert_equal(m.shift_sites(sites), [5, 10, 15, 18, 26, 28])

    def test_shift_sites_acc_ins_after(self):
        """
        Insertion 1 position after intron/exon boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(14, 'A')   # g.14_15insA
        assert_equal(m.shift_sites(sites), [4, 9, 14, 18, 26, 28])

    def test_shift_sites_don_ins_before(self):
        """
        Insertion 1 position before exon/intron boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(16, 'A')   # g.16_17insA
        assert_equal(m.shift_sites(sites), [4, 9, 14, 18, 26, 28])

    def test_shift_sites_don_ins_on(self):
        """
        Insertion in exon/intron boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(17, 'A')   # g.17_18insA
        assert_equal(m.shift_sites(sites), [4, 9, 14, 18, 26, 28])

    def test_shift_sites_last_don_ins_on(self):
        """
        Insertion in last exon/intron boundary should not be included.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(27, 'A')   # g.27_28insA
        assert_equal(m.shift_sites(sites), [4, 9, 14, 17, 25, 27])

    def test_shift_sites_don_ins_after(self):
        """
        Insertion 1 position after exon/intron boundary.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(18, 'A')   # g.18_19insA
        assert_equal(m.shift_sites(sites), [4, 9, 14, 17, 26, 28])

    def test_shift_sites_acc_ins2_before(self):
        """
        Insertion of 2 1 position before intron/exon boundary.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(12, 'AT')   # g.12_13insAT
        assert_equal(m.shift_sites(sites), [4, 9, 16, 19, 27, 29])

    def test_shift_sites_first_acc_ins2_on(self):
        """
        Insertion of 2 in last exon/intron boundary should not be included.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(3, 'AT')   # g.3_4insAT
        assert_equal(m.shift_sites(sites), [6, 11, 16, 19, 27, 29])

    def test_shift_sites_acc_ins2_after(self):
        """
        Insertion of 2 1 position after intron/exon boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(14, 'AT')   # g.14_15insAT
        assert_equal(m.shift_sites(sites), [4, 9, 14, 19, 27, 29])

    def test_shift_sites_don_ins2_before(self):
        """
        Insertion of 2 1 position before exon/intron boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(16, 'AT')   # g.16_17insAT
        assert_equal(m.shift_sites(sites), [4, 9, 14, 19, 27, 29])

    def test_shift_sites_last_don_ins2_on(self):
        """
        Insertion of 2 in last exon/intron boundary should not be included.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(27, 'AT')   # g.27_28insAT
        assert_equal(m.shift_sites(sites), [4, 9, 14, 17, 25, 27])

    def test_shift_sites_don_ins2_after(self):
        """
        Insertion of 2 1 position after exon/intron boundary.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(18, 'AT')   # g.18_19insAT
        assert_equal(m.shift_sites(sites), [4, 9, 14, 17, 27, 29])

    def test_shift_sites_acc_ins3_before(self):
        """
        Insertion of 3 1 position before intron/exon boundary.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(12, 'ATT')   # g.12_13insATT
        assert_equal(m.shift_sites(sites), [4, 9, 17, 20, 28, 30])

    def test_shift_sites_acc_ins3_on(self):
        """
        Insertion of 3 in intron/exon boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(13, 'ATT')   # g.13_14insATT
        assert_equal(m.shift_sites(sites), [4, 9, 14, 20, 28, 30])

    def test_shift_sites_first_acc_ins3_on(self):
        """
        Insertion of 3 in first intron/exon boundary should not be included.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(3, 'ATT')   # g.3_4insATT
        assert_equal(m.shift_sites(sites), [7, 12, 17, 20, 28, 30])

    def test_shift_sites_acc_ins3_after(self):
        """
        Insertion of 3 1 position after intron/exon boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(14, 'ATT')   # g.14_15insATT
        assert_equal(m.shift_sites(sites), [4, 9, 14, 20, 28, 30])

    def test_shift_sites_don_ins3_before(self):
        """
        Insertion of 3 1 position before exon/intron boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(16, 'ATT')   # g.16_17insATT
        assert_equal(m.shift_sites(sites), [4, 9, 14, 20, 28, 30])

    def test_shift_sites_don_ins3_on(self):
        """
        Insertion of 3 in exon/intron boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(17, 'ATT')   # g.17_18insATT
        assert_equal(m.shift_sites(sites), [4, 9, 14, 20, 28, 30])

    def test_shift_sites_last_don_ins3_on(self):
        """
        Insertion of 3 in last exon/intron boundary should not be included.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(27, 'ATT')   # g.27_28insATT
        assert_equal(m.shift_sites(sites), [4, 9, 14, 17, 25, 27])

    def test_shift_sites_don_ins3_after(self):
        """
        Insertion of 3 1 position after exon/intron boundary.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insertion(18, 'ATT')   # g.18_19insATT
        assert_equal(m.shift_sites(sites), [4, 9, 14, 17, 28, 30])

    def test_shift_sites_adj_del_before1(self):
        """
        Adjacent exons: deletion at second-last position of first exon.

        @note: In this example we have adjacent exons (like e.g. in RNA),
        which looks like this (the square brackets [ and ] are part of the
        exons):
                    ---[====][======][========]---
                       |   /  \    /  \       |
                       4  9   10  17  18      27
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.deletion(16, 16)   # g.16del
        assert_equal(m.shift_sites(sites), [4, 9, 10, 16, 17, 26])

    def test_shift_sites_adj_del_before(self):
        """
        Adjacent exons: deletion at last position of first exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.deletion(17, 17)   # g.17del
        assert_equal(m.shift_sites(sites), [4, 9, 10, 16, 17, 26])

    def test_shift_sites_adj_del_after(self):
        """
        Adjacent exons: deletion at first position of second exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.deletion(18, 18)   # g.18del
        assert_equal(m.shift_sites(sites), [4, 9, 10, 17, 18, 26])

    def test_shift_sites_adj_del_after1(self):
        """
        Adjacent exons: deletion at second position of second exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.deletion(19, 19)   # g.19del
        assert_equal(m.shift_sites(sites), [4, 9, 10, 17, 18, 26])

    def test_shift_sites_adj_ins_before(self):
        """
        Adjacent exons: insertion 1 position before exon/exon boundary.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.insertion(16, 'A')   # g.16_17insA
        assert_equal(m.shift_sites(sites), [4, 9, 10, 18, 19, 28])

    def test_shift_sites_adj_ins_on(self):
        """
        Adjacent exons: insertion at exon/exon boundary.

        @note: This insertion could be seen as being
               1) at the end of the first exon, or
               2) at the start of the second exon.
               Both would probably be 'correct', but we would like consistent
               results. Therefore, we stick to the first option.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.insertion(17, 'A')   # g.17_18insA
        assert_equal(m.shift_sites(sites), [4, 9, 10, 18, 19, 28])

    def test_shift_sites_adj_ins_after(self):
        """
        Adjacent exons: insertion 1 position after exon/exon boundary.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.insertion(18, 'A')   # g.18_19insA
        assert_equal(m.shift_sites(sites), [4, 9, 10, 17, 18, 28])

    def test_shift_sites_adj_del2_before1(self):
        """
        Adjacent exons: deletion of 2 at second-last position of first exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.deletion(15, 16)   # g.15_16del
        assert_equal(m.shift_sites(sites), [4, 9, 10, 15, 16, 25])

    def test_shift_sites_adj_del2_before(self):
        """
        Adjacent exons: deletion of 2 at last position of first exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.deletion(16, 17)   # g.16_17del
        assert_equal(m.shift_sites(sites), [4, 9, 10, 15, 16, 25])

    def test_shift_sites_adj_del2_on(self):
        """
        Adjacent exons: deletion of 2 at exon/exon boundary.

        @todo: This is a special case of bug #????. Once fixed, the two
               exons will be joined to one new exon.
        """
        return

        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.deletion(17, 18)   # g.17_18del
        assert_equal(m.shift_sites(sites), [4, 9, 10, 16, 17, 25])

    def test_shift_sites_adj_del2_after(self):
        """
        Adjacent exons: deletion of 2 at first position of second exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.deletion(18, 19)   # g.18_19del
        assert_equal(m.shift_sites(sites), [4, 9, 10, 17, 18, 25])

    def test_shift_sites_adj_del2_after1(self):
        """
        Adjacent exons: deletion of 2 at second position of second exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.deletion(19, 20)   # g.19_20del
        assert_equal(m.shift_sites(sites), [4, 9, 10, 17, 18, 25])

    def test_shift_sites_adj_ins2_before(self):
        """
        Adjacent exons: insertion of 2 1 position before exon/exon boundary.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.insertion(16, 'AT')   # g.16_17insAT
        assert_equal(m.shift_sites(sites), [4, 9, 10, 19, 20, 29])

    def test_shift_sites_adj_ins2_on(self):
        """
        Adjacent exons: insertion of 2 at exon/exon boundary.

        @note: This insertion could be seen as being
               1) at the end of the first exon, or
               2) at the start of the second exon.
               Both would probably be 'correct', but we would like consistent
               results. Therefore, we stick to the first option.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.insertion(17, 'AT')   # g.17_18insAT
        assert_equal(m.shift_sites(sites), [4, 9, 10, 19, 20, 29])

    def test_shift_sites_adj_ins2_after(self):
        """
        Adjacent exons: insertion of 2 1 position after exon/exon boundary.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.insertion(18, 'AT')   # g.18_19insAT
        assert_equal(m.shift_sites(sites), [4, 9, 10, 17, 18, 29])

    def test_del(self):
        """
        Simple deletion 2del.
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(2, 2)
        assert_equal(str(m.mutated), str(Seq('ACGATCG')))

    def test_largedel(self):
        """
        Simple large deletion 2_7del.
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(2, 7)
        assert_equal(str(m.mutated), str(Seq('AG')))

    def test_ins(self):
        """
        Simple insertion 2_3insA.
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'A')
        assert_equal(str(m.mutated), str(Seq('ATACGATCG')))

    def test_largeins(self):
        """
        Simple large insertion 2_3insATCG.
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'ATCG')
        assert_equal(str(m.mutated), str(Seq('ATATCGCGATCG')))

    def test_sub(self):
        """
        Simple substitution 3C>G.
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.substitution(3, 'G')
        assert_equal(str(m.mutated), str(Seq('ATGGATCG')))

    def test_adjecent_del_sub_1(self):
        """
        Deletion and substitution directly adjecent to each other [2del;3C>G].

        See Trac #83.
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(2, 2)
        m.substitution(3, 'G')
        assert_equal(str(m.mutated), str(Seq('AGGATCG')))

    def test_adjecent_del_sub_2(self):
        """
        Deletion and substitution directly adjecent to each other [3del;2T>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(3, 3)
        m.substitution(2, 'G')
        assert_equal(str(m.mutated), str(Seq('AGGATCG')))

    def test_near_adjecent_del_sub_1(self):
        """
        Deletion and substitution almost adjecent to each other [2del;4G>T].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(2, 2)
        m.substitution(4, 'T')
        assert_equal(str(m.mutated), str(Seq('ACTATCG')))

    def test_near_adjecent_del_sub_2(self):
        """
        Deletion and substitution almost adjecent to each other [4del;2T>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(4, 4)
        m.substitution(2, 'G')
        assert_equal(str(m.mutated), str(Seq('AGCATCG')))

    def test_adjecent_largedel_sub_1(self):
        """
        Large deletion and substitution directly adjecent to each other
        [2_6del;7C>T].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(2, 6)
        m.substitution(7, 'T')
        assert_equal(str(m.mutated), str(Seq('ATG')))

    def test_adjecent_largedel_sub_2(self):
        """
        Large deletion and substitution directly adjecent to each other
        [3_7del;2T>C].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(3, 7)
        m.substitution(2, 'C')
        assert_equal(str(m.mutated), str(Seq('ACG')))

    def test_near_adjecent_largedel_sub_1(self):
        """
        Large deletion and substitution almost adjecent to each other [2_5del;7C>T].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(2, 5)
        m.substitution(7, 'T')
        assert_equal(str(m.mutated), str(Seq('ATTG')))

    def test_near_adjecent_largedel_sub_2(self):
        """
        Large deletion and substitution almost adjecent to each other [4_7del;2T>C].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(4, 7)
        m.substitution(2, 'C')
        assert_equal(str(m.mutated), str(Seq('ACCG')))

    def test_adjectent_del_ins_1(self):
        """
        Deletion and insertion adjecent to each other [2del;2_3insG].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(2, 2)
        m.insertion(2, 'G')
        assert_equal(str(m.mutated), str(Seq('AGCGATCG')))

    def test_adjectent_del_ins_2(self):
        """
        Deletion and insertion adjecent to each other [3del;2_3insA].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(3, 3)
        m.insertion(2, 'A')
        assert_equal(str(m.mutated), str(Seq('ATAGATCG')))

    def test_near_adjectent_del_ins(self):
        """
        Deletion and insertion almost adjecent to each other [2del;3_4insG].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(2, 2)
        m.insertion(3, 'T')
        assert_equal(str(m.mutated), str(Seq('ACTGATCG')))

    def test_adjecent_ins_sub_1(self):
        """
        Insertion and substitution directly adjecent to each other
        [2_3insA;3C>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'A')
        m.substitution(3, 'G')
        assert_equal(str(m.mutated), str(Seq('ATAGGATCG')))

    def test_adjecent_ins_sub_2(self):
        """
        Insertion and substitution directly adjecent to each other
        [2_3insA;2T>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'A')
        m.substitution(2, 'G')
        assert_equal(str(m.mutated), str(Seq('AGACGATCG')))

    def test_near_adjecent_ins_sub(self):
        """
        Insertion and substitution almost adjecent to each other
        [2_3insA;4C>T].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'A')
        m.substitution(4, 'T')
        assert_equal(str(m.mutated), str(Seq('ATACTATCG')))

    def test_adjecent_largeins_sub_1(self):
        """
        Large insertion and substitution directly adjecent to each other
        [2_3insATCG;3C>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'ATCG')
        m.substitution(3, 'G')
        assert_equal(str(m.mutated), str(Seq('ATATCGGGATCG')))

    def test_adjecent_largeins_sub_2(self):
        """
        Large insertion and substitution directly adjecent to each other
        [2_3insATCG;2T>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'ATCG')
        m.substitution(2, 'G')
        assert_equal(str(m.mutated), str(Seq('AGATCGCGATCG')))

    def test_near_adjecent_largeins_sub(self):
        """
        Large insertion and substitution almost adjecent to each other
        [2_3insATCG;4C>T].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'ATCG')
        m.substitution(4, 'T')
        assert_equal(str(m.mutated), str(Seq('ATATCGCTATCG')))

    def test_adjecent_del_del_1(self):
        """
        Deletion and deletion directly adjecent to each other [2del;3del].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(2, 2)
        m.deletion(3, 3)
        assert_equal(str(m.mutated), str(Seq('AGATCG')))

    def test_adjecent_del_del_2(self):
        """
        Deletion and deletion directly adjecent to each other [3del;2del].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(3, 3)
        m.deletion(2, 2)
        assert_equal(str(m.mutated), str(Seq('AGATCG')))

    def test_adjecent_delins_snp_1(self):
        """
        Delins and deletion directly adjecent to each other [2delinsA;3C>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 2, 'A')
        m.substitution(3, 'G')
        assert_equal(str(m.mutated), str(Seq('AAGGATCG')))

    def test_adjecent_delins_snp_2(self):
        """
        Delins and deletion directly adjecent to each other [3delinsA;2T>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 3, 'A')
        m.substitution(2, 'G')
        assert_equal(str(m.mutated), str(Seq('AGAGATCG')))

    def test_adjecent_largedelins_eq_snp_1(self):
        """
        Large delins and deletion directly adjecent to each other
        [2_6delinsAAAAA;7C>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 6, 'AAAAA')
        m.substitution(7, 'G')
        assert_equal(str(m.mutated), str(Seq('AAAAAAGG')))

    def test_adjecent_largedelins_min_snp_1(self):
        """
        Large delins (min) and deletion directly adjecent to each other
        [2_6delinsAAA;7C>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 6, 'AAA')
        m.substitution(7, 'G')
        assert_equal(str(m.mutated), str(Seq('AAAAGG')))

    def test_adjecent_largedelins_plus_snp_1(self):
        """
        Large delins (plus) and deletion directly adjecent to each other
        [2_6delinsAAAAAAA;7C>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 6, 'AAAAAAA')
        m.substitution(7, 'G')
        assert_equal(str(m.mutated), str(Seq('AAAAAAAAGG')))

    def test_adjecent_largedelins_eq_snp_2(self):
        """
        Large delins and deletion directly adjecent to each other
        [3_7delinsAAAAA;2T>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 7, 'AAAAA')
        m.substitution(2, 'G')
        assert_equal(str(m.mutated), str(Seq('AGAAAAAG')))

    def test_adjecent_largedelins_min_snp_2(self):
        """
        Large delins (min) and deletion directly adjecent to each other
        [3_7delinsAAA;2T>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 7, 'AAA')
        m.substitution(2, 'G')
        assert_equal(str(m.mutated), str(Seq('AGAAAG')))

    def test_adjecent_largedelins_plus_snp_2(self):
        """
        Large delins (plus) and deletion directly adjecent to each other
        [3_7delinsAAAAAAA;2T>G].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 7, 'AAAAAAA')
        m.substitution(2, 'G')
        assert_equal(str(m.mutated), str(Seq('AGAAAAAAAG')))

    def test_adjecent_delins_del_1(self):
        """
        Delins and deletion directly adjecent to each other [2delinsA;3del].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 2, 'A')
        m.deletion(3, 3)
        assert_equal(str(m.mutated), str(Seq('AAGATCG')))

    def test_adjecent_delins_del_2(self):
        """
        Delins and deletion directly adjecent to each other [3delinsA;2del].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 3, 'A')
        m.deletion(2, 2)
        assert_equal(str(m.mutated), str(Seq('AAGATCG')))

    def test_adjecent_largedelins_eq_del_1(self):
        """
        Large delins and deletion directly adjecent to each other
        [2_6delinsAAAAA;7del].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 6, 'AAAAA')
        m.deletion(7, 7)
        assert_equal(str(m.mutated), str(Seq('AAAAAAG')))

    def test_adjecent_largedelins_min_del_1(self):
        """
        Large delins (min) and deletion directly adjecent to each other
        [2_6delinsAAA;7del].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 6, 'AAA')
        m.deletion(7, 7)
        assert_equal(str(m.mutated), str(Seq('AAAAG')))

    def test_adjecent_largedelins_plus_del_1(self):
        """
        Large delins (plus) and deletion directly adjecent to each other
        [2_6delinsAAAAAAA;7del].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 6, 'AAAAAAA')
        m.deletion(7, 7)
        assert_equal(str(m.mutated), str(Seq('AAAAAAAAG')))

    def test_adjecent_largedelins_eq_del_2(self):
        """
        Large delins and deletion directly adjecent to each other
        [3_7delinsAAAAA;2del].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 7, 'AAAAA')
        m.deletion(2, 2)
        assert_equal(str(m.mutated), str(Seq('AAAAAAG')))

    def test_adjecent_largedelins_min_del_2(self):
        """
        Large delins (min) and deletion directly adjecent to each other
        [3_7delinsAAA;2del].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 7, 'AAA')
        m.deletion(2, 2)
        assert_equal(str(m.mutated), str(Seq('AAAAG')))

    def test_adjecent_largedelins_plus_del_2(self):
        """
        Large delins (plus) and deletion directly adjecent to each other
        [3_7delinsAAAAAAA;2del].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 7, 'AAAAAAA')
        m.deletion(2, 2)
        assert_equal(str(m.mutated), str(Seq('AAAAAAAAG')))

    def test_adjectent_delins_ins_1(self):
        """
        Delins and insertion adjecent to each other [2delinsA;2_3insG].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 2, 'A')
        m.insertion(2, 'G')
        assert_equal(str(m.mutated), str(Seq('AAGCGATCG')))

    def test_adjectent_delins_ins_2(self):
        """
        Delins and insertion adjecent to each other [3delinsA;2_3insG].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 3, 'A')
        m.insertion(2, 'G')
        assert_equal(str(m.mutated), str(Seq('ATGAGATCG')))

    def test_adjectent_largedelins_eq_ins_1(self):
        """
        Large delins and insertion adjecent to each other [2_6delinsAAAAA;6_7insG].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 6, 'AAAAA')
        m.insertion(6, 'G')
        assert_equal(str(m.mutated), str(Seq('AAAAAAGCG')))

    def test_adjectent_largedelins_min_ins_1(self):
        """
        Large delins (min) and insertion adjecent to each other [2_6delinsAAA;6_7insG].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 6, 'AAA')
        m.insertion(6, 'G')
        assert_equal(str(m.mutated), str(Seq('AAAAGCG')))

    def test_adjectent_largedelins_plus_ins_1(self):
        """
        Large delins (plus) and insertion adjecent to each other [2_6delinsAAAAAAA;6_7insG].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 6, 'AAAAAAA')
        m.insertion(6, 'G')
        assert_equal(str(m.mutated), str(Seq('AAAAAAAAGCG')))

    def test_adjectent_largedelins_eq_ins_2(self):
        """
        Large delins and insertion adjecent to each other [3_7delinsAAAAA;2_3insG].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 7, 'AAAAA')
        m.insertion(2, 'G')
        assert_equal(str(m.mutated), str(Seq('ATGAAAAAG')))

    def test_adjectent_largedelins_min_ins_2(self):
        """
        Large delins (min) and insertion adjecent to each other [3_7delinsAAA;2_3insG].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 7, 'AAA')
        m.insertion(2, 'G')
        assert_equal(str(m.mutated), str(Seq('ATGAAAG')))

    def test_adjectent_largedelins_plus_ins_2(self):
        """
        Large delins (plus) and insertion adjecent to each other [3_7delinsAAAAAAA;2_3insG].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 7, 'AAAAAAA')
        m.insertion(2, 'G')
        assert_equal(str(m.mutated), str(Seq('ATGAAAAAAAG')))

    def test_adjectent_delins_del_delins(self):
        """
        Delins (deletion) and delins (SNP) adjecent to each other [2_3delinsA;4delinsT].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 3, 'A')
        m.delins(4, 4, 'T')
        assert_equal(str(m.mutated), str(Seq('AATATCG')))

    def test_adjectent_largedelins_plus_delins_1(self):
        """
        Large delins (plus) and delins adjecent to each other [2_6delinsAAAAAAA;7delinsT].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 6, 'AAAAAAA')
        m.delins(7, 7, 'T')
        assert_equal(str(m.mutated), str(Seq('AAAAAAAATG')))

    def test_adjectent_largedelins_plus_delins_2(self):
        """
        Large delins (plus) and delins adjecent to each other [3_7delinsAAAAAAA;2delinsC].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 7, 'AAAAAAA')
        m.delins(2, 2, 'C')
        assert_equal(str(m.mutated), str(Seq('ACAAAAAAAG')))

    def test_adjectent_largedelins_min_delins_1(self):
        """
        Large delins (min) and delins adjecent to each other [2_6delinsAAA;7delinsT].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(2, 6, 'AAA')
        m.delins(7, 7, 'T')
        assert_equal(str(m.mutated), str(Seq('AAAATG')))

    def test_adjectent_largedelins_min_delins_2(self):
        """
        Large delins (min) and delins adjecent to each other [3_7delinsAAA;2delinsC].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delins(3, 7, 'AAA')
        m.delins(2, 2, 'C')
        assert_equal(str(m.mutated), str(Seq('ACAAAG')))

    def test_adjectent_del_dup_1(self):
        """
        Deletion and duplication adjecent to each other [2del;3dup].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(2, 2)
        m.duplication(3, 3)
        assert_equal(str(m.mutated), str(Seq('ACCGATCG')))

    def test_adjectent_del_dup_2(self):
        """
        Deletion and duplication adjecent to each other [3del;2dup].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(3, 3)
        m.duplication(2, 2)
        assert_equal(str(m.mutated), str(Seq('ATTGATCG')))

    def test_adjectent_ins_dup_1(self):
        """
        Insertion and duplication adjecent to each other [2_3insG;3dup].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'G')
        m.duplication(3, 3)
        assert_equal(str(m.mutated), str(Seq('ATGCCGATCG')))

    def test_adjectent_ins_dup_2(self):
        """
        Insertion and duplication adjecent to each other [2_3insG;2dup].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'G')
        m.duplication(2, 2)
        assert_equal(str(m.mutated), str(Seq('ATTGCGATCG')))

    def test_adjectent_ins_ins_1(self):
        """
        Insertion and insertion adjecent to each other [2_3insG;3_4insA].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'G')
        m.insertion(3, 'A')
        assert_equal(str(m.mutated), str(Seq('ATGCAGATCG')))

    def test_adjectent_ins_ins_2(self):
        """
        Insertion and insertion adjecent to each other [3_4insA;2_3insG].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(3, 'A')
        m.insertion(2, 'G')
        assert_equal(str(m.mutated), str(Seq('ATGCAGATCG')))

    def test_ins_ins(self):
        """
        Insertion and insertion at same position [2_3insG;2_3insA].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'G')
        m.insertion(2, 'A')
        assert str(m.mutated) in (str(Seq('ATGACGATCG')), str(Seq('ATAGCGATCG')))

    def test_adjecent_inv_inv_1(self):
        """
        Inversion and inversion directly adjecent to each other [2inv;3inv].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.inversion(2, 2)
        m.inversion(3, 3)
        assert_equal(str(m.mutated), str(Seq('AAGGATCG')))

    def test_adjecent_inv_inv_2(self):
        """
        Inversion and inversion directly adjecent to each other [3inv;2inv].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.inversion(3, 3)
        m.inversion(2, 2)
        assert_equal(str(m.mutated), str(Seq('AAGGATCG')))

    def test_adjecent_dup_dup_1(self):
        """
        Duplication and duplication directly adjecent to each other [2dup;3dup].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.duplication(2, 2)
        m.duplication(3, 3)
        assert_equal(str(m.mutated), str(Seq('ATTCCGATCG')))

    def test_adjecent_dup_dup_2(self):
        """
        Duplication and duplication directly adjecent to each other [3dup;2dup].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.duplication(3, 3)
        m.duplication(2, 2)
        assert_equal(str(m.mutated), str(Seq('ATTCCGATCG')))

    def test_adjecent_del_inv_1(self):
        """
        Deletion and inversion directly adjecent to each other [2del;3inv].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(2, 2)
        m.inversion(3, 3)
        assert_equal(str(m.mutated), str(Seq('AGGATCG')))

    def test_adjecent_del_inv_2(self):
        """
        Deletion and inversion directly adjecent to each other [3del;2inv].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.deletion(3, 3)
        m.inversion(2, 2)
        assert_equal(str(m.mutated), str(Seq('AAGATCG')))

    def test_adjecent_ins_inv_1(self):
        """
        Insertion and inversion directly adjecent to each other [2_3insG;3inv].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'G')
        m.inversion(3, 3)
        assert_equal(str(m.mutated), str(Seq('ATGGGATCG')))

    def test_adjecent_ins_inv_2(self):
        """
        Insertion and inversion directly adjecent to each other [2_3insG;2inv].
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insertion(2, 'G')
        m.inversion(2, 2)
        assert_equal(str(m.mutated), str(Seq('AAGCGATCG')))
