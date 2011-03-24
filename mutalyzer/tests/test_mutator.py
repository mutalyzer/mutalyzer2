#!/usr/bin/env python

"""
Tests for the mutalyzer.mutator module.
"""

#import logging; logging.basicConfig()
import re
import os
import random
import unittest
import site
from Bio.Seq import Seq

# Todo: Can this be done in a more elegant way?
os.chdir('../..')
site.addsitedir('.')

from mutalyzer import Config
from mutalyzer import Output
from mutalyzer import mutator


def _seq(length):
    """
    Return random DNA sequence of given length.
    """
    sequence = ''
    for i in range(length):
        sequence += random.choice('ACGT')
    return Seq(sequence)


class TestMutator(unittest.TestCase):
    """
    Test the mutator module.
    """

    def setUp(self):
        """
        Initialize test mutator module.
        """
        self.config = Config.Config()
        self.output = Output.Output(__file__, self.config.Output)

    def _mutator(self, sequence):
        """
        Create a Mutator instance for a given sequence.
        """
        return mutator.Mutator(sequence,
                               self.config.Mutator,
                               self.output)

    def test_shiftpos_no_change(self):
        """
        No change, no shifts.
        """
        l = 10
        m = self._mutator(_seq(l))
        # Numbering is 1-based
        for i in range(1, l + 1):
            self.assertEqual(m.shiftpos(i), i)

    def test_shiftpos_del_example(self):
        """
        Example of g.2del.
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.delM(2, 2)
        self.assertEqual(m.shiftpos(1), 1)
        self.assertEqual(m.shiftpos(2), 2)
        self.assertEqual(m.shiftpos(3), 2)

    def test_shiftpos_del(self):
        """
        Starting from the deleted position (not included), shift -1.
        """
        l = 10
        for d in range(1, l + 1):
            m = self._mutator(_seq(l))
            m.delM(d, d)
            for p in range(1, d + 1):
                self.assertEqual(m.shiftpos(p), p)
            for p in range(d + 1, l + 1):
                self.assertEqual(m.shiftpos(p), p - 1)

    def test_shiftpos_del2(self):
        """
        Starting from the deleted positions (not included), shift -2.
        """
        l = 10
        for d in range(1, l):
            m = self._mutator(_seq(l))
            m.delM(d, d + 1)
            for p in range(1, d + 2):
                self.assertEqual(m.shiftpos(p), p)
            for p in range(d + 2, l + 1):
                self.assertEqual(m.shiftpos(p), p - 2)

    def test_shiftpos_ins_example(self):
        """
        Example of g.2_3insA.
        """
        m = self._mutator(Seq('ATCGATCG'))
        m.insM(2, 'A')
        self.assertEqual(m.shiftpos(1), 1)
        self.assertEqual(m.shiftpos(2), 2)
        self.assertEqual(m.shiftpos(3), 4)

    def test_shiftpos_ins(self):
        """
        Starting from the interbase insertion position, shift +1.
        """
        l = 10
        for i in range(0, l + 1):
            m = self._mutator(_seq(l))
            m.insM(i, 'T')
            for p in range(1, i + 1):
                self.assertEqual(m.shiftpos(p), p)
            for p in range(i + 1, l + 1):
                self.assertEqual(m.shiftpos(p), p + 1)

    def test_shiftpos_ins2(self):
        """
        Starting from the interbase insertion position, shift +2.
        """
        l = 10
        for i in range(0, l + 1):
            m = self._mutator(_seq(l))
            m.insM(i, 'TT')
            for p in range(1, i + 1):
                self.assertEqual(m.shiftpos(p), p)
            for p in range(i + 1, l + 1):
                self.assertEqual(m.shiftpos(p), p + 2)

    def test_newSplice_no_change(self):
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
        self.assertEqual(m.newSplice(sites), sites)

    def test_newSplice_acc_del_before(self):
        """
        Deletion in intron directly before exon.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.delM(13, 13)   # g.13del
        self.assertEqual(m.newSplice(sites), [4, 9, 13, 16, 24, 26])

    def test_newSplice_acc_del_after(self):
        """
        Deletion at first exon position.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.delM(14, 14)   # g.14del
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 16, 24, 26])

    def test_newSplice_don_del_before(self):
        """
        Deletion at last exon position.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.delM(17, 17)   # g.17del
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 16, 24, 26])

    def test_newSplice_don_del_after(self):
        """
        Deletion in intron directly after exon.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.delM(18, 18)   # g.18del
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 17, 24, 26])

    def test_newSplice_acc_del2_before(self):
        """
        Deletion of 2 in intron directly before exon.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.delM(12, 13)   # g.12_13del
        self.assertEqual(m.newSplice(sites), [4, 9, 12, 15, 23, 25])

    def test_newSplice_acc_del2_on(self):
        """
        Deletion of 2 in intron/exon.

        @note: This hits a splice site, so we don't really support it.
        """
        return   # Disabled (see docstring)
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.delM(13, 14)   # g.13_14del
        self.assertEqual(m.newSplice(sites), [4, 9, 13, 15, 23, 25])

    def test_newSplice_acc_del2_after(self):
        """
        Deletion of 2 at first exon position.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.delM(14, 15)   # g.14_15del
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 15, 23, 25])

    def test_newSplice_don_del2_before(self):
        """
        Deletion of 2 at last exon positions.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.delM(16, 17)   # g.16_17del
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 15, 23, 25])

    def test_newSplice_don_del2_on(self):
        """
        Deletion of 2 in exon/intron.

        @note: This hits a splice site, so we don't really support it.
        """
        return   # Disabled (see docstring)
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.delM(17, 18)   # g.17_18del
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 16, 23, 25])

    def test_newSplice_don_del2_after(self):
        """
        Deletion of 2 in intron directly after exon.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.delM(18, 19)   # g.18_19del
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 17, 23, 25])

    def test_newSplice_acc_ins_before(self):
        """
        Insertion 1 position before intron/exon boundary.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(12, 'A')   # g.12_13insA
        self.assertEqual(m.newSplice(sites), [4, 9, 15, 18, 26, 28])

    def test_newSplice_acc_ins_on(self):
        """
        Insertion in intron/exon boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(13, 'A')   # g.13_14insA
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 18, 26, 28])

    def test_newSplice_first_acc_ins_on(self):
        """
        Insertion in first intron/exon boundary not be included.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(3, 'A')   # g.3_4insA
        self.assertEqual(m.newSplice(sites), [5, 10, 15, 18, 26, 28])

    def test_newSplice_acc_ins_after(self):
        """
        Insertion 1 position after intron/exon boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(14, 'A')   # g.14_15insA
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 18, 26, 28])

    def test_newSplice_don_ins_before(self):
        """
        Insertion 1 position before exon/intron boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(16, 'A')   # g.16_17insA
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 18, 26, 28])

    def test_newSplice_don_ins_on(self):
        """
        Insertion in exon/intron boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(17, 'A')   # g.17_18insA
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 18, 26, 28])

    def test_newSplice_last_don_ins_on(self):
        """
        Insertion in last exon/intron boundary should not be included.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(27, 'A')   # g.27_28insA
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 17, 25, 27])

    def test_newSplice_don_ins_after(self):
        """
        Insertion 1 position after exon/intron boundary.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(18, 'A')   # g.18_19insA
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 17, 26, 28])

    def test_newSplice_acc_ins2_before(self):
        """
        Insertion of 2 1 position before intron/exon boundary.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(12, 'AT')   # g.12_13insAT
        self.assertEqual(m.newSplice(sites), [4, 9, 16, 19, 27, 29])

    def test_newSplice_first_acc_ins2_on(self):
        """
        Insertion of 2 in last exon/intron boundary should not be included.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(3, 'AT')   # g.3_4insAT
        self.assertEqual(m.newSplice(sites), [6, 11, 16, 19, 27, 29])

    def test_newSplice_acc_ins2_after(self):
        """
        Insertion of 2 1 position after intron/exon boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(14, 'AT')   # g.14_15insAT
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 19, 27, 29])

    def test_newSplice_don_ins2_before(self):
        """
        Insertion of 2 1 position before exon/intron boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(16, 'AT')   # g.16_17insAT
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 19, 27, 29])

    def test_newSplice_last_don_ins2_on(self):
        """
        Insertion of 2 in last exon/intron boundary should not be included.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(27, 'AT')   # g.27_28insAT
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 17, 25, 27])

    def test_newSplice_don_ins2_after(self):
        """
        Insertion of 2 1 position after exon/intron boundary.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(18, 'AT')   # g.18_19insAT
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 17, 27, 29])

    def test_newSplice_acc_ins3_before(self):
        """
        Insertion of 3 1 position before intron/exon boundary.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(12, 'ATT')   # g.12_13insATT
        self.assertEqual(m.newSplice(sites), [4, 9, 17, 20, 28, 30])

    def test_newSplice_acc_ins3_on(self):
        """
        Insertion of 3 in intron/exon boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(13, 'ATT')   # g.13_14insATT
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 20, 28, 30])

    def test_newSplice_first_acc_ins3_on(self):
        """
        Insertion of 3 in first intron/exon boundary should not be included.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(3, 'ATT')   # g.3_4insATT
        self.assertEqual(m.newSplice(sites), [7, 12, 17, 20, 28, 30])

    def test_newSplice_acc_ins3_after(self):
        """
        Insertion of 3 1 position after intron/exon boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(14, 'ATT')   # g.14_15insATT
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 20, 28, 30])

    def test_newSplice_don_ins3_before(self):
        """
        Insertion of 3 1 position before exon/intron boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(16, 'ATT')   # g.16_17insATT
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 20, 28, 30])

    def test_newSplice_don_ins3_on(self):
        """
        Insertion of 3 in exon/intron boundary.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(17, 'ATT')   # g.17_18insATT
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 20, 28, 30])

    def test_newSplice_last_don_ins3_on(self):
        """
        Insertion of 3 in last exon/intron boundary should not be included.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(27, 'ATT')   # g.27_28insATT
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 17, 25, 27])

    def test_newSplice_don_ins3_after(self):
        """
        Insertion of 3 1 position after exon/intron boundary.

        @note: This hits a splice site, so we don't really support it.
        """
        l = 30
        sites = [4, 9, 14, 17, 25, 27]
        m = self._mutator(_seq(l))
        m.insM(18, 'ATT')   # g.18_19insATT
        self.assertEqual(m.newSplice(sites), [4, 9, 14, 17, 28, 30])

    def test_newSplice_adj_del_before1(self):
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
        m.delM(16, 16)   # g.16del
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 16, 17, 26])

    def test_newSplice_adj_del_before(self):
        """
        Adjacent exons: deletion at last position of first exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.delM(17, 17)   # g.17del
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 16, 17, 26])

    def test_newSplice_adj_del_after(self):
        """
        Adjacent exons: deletion at first position of second exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.delM(18, 18)   # g.18del
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 17, 18, 26])

    def test_newSplice_adj_del_after1(self):
        """
        Adjacent exons: deletion at second position of second exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.delM(19, 19)   # g.19del
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 17, 18, 26])

    def test_newSplice_adj_ins_before(self):
        """
        Adjacent exons: insertion 1 position before exon/exon boundary.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.insM(16, 'A')   # g.16_17insA
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 18, 19, 28])

    def test_newSplice_adj_ins_on(self):
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
        m.insM(17, 'A')   # g.17_18insA
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 18, 19, 28])

    def test_newSplice_adj_ins_after(self):
        """
        Adjacent exons: insertion 1 position after exon/exon boundary.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.insM(18, 'A')   # g.18_19insA
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 17, 18, 28])

    def test_newSplice_adj_del2_before1(self):
        """
        Adjacent exons: deletion of 2 at second-last position of first exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.delM(15, 16)   # g.15_16del
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 15, 16, 25])

    def test_newSplice_adj_del2_before(self):
        """
        Adjacent exons: deletion of 2 at last position of first exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.delM(16, 17)   # g.16_17del
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 15, 16, 25])

    def test_newSplice_adj_del2_on(self):
        """
        Adjacent exons: deletion of 2 at exon/exon boundary.

        @todo: This is a special case of bug #????. Once fixed, the two
               exons will be joined to one new exon.
        """
        return   # Disabled (see docstring)
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.delM(17, 18)   # g.17_18del
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 16, 17, 25])

    def test_newSplice_adj_del2_after(self):
        """
        Adjacent exons: deletion of 2 at first position of second exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.delM(18, 19)   # g.18_19del
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 17, 18, 25])

    def test_newSplice_adj_del2_after1(self):
        """
        Adjacent exons: deletion of 2 at second position of second exon.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.delM(19, 20)   # g.19_20del
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 17, 18, 25])

    def test_newSplice_adj_ins2_before(self):
        """
        Adjacent exons: insertion of 2 1 position before exon/exon boundary.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.insM(16, 'AT')   # g.16_17insAT
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 19, 20, 29])

    def test_newSplice_adj_ins2_on(self):
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
        m.insM(17, 'AT')   # g.17_18insAT
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 19, 20, 29])

    def test_newSplice_adj_ins2_after(self):
        """
        Adjacent exons: insertion of 2 1 position after exon/exon boundary.
        """
        l = 30
        sites = [4, 9, 10, 17, 18, 27]
        m = self._mutator(_seq(l))
        m.insM(18, 'AT')   # g.18_19insAT
        self.assertEqual(m.newSplice(sites), [4, 9, 10, 17, 18, 29])


if __name__ == '__main__':
    # Usage:
    #   ./test_mutator.py -v
    # Or, selecting a specific test:
    #   ./test_mutator.py -v TestMutator.test_mutated
    unittest.main()
