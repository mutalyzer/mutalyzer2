#!/usr/bin/env python

"""
Tests for the Mutalyzer module.
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
site.addsitedir('src')

from Modules import Config
from Modules import Output
import Mutalyzer


class TestMutalyzer(unittest.TestCase):
    """
    Test the Mutalyzer module.
    """

    def setUp(self):
        """
        Initialize test Mutalyzer module.
        """
        self.config = Config.Config()
        self.output = Output.Output(__file__, self.config.Output)

    def test_roll(self):
        """
        Just a variant where we should roll.
        """
        Mutalyzer.process('NM_003002.2:c.273del', self.config, self.output)
        wroll = self.output.getMessagesWithErrorCode('WROLL')
        self.assertTrue(len(wroll) > 0)

    def test_no_roll(self):
        """
        Just a variant where we cannot roll.
        """
        Mutalyzer.process('NM_003002.2:c.274del', self.config, self.output)
        wroll = self.output.getMessagesWithErrorCode('WROLL')
        self.assertTrue(len(wroll) == 0)

    def test_no_roll_splice(self):
        """
        Here we can roll but should not, because it is over a splice site.
        """
        Mutalyzer.process('NM_000088.3:g.459del', self.config, self.output)
        wrollback = self.output.getMessagesWithErrorCode('IROLLBACK')
        self.assertTrue(len(wrollback) > 0)
        wroll = self.output.getMessagesWithErrorCode('WROLL')
        self.assertTrue(len(wroll) == 0)

    def test_partial_roll_splice(self):
        """
        Here we can roll two positions, but should roll only one because
        otherwise it is over a splice site.
        """
        Mutalyzer.process('NM_000088.3:g.494del', self.config, self.output)
        wrollback = self.output.getMessagesWithErrorCode('IROLLBACK')
        self.assertTrue(len(wrollback) > 0)
        wroll = self.output.getMessagesWithErrorCode('WROLL')
        self.assertTrue(len(wroll) > 0)

    def test_roll_after_splice(self):
        """
        Here we can roll and should, we stay in the same exon.
        """
        Mutalyzer.process('NM_000088.3:g.460del', self.config, self.output)
        wroll = self.output.getMessagesWithErrorCode('WROLL')
        self.assertTrue(len(wroll) > 0)

    def test_ins_cds_start(self):
        """
        Insertion on CDS start boundary should not be included in CDS.
        """
        Mutalyzer.process('NM_000143.3:c.-1_1insCAT', self.config, self.output)
        self.assertEqual(self.output.getIndexedOutput("newprotein", 0), None)
        # Todo: is this a good test?

    def test_ins_cds_start_after(self):
        """
        Insertion after CDS start boundary should be included in CDS.
        """
        Mutalyzer.process('NM_000143.3:c.1_2insCAT', self.config, self.output)
        self.assertEqual(self.output.getIndexedOutput("newprotein", 0), '?')
        # Todo: is this a good test?

    def test_del_splice_site(self):
        """
        Deletion hitting one splice site should not be possible.
        """
        Mutalyzer.process('NG_012772.1(BRCA2_v001):c.632-5_670del',
                          self.config, self.output)
        woversplice = self.output.getMessagesWithErrorCode('WOVERSPLICE')
        self.assertTrue(len(woversplice) > 0)
        idelsplice = self.output.getMessagesWithErrorCode('IDELSPLICE')
        self.assertTrue(len(idelsplice) == 0)
        # Todo: For now, the following is how to check if no proteins
        # prediction is done.
        self.assertFalse(self.output.getOutput('newprotein'))

    def test_del_exon(self):
        """
        Deletion of an entire exon should be possible.
        """
        Mutalyzer.process('NG_012772.1(BRCA2_v001):c.632-5_681+7del',
                          self.config, self.output)
        woversplice = self.output.getMessagesWithErrorCode('WOVERSPLICE')
        self.assertTrue(len(woversplice) > 0)
        idelsplice = self.output.getMessagesWithErrorCode('IDELSPLICE')
        self.assertTrue(len(idelsplice) > 0)
        # Todo: For now, the following is how to check if no proteins
        # prediction is done.
        self.assertTrue(self.output.getOutput('newprotein'))

    def test_del_exon_in_frame(self):
        """
        Deletion of an entire exon with length a triplicate should give a
        proteine product with just this deletion (and possibly substitutions
        directly before and after).

        NG_012772.1(BRCA2_v001):c.68-7_316+7del is such a variant, since
        positions 68 through 316 are exactly one exon and (316-68+1)/3 = 83.
        """
        Mutalyzer.process('NG_012772.1(BRCA2_v001):c.68-7_316+7del',
                          self.config, self.output)
        woversplice = self.output.getMessagesWithErrorCode('WOVERSPLICE')
        self.assertTrue(len(woversplice) > 0)
        idelsplice = self.output.getMessagesWithErrorCode('IDELSPLICE')
        self.assertTrue(len(idelsplice) > 0)
        # Todo: For now, the following is how to check if no proteins
        # prediction is done.
        self.assertTrue(self.output.getOutput('newprotein'))
        # Todo: assert that protein products indeed have only this difference.

    def test_del_exons(self):
        """
        Deletion of two entire exons should be possible.
        """
        Mutalyzer.process('NG_012772.1(BRCA2_v001):c.632-5_793+7del',
                          self.config, self.output)
        woversplice = self.output.getMessagesWithErrorCode('WOVERSPLICE')
        self.assertTrue(len(woversplice) > 0)
        idelsplice = self.output.getMessagesWithErrorCode('IDELSPLICE')
        self.assertTrue(len(idelsplice) > 0)
        # Todo: For now, the following is how to check if no proteins
        # prediction is done.
        self.assertTrue(self.output.getOutput('newprotein'))

    def test_del_intron(self):
        """
        Deletion of an entire intron should be possible (fusion of remaining
        exonic parts).
        """
        Mutalyzer.process('NG_012772.1(BRCA2_v001):c.622_674del',
                          self.config, self.output)
        woversplice = self.output.getMessagesWithErrorCode('WOVERSPLICE')
        self.assertTrue(len(woversplice) > 0)
        idelsplice = self.output.getMessagesWithErrorCode('IDELSPLICE')
        self.assertTrue(len(idelsplice) > 0)
        # Todo: For now, the following is how to check if no proteins
        # prediction is done.
        self.assertTrue(self.output.getOutput('newprotein'))

    def test_del_intron_in_frame(self):
        """
        Deletion of an entire intron should be possible (fusion of remaining
        exonic parts).
        """
        Mutalyzer.process('NG_012772.1(BRCA2_v001):c.622_672del',
                          self.config, self.output)
        woversplice = self.output.getMessagesWithErrorCode('WOVERSPLICE')
        self.assertTrue(len(woversplice) > 0)
        idelsplice = self.output.getMessagesWithErrorCode('IDELSPLICE')
        self.assertTrue(len(idelsplice) > 0)
        # Todo: For now, the following is how to check if no proteins
        # prediction is done.
        self.assertTrue(self.output.getOutput('newprotein'))
        # Todo: assert that protein products indeed have only this difference.


if __name__ == '__main__':
    # Usage:
    #   ./test_mutalyzer.py -v
    # Or, selecting a specific test:
    #   ./test_mutalyzer.py -v TestMutalyzer.test_ins_cds_start
    unittest.main()
