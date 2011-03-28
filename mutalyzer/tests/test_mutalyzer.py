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
site.addsitedir('.')

from mutalyzer.config import Config
from mutalyzer.output import Output
from mutalyzer.variantchecker import check_variant


class TestMutalyzer(unittest.TestCase):
    """
    Test the Mutalyzer module.
    """

    def setUp(self):
        """
        Initialize test Mutalyzer module.
        """
        self.config = Config()
        self.output = Output(__file__, self.config.Output)

    def test_roll(self):
        """
        Just a variant where we should roll.
        """
        check_variant('NM_003002.2:c.273del', self.config, self.output)
        wroll = self.output.getMessagesWithErrorCode('WROLL')
        self.assertTrue(len(wroll) > 0)

    def test_no_roll(self):
        """
        Just a variant where we cannot roll.
        """
        check_variant('NM_003002.2:c.274del', self.config, self.output)
        wroll = self.output.getMessagesWithErrorCode('WROLL')
        self.assertTrue(len(wroll) == 0)

    def test_no_roll_splice(self):
        """
        Here we can roll but should not, because it is over a splice site.
        """
        check_variant('NM_000088.3:g.459del', self.config, self.output)
        wrollback = self.output.getMessagesWithErrorCode('IROLLBACK')
        self.assertTrue(len(wrollback) > 0)
        wroll = self.output.getMessagesWithErrorCode('WROLL')
        self.assertTrue(len(wroll) == 0)

    def test_partial_roll_splice(self):
        """
        Here we can roll two positions, but should roll only one because
        otherwise it is over a splice site.
        """
        check_variant('NM_000088.3:g.494del', self.config, self.output)
        wrollback = self.output.getMessagesWithErrorCode('IROLLBACK')
        self.assertTrue(len(wrollback) > 0)
        wroll = self.output.getMessagesWithErrorCode('WROLL')
        self.assertTrue(len(wroll) > 0)

    def test_roll_after_splice(self):
        """
        Here we can roll and should, we stay in the same exon.
        """
        check_variant('NM_000088.3:g.460del', self.config, self.output)
        wroll = self.output.getMessagesWithErrorCode('WROLL')
        self.assertTrue(len(wroll) > 0)

    def test_ins_cds_start(self):
        """
        Insertion on CDS start boundary should not be included in CDS.
        """
        check_variant('NM_000143.3:c.-1_1insCAT', self.config, self.output)
        self.assertEqual(self.output.getIndexedOutput("newprotein", 0), None)

    def test_ins_cds_start_after(self):
        """
        Insertion after CDS start boundary should be included in CDS.
        """
        check_variant('NM_000143.3:c.1_2insCAT', self.config, self.output)
        self.assertEqual(self.output.getIndexedOutput("newprotein", 0), '?')


if __name__ == '__main__':
    # Usage:
    #   ./test_mutalyzer.py -v
    # Or, selecting a specific test:
    #   ./test_mutalyzer.py -v TestMutalyzer.test_ins_cds_start
    unittest.main()
