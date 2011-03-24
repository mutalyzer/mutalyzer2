#!/usr/bin/env python

"""
Tests for the mutalyzer.grammar module.
"""

#import logging; logging.basicConfig()
import os
import unittest
import site

# Todo: Can this be done in a more elegant way?
os.chdir('../..')
site.addsitedir('.')

from mutalyzer.grammar import Grammar
from mutalyzer import Config
from mutalyzer import Output


class TestGrammar(unittest.TestCase):
    """
    Test the mytalyzer.grammar module.
    """

    def setUp(self):
        """
        Initialize test Grammar instance.
        """
        self.config = Config.Config()
        self.output = Output.Output(__file__, self.config.Output)
        self.grammar = Grammar(self.output)

    def test_some_variants(self):
        """
        No change, no shifts.
        """
        self.grammar.parse('NM_002001.2:c.[12del]')
        self.grammar.parse('NM_002001.2:c.[(12del)]')
        self.grammar.parse('NM_002001.2:c.[(12del)?]')
        self.grammar.parse('NM_002001.2:c.[(12del);(12del)]')
        self.grammar.parse('NM_002001.2:c.[(12del;12del)]')
        self.grammar.parse('NM_002001.2:c.[((12del)?;12del)?]')


if __name__ == '__main__':
    # Usage:
    #   ./test_parser.py -v
    # Or, selecting a specific test:
    #   ./test_parser.py -v TestParser.test_some_variants
    unittest.main()
