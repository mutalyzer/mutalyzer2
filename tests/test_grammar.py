"""
Tests for the mutalyzer.grammar module.
"""


#import logging; logging.basicConfig()
import os
from nose.tools import *

import mutalyzer
from mutalyzer.grammar import Grammar
from mutalyzer.output import Output


class TestGrammar():
    """
    Test the mytalyzer.grammar module.
    """

    def setUp(self):
        """
        Initialize test Grammar instance.
        """
        self.output = Output(__file__)
        self.grammar = Grammar(self.output)

    def test_some_variants(self):
        """
        Some example variants.
        """
        self.grammar.parse('NM_002001.2:c.[12del]')
        self.grammar.parse('NM_002001.2:c.[(12del)]')
        self.grammar.parse('NM_002001.2:c.[(12del)?]')
        self.grammar.parse('NM_002001.2:c.[(12del);(12del)]')
        self.grammar.parse('NM_002001.2:c.[(12del;12del)]')
        self.grammar.parse('NM_002001.2:c.[((12del)?;12del)?]')
