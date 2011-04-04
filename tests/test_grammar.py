"""
Tests for the mutalyzer.grammar module.
"""

#import logging; logging.basicConfig()
import os
import site
from nose.tools import *

# Todo: Get this from the configuration file
root_dir = os.path.split(os.path.dirname(__file__))[0]
site.addsitedir(root_dir)
# Todo: Fix Mutalyzer to not depend on working directory
if not __name__ == '__main__':
    os.chdir(root_dir)

from mutalyzer.config import Config
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
        self.config = Config()
        self.output = Output(__file__, self.config.Output)
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
