"""
Tests for the mutalyzer.parsers.genbank module.
"""


from __future__ import unicode_literals

#import logging; logging.basicConfig()
import os

from mutalyzer.parsers import genbank
from mutalyzer.config import settings

from fixtures import REFERENCES
from fixtures import database, cache
from utils import MutalyzerTest
from utils import fix


class TestMutator(MutalyzerTest):
    """
    Test the mutator module.
    """
    fixtures = (database, )

    def setup(self):
        super(TestMutator, self).setup()
        self.gb_parser = genbank.GBparser()

    def test_product_lists_mismatch(self):
        """
        Test finding mismatches in some product lists.
        """
        tests = [(['a b c d e', 'a b C D e', 'a b c d e'], (2, 1)),
                 (['a b c d e', 'a b C d e', 'a B c d e'], (1, 2)),
                 (['a c d a', 'a b a', 'a a', 'a'], (1, 1)),
                 ([''], (-1, -1)),
                 (['', ''], (-1, -1)),
                 (['a', 'a'], (-1, -1)),
                 (['a', 'b'], (0, 0)),
                 (['a b c', 'a b c'], (-1, -1)),
                 (['a b c d a b', 'a b'], (2, 2))]
        for test in tests:
            assert self.gb_parser._find_mismatch(test[0]) == test[1]

    @fix(cache('A1BG'))
    def test_only_complete_genes_included(self):
        """
        Incomplete genes from the reference file should be ignored.
        """
        # contains A1BG (complete) and A1BG-AS1, ZNF497, LOC100419840
        # (incomplete).
        genbank_filename = os.path.join(settings.CACHE_DIR,
                                        REFERENCES['A1BG']['filename'])
        record = self.gb_parser.create_record(genbank_filename)
        assert [g.name for g in record.geneList] == ['A1BG']
