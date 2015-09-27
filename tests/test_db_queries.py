"""
Tests for the mutalyzer.db.queries module.
"""


from __future__ import unicode_literals

#import logging; logging.basicConfig()

from mutalyzer.db import queries

from fixtures import database, cache
from utils import MutalyzerTest
from utils import fix


class TestMutator(MutalyzerTest):
    """
    Test the queries module.
    """
    fixtures = (database, )

    def setup(self):
        super(TestMutator, self).setup()

    @fix(cache('MARK1'))
    def test_get_transcript_protein_link(self):
        """
        Query a transcript-protein link by transcript.
        """
        link = queries.get_transcript_protein_link('NM_018650')
        assert link.transcript_accession == 'NM_018650'
        assert link.protein_accession == 'NP_061120'

    @fix(cache('MARK1'))
    def test_get_transcript_protein_link_negative(self):
        """
        Query a negative transcript-protein link by transcript.
        """
        link = queries.get_transcript_protein_link('XM_005273133')
        assert link.transcript_accession == 'XM_005273133'
        assert link.protein_accession is None

    @fix(cache('MARK1'))
    def test_get_transcript_protein_link_missing(self):
        """
        Query a missing transcript-protein link by transcript.
        """
        link = queries.get_transcript_protein_link('NM_123456')
        assert link is None

    @fix(cache('MARK1'))
    def test_get_transcript_protein_link_reverse(self):
        """
        Query a transcript-protein link by protein.
        """
        link = queries.get_transcript_protein_link('NP_061120', reverse=True)
        assert link.transcript_accession == 'NM_018650'
        assert link.protein_accession == 'NP_061120'

    @fix(cache('MARK1'))
    def test_get_transcript_protein_link_reverse_missing(self):
        """
        Query a missing transcript-protein link by protein.
        """
        link = queries.get_transcript_protein_link('NP_123456')
        assert link is None
