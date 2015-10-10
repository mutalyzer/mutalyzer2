"""
Tests for the mutalyzer.db.queries module.
"""


from __future__ import unicode_literals

import pytest

from mutalyzer.db import queries


pytestmark = [
    pytest.mark.usefixtures('references'),
    pytest.mark.parametrize('references', [['MARK1']], indirect=True)
]


def test_get_transcript_protein_link():
    """
    Query a transcript-protein link by transcript.
    """
    link = queries.get_transcript_protein_link('NM_018650')
    assert link.transcript_accession == 'NM_018650'
    assert link.protein_accession == 'NP_061120'


def test_get_transcript_protein_link_negative():
    """
    Query a negative transcript-protein link by transcript.
    """
    link = queries.get_transcript_protein_link('XM_005273133')
    assert link.transcript_accession == 'XM_005273133'
    assert link.protein_accession is None


def test_get_transcript_protein_link_missing():
    """
    Query a missing transcript-protein link by transcript.
    """
    link = queries.get_transcript_protein_link('NM_123456')
    assert link is None


def test_get_transcript_protein_link_reverse():
    """
    Query a transcript-protein link by protein.
    """
    link = queries.get_transcript_protein_link('NP_061120', reverse=True)
    assert link.transcript_accession == 'NM_018650'
    assert link.protein_accession == 'NP_061120'


def test_get_transcript_protein_link_reverse_missing():
    """
    Query a missing transcript-protein link by protein.
    """
    link = queries.get_transcript_protein_link('NP_123456')
    assert link is None
