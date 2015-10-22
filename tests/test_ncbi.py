"""
Tests for the mutalyzer.ncbi module.
"""


from __future__ import unicode_literals

from mutalyzer import ncbi

from fixtures import with_links


@with_links(('NM_018650', 'NP_061120'))
def test_transcript_to_protein():
    """
    Get protein for transcript.
    """
    assert ncbi.transcript_to_protein('NM_018650') == 'NP_061120'


@with_links(('XM_005273133', None))
def test_transcript_to_protein_negative():
    """
    Get no protein for transcript.
    """
    assert ncbi.transcript_to_protein('XM_005273133') is None


@with_links(('NM_018650', 'NP_061120'))
def test_protein_to_transcript():
    """
    Get transcript for protein.
    """
    assert ncbi.protein_to_transcript('NP_061120') == 'NM_018650'
