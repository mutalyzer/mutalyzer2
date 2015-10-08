"""
Tests for the mutalyzer.backtranslator module.
"""


from __future__ import unicode_literals

import pytest

from mutalyzer.backtranslator import backtranslate
import mutalyzer.Retriever

from fixtures import with_links, with_references


DIVERSIONS = {
    'NM_003002': 'NM_003002.2'
}


@pytest.fixture
def divert(monkeypatch, request):
    """
    Ad-hoc fixture diverting a retriever call for NM_003002 to NM_003002.2.

    This is useful to mimic behaviour of the retriever for the accession
    without version number, which we cannot do with a normal fixture (the
    retriever does not support non-versioned references in the database).
    """
    loadrecord = mutalyzer.Retriever.GenBankRetriever.loadrecord
    def mock_loadrecord(self, identifier):
        identifier = DIVERSIONS.get(identifier, identifier)
        return loadrecord(self, identifier)
    monkeypatch.setattr(mutalyzer.Retriever.GenBankRetriever, 'loadrecord', mock_loadrecord)


@with_references('NM_003002.2')
def test_nm(output):
    """
    Back translate on NM.
    """
    variants = backtranslate(output, 'NM_003002.2:p.Asp92Tyr')
    assert sorted(variants) == ['NM_003002.2:c.274G>T']
    assert len(output.getMessagesWithErrorCode('WNOVERSION')) == 0
    assert len(output.getMessagesWithErrorCode('WNODNA')) == 0
    assert len(output.getMessagesWithErrorCode('WIMPROVE')) == 0


@with_references('NM_003002.2')
@pytest.mark.usefixtures('divert')
def test_nm_no_version(output):
    """
    Back translate on NM.
    """
    variants = backtranslate(output, 'NM_003002:p.Asp92Tyr')
    assert sorted(variants) == ['NM_003002:c.274G>T']
    assert len(output.getMessagesWithErrorCode('WNOVERSION')) == 0
    assert len(output.getMessagesWithErrorCode('WNODNA')) == 0
    assert len(output.getMessagesWithErrorCode('WIMPROVE')) == 0


@with_links(('NM_003002.2', 'NP_002993.1'))
@with_references('NM_003002.2')
def test_np(output):
    """
    Back translate on NP.
    """
    variants = backtranslate(output, 'NP_002993.1:p.Asp92Glu')
    assert sorted(variants) == ['NM_003002.2:c.276C>A',
                                'NM_003002.2:c.276C>G']
    assert len(output.getMessagesWithErrorCode('WNOVERSION')) == 0
    assert len(output.getMessagesWithErrorCode('WNODNA')) == 0
    assert len(output.getMessagesWithErrorCode('WIMPROVE')) == 0


@with_links(('NM_003002', 'NP_002993'))
@with_references('NM_003002.2')
@pytest.mark.usefixtures('divert')
def test_np_no_version(output):
    """
    Back translate on NP.
    """
    variants = backtranslate(output, 'NP_002993:p.Asp92Glu')
    assert sorted(variants) == ['NM_003002:c.276C>A',
                                'NM_003002:c.276C>G']
    assert len(output.getMessagesWithErrorCode('WNOVERSION')) == 1
    assert len(output.getMessagesWithErrorCode('WNODNA')) == 0
    assert len(output.getMessagesWithErrorCode('WIMPROVE')) == 0


@with_links((None, 'NP_000000.0'))
def test_np_no_nm(output):
    """
    Back translate on NP without NM.
    """
    variants = backtranslate(output, 'NP_000000.0:p.Asp92Tyr')
    assert sorted(variants) == ['UNKNOWN:c.274G>T']
    assert len(output.getMessagesWithErrorCode('WNOVERSION')) == 0
    assert len(output.getMessagesWithErrorCode('WNODNA')) == 1
    assert len(output.getMessagesWithErrorCode('WIMPROVE')) == 0


@with_links((None, 'NP_000000.0'))
def test_np_no_nm_improvable(output):
    """
    Back translate on NP without NM improvable.
    """
    variants = backtranslate(output, 'NP_000000.0:p.Leu92Phe')
    assert sorted(variants) == ['UNKNOWN:c.274C>T',
                                'UNKNOWN:c.276A>C',
                                'UNKNOWN:c.276A>T',
                                'UNKNOWN:c.276G>C',
                                'UNKNOWN:c.276G>T']
    assert len(output.getMessagesWithErrorCode('WNOVERSION')) == 0
    assert len(output.getMessagesWithErrorCode('WNODNA')) == 1
    assert len(output.getMessagesWithErrorCode('WIMPROVE')) == 1
