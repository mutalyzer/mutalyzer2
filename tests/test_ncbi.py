"""
Tests for the mutalyzer.ncbi module.
"""


from __future__ import unicode_literals

import Bio.Entrez
import pytest

from mutalyzer import ncbi
from mutalyzer.redisclient import client as redis

from fixtures import with_links


@pytest.fixture
def entrez(request, monkeypatch):
    """
    Fixture monkey-patching the NCBI Entrez API to return transcript-protein
    links defined in the fixture parameter.

    The fixture is similar to the :func:`fixtures.links` fixture, but instead
    of storing the links in the cache, the API is monkey-patched.
    """
    try:
        links = request.param
    except AttributeError:
        return []

    # We need two-way lookup.
    transcript_to_protein = dict(links)
    protein_to_transcript = dict((p, t) for t, p in links)

    # Store original methods which should be called as a fallback.
    esearch = Bio.Entrez.esearch
    elink = Bio.Entrez.elink
    efetch = Bio.Entrez.efetch

    # Intermediate Entrez result object which can be parsed with Entrez.read.
    class EntrezResult(object):
        def __init__(self, result):
            self.result = result

        def read(self):
            return self.result

        def close(self):
            pass

    def mock_esearch(db=None, term=None):
        if ((db == 'nucleotide' and term in transcript_to_protein)
                or (db == 'protein' and term in protein_to_transcript)):
            return EntrezResult({
                'Count': '1',
                'RetMax': '1',
                'IdList': [term],
                'TranslationSet': [],
                'RetStart': '0',
                'QueryTranslation': ''
            })
        return esearch(db=db, term=term)

    def mock_elink(dbfrom=None, db=None, id=None):
        if dbfrom == 'nucleotide' and id in transcript_to_protein:
            if transcript_to_protein[id] is None:
                linkset = []
            else:
                linkset = [{'DbTo': 'protein',
                            'Link': [{'Id': transcript_to_protein[id]}],
                            'LinkName': 'nuccore_protein'}]
            return EntrezResult([{
                'LinkSetDb': linkset,
                'DbFrom': 'nuccore',
                'IdList': [id],
                'LinkSetDbHistory': [],
                'ERROR': []
            }])
        if dbfrom == 'protein' and id in protein_to_transcript:
            if protein_to_transcript[id] is None:
                linkset = []
            else:
                linkset = [{'DbTo': 'nuccore',
                            'Link': [{'Id': '568815587'},
                                     {'Id': '528476600'},
                                     {'Id': '568815270'},
                                     {'Id': '528474155'},
                                     {'Id': '452415518'},
                                     {'Id': '452405284'},
                                     {'Id': '383209650'}],
                            'LinkName': 'protein_nuccore'},
                           {'DbTo': 'nuccore',
                            'Link': [{'Id': '4506864'}],
                            'LinkName': 'protein_nuccore_cds'},
                           {'DbTo': 'nuccore',
                            'Link': [{'Id': '48735311'},
                                     {'Id': '48734961'},
                                     {'Id': '47682402'},
                                     {'Id': '18490203'},
                                     {'Id': '16359050'},
                                     {'Id': '16306997'},
                                     {'Id': '15929518'},
                                     {'Id': '15214938'},
                                     {'Id': '13528941'}],
                            'LinkName': 'protein_nuccore_mgc_refseq'},
                           {'DbTo': 'nuccore',
                            'Link': [{'Id': protein_to_transcript[id]}],
                            'LinkName': 'protein_nuccore_mrna'}]
            return EntrezResult([{
                'LinkSetDb': linkset,
                'DbFrom': 'protein',
                'IdList': [id],
                'LinkSetDbHistory': [],
                'ERROR': []
            }])
        return elink(dbfrom=dbfrom, db=db, id=id)

    def mock_efetch(db=None, id=None, rettype=None, retmode=None):
        if ((db == 'nucleotide' and id in transcript_to_protein)
                or (db == 'protein' and id in protein_to_transcript)):
            if '.' not in id:
                id += '.9999'
            return EntrezResult(id + '\n')
        return efetch(db=db, id=id, rettype=rettype, retmode=retmode)

    def mock_read(result):
        return result.read()

    monkeypatch.setattr(Bio.Entrez, 'esearch', mock_esearch)
    monkeypatch.setattr(Bio.Entrez, 'elink', mock_elink)
    monkeypatch.setattr(Bio.Entrez, 'efetch', mock_efetch)
    monkeypatch.setattr(Bio.Entrez, 'read', mock_read)
    return links


def with_entrez(*links):
    """
    Convenience decorator for parameterizing tests with transcript-protein
    link fixtures in the Entrez API.

    Similar to :func:`fixtures.with_links`.
    """
    def test_with_entrez(test):
        return pytest.mark.usefixtures('entrez')(
            pytest.mark.parametrize(
                'entrez', [links], indirect=True,
                ids=[','.join('/'.join(a or '*' for a in l)
                              for l in links)])(test))
    return test_with_entrez


@with_entrez(('NM_11111.1', None),
             ('NM_11111.2', 'NP_11111.2'),
             ('NM_22222.2', None),
             ('NM_22222.3', 'NP_22222.3'),
             ('NM_33333.4', None),
             ('NM_33333.5', 'NP_33333.5'),
             ('NM_44444', None),
             ('NM_44444.5', None),
             ('NM_44444.6', None),
             ('NM_55555', 'NP_55555'),
             ('NM_55555.6', None),
             ('NM_66666', 'NP_66666'),
             ('NM_66666.6', 'NP_66666.6'),
             ('NM_66666.7', 'NP_66666.7'),
             ('NM_66666.8', None),
             ('NM_77777', 'NP_77777'),
             ('NM_77777.7', 'NP_77777.7'),
             ('NM_77777.8', None),
             ('NM_88888', None),
             ('NM_88888.8', None),
             ('NM_88888.9', 'NP_88888.9'))
@with_links(('NM_11111', 'NP_11111'),
            ('NM_22222', None),
            ('NM_33333.3', 'NP_33333.3'),
            ('NM_44444.4', None),
            ('NM_55555.5', None),
            ('NM_66666.6', None))
@pytest.mark.parametrize('accession,version,match_version,expected', [
    ('NM_11111', None, False, ('NP_11111', None)),
    ('NM_11111', 1, False, ('NP_11111', None)),
    ('NM_11111', 1, True, None),
    ('NM_11111', 2, False, ('NP_11111', None)),
    ('NM_11111', 2, True, ('NP_11111', 2)),
    ('NM_22222', None, False, None),
    ('NM_22222', 2, False, None),
    ('NM_22222', 2, True, None),
    ('NM_22222', 3, False, None),
    ('NM_22222', 3, True, ('NP_22222', 3)),
    ('NM_33333', None, False, ('NP_33333', None)),
    ('NM_33333', 3, True, ('NP_33333', 3)),
    ('NM_33333', 3, False, ('NP_33333', 3)),
    ('NM_33333', 4, True, None),
    ('NM_33333', 4, False, ('NP_33333', None)),
    ('NM_33333', 5, True, ('NP_33333', 5)),
    ('NM_33333', 5, False, ('NP_33333', None)),
    ('NM_44444', None, False, None),
    ('NM_44444', 4, True, None),
    ('NM_44444', 4, False, None),
    ('NM_44444', 5, True, None),
    ('NM_44444', 5, False, None),
    ('NM_44444', 6, True, None),
    ('NM_44444', 6, False, None),
    ('NM_55555', None, False, ('NP_55555', None)),
    ('NM_55555', 5, True, None),
    ('NM_55555', 5, False, None),
    ('NM_55555', 6, True, None),
    ('NM_55555', 6, False, ('NP_55555', None)),
    ('NM_66666', None, False, ('NP_66666', None)),
    ('NM_66666', 6, True, None),
    ('NM_66666', 6, False, None),
    ('NM_66666', 7, True, ('NP_66666', 7)),
    ('NM_66666', 7, False, ('NP_66666', 7)),
    ('NM_66666', 8, True, None),
    ('NM_66666', 8, False, ('NP_66666', None)),
    ('NM_77777', None, False, ('NP_77777', None)),
    ('NM_77777', 7, False, ('NP_77777', 7)),
    ('NM_77777', 7, True, ('NP_77777', 7)),
    ('NM_77777', 8, False, ('NP_77777', None)),
    ('NM_77777', 8, True, None),
    ('NM_88888', None, False, None),
    ('NM_88888', 8, False, None),
    ('NM_88888', 8, True, None),
    ('NM_88888', 9, False, ('NP_88888', 9)),
    ('NM_88888', 9, True, ('NP_88888', 9))])
def test_transcript_to_protein(accession, version, match_version, expected):
    """
    Get protein for transcript.

    Both the Entrez API and our cache are fixed with a set of
    transcript-protein links. This test is parametrized with a list of
    arguments for the :func:`ncbi.transcript_to_protein` function and the
    corresponding expected result.
    """
    assert ncbi.transcript_to_protein(
        accession, version, match_version) == expected


@with_entrez((None, 'NP_11111.1'),
             ('NM_11111.2', 'NP_11111.2'),
             (None, 'NP_22222.2'),
             ('NM_22222.3', 'NP_22222.3'),
             (None, 'NP_33333.4'),
             ('NM_33333.5', 'NP_33333.5'),
             (None, 'NP_44444'),
             (None, 'NP_44444.5'),
             (None, 'NP_44444.6'),
             ('NM_55555', 'NP_55555'),
             (None, 'NP_55555.6'),
             ('NM_66666', 'NP_66666'),
             ('NM_66666.6', 'NP_66666.6'),
             ('NM_66666.7', 'NP_66666.7'),
             (None, 'NP_66666.8'),
             ('NM_77777', 'NP_77777'),
             ('NM_77777.7', 'NP_77777.7'),
             (None, 'NP_77777.8'),
             (None, 'NP_88888'),
             (None, 'NP_88888.8'),
             ('NM_88888.9', 'NP_88888.9'))
@with_links(('NM_11111', 'NP_11111'),
            (None, 'NP_22222'),
            ('NM_33333.3', 'NP_33333.3'),
            (None, 'NP_44444.4'),
            (None, 'NP_55555.5'),
            (None, 'NP_66666.6'))
@pytest.mark.parametrize('accession,version,match_version,expected', [
    ('NP_11111', None, False, ('NM_11111', None)),
    ('NP_11111', 1, False, ('NM_11111', None)),
    ('NP_11111', 1, True, None),
    ('NP_11111', 2, False, ('NM_11111', None)),
    ('NP_11111', 2, True, ('NM_11111', 2)),
    ('NP_22222', None, False, None),
    ('NP_22222', 2, False, None),
    ('NP_22222', 2, True, None),
    ('NP_22222', 3, False, None),
    ('NP_22222', 3, True, ('NM_22222', 3)),
    ('NP_33333', None, False, ('NM_33333', None)),
    ('NP_33333', 3, True, ('NM_33333', 3)),
    ('NP_33333', 3, False, ('NM_33333', 3)),
    ('NP_33333', 4, True, None),
    ('NP_33333', 4, False, ('NM_33333', None)),
    ('NP_33333', 5, True, ('NM_33333', 5)),
    ('NP_33333', 5, False, ('NM_33333', None)),
    ('NP_44444', None, False, None),
    ('NP_44444', 4, True, None),
    ('NP_44444', 4, False, None),
    ('NP_44444', 5, True, None),
    ('NP_44444', 5, False, None),
    ('NP_44444', 6, True, None),
    ('NP_44444', 6, False, None),
    ('NP_55555', None, False, ('NM_55555', None)),
    ('NP_55555', 5, True, None),
    ('NP_55555', 5, False, None),
    ('NP_55555', 6, True, None),
    ('NP_55555', 6, False, ('NM_55555', None)),
    ('NP_66666', None, False, ('NM_66666', None)),
    ('NP_66666', 6, True, None),
    ('NP_66666', 6, False, None),
    ('NP_66666', 7, True, ('NM_66666', 7)),
    ('NP_66666', 7, False, ('NM_66666', 7)),
    ('NP_66666', 8, True, None),
    ('NP_66666', 8, False, ('NM_66666', None)),
    ('NP_77777', None, False, ('NM_77777', None)),
    ('NP_77777', 7, False, ('NM_77777', 7)),
    ('NP_77777', 7, True, ('NM_77777', 7)),
    ('NP_77777', 8, False, ('NM_77777', None)),
    ('NP_77777', 8, True, None),
    ('NP_88888', None, False, None),
    ('NP_88888', 8, False, None),
    ('NP_88888', 8, True, None),
    ('NP_88888', 9, False, ('NM_88888', 9)),
    ('NP_88888', 9, True, ('NM_88888', 9))])
def test_protein_to_transcript(accession, version, match_version, expected):
    """
    Get transcript for protein.

    Both the Entrez API and our cache are fixed with a set of
    transcript-protein links. This test is parametrized with a list of
    arguments for the :func:`ncbi.transcript_to_protein` function and the
    corresponding expected result.

    Fixtures and parameters of this test mirror those of the
    `test_transcript_to_protein` test.
    """
    assert ncbi.protein_to_transcript(
        accession, version, match_version) == expected


@with_entrez(('NM_11111', None),
             ('NM_22222', 'NP_22222'),
             ('NM_33333', None),
             ('NM_33333.3', None),
             ('NM_44444', None),
             ('NM_44444.4', 'NP_44444.4'))
@pytest.mark.parametrize('accession,version,match_version,expected_forward,expected_reverse', [
    ('NM_11111', None, False, [('NM_11111', None)], []),
    ('NM_22222', None, False,
     [('NM_22222', 'NP_22222')], [('NM_22222', 'NP_22222')]),
    ('NM_33333', None, False, [('NM_33333', None)], []),
    ('NM_33333', 3, False, [('NM_33333', None), ('NM_33333.3', None)], []),
    ('NM_33333', 3, True, [('NM_33333.3', None)], []),
    ('NM_44444', None, False, [('NM_44444', None)], []),
    ('NM_44444', 4, False,
     [('NM_44444', 'NP_44444'), ('NM_44444.4', 'NP_44444.4')],
     [('NM_44444', 'NP_44444'), ('NM_44444.4', 'NP_44444.4')]),
    ('NM_44444', 4, True,
     [('NM_44444', 'NP_44444'), ('NM_44444.4', 'NP_44444.4')],
     [('NM_44444', 'NP_44444'), ('NM_44444.4', 'NP_44444.4')])])
def test_transcript_to_protein_cache(accession, version, match_version,
                                     expected_forward, expected_reverse):
    """
    Get protein for transcript and check the resulting cache state.
    """
    ncbi.transcript_to_protein(accession, version, match_version)

    forward = [(key.split(':')[-1], redis.get(key) or None)
               for key in redis.keys('ncbi:transcript-to-protein:*')]
    assert sorted(forward) == sorted(expected_forward)

    reverse = [(redis.get(key) or None, key.split(':')[-1])
               for key in redis.keys('ncbi:protein-to-transcript:*')]
    assert sorted(reverse) == sorted(expected_reverse)


@with_entrez((None, 'NP_11111'),
             ('NM_22222', 'NP_22222'),
             (None, 'NP_33333'),
             (None, 'NP_33333.3'),
             (None, 'NP_44444'),
             ('NM_44444.4', 'NP_44444.4'))
@pytest.mark.parametrize('accession,version,match_version,expected_forward,expected_reverse', [
    ('NP_11111', None, False, [], [(None, 'NP_11111')]),
    ('NP_22222', None, False,
     [('NM_22222', 'NP_22222')], [('NM_22222', 'NP_22222')]),
    ('NP_33333', None, False, [], [(None, 'NP_33333')]),
    ('NP_33333', 3, False, [], [(None, 'NP_33333'), (None, 'NP_33333.3')]),
    ('NP_33333', 3, True, [], [(None, 'NP_33333.3')]),
    ('NP_44444', None, False, [], [(None, 'NP_44444')]),
    ('NP_44444', 4, False,
     [('NM_44444', 'NP_44444'), ('NM_44444.4', 'NP_44444.4')],
     [('NM_44444', 'NP_44444'), ('NM_44444.4', 'NP_44444.4')]),
    ('NP_44444', 4, True,
     [('NM_44444', 'NP_44444'), ('NM_44444.4', 'NP_44444.4')],
     [('NM_44444', 'NP_44444'), ('NM_44444.4', 'NP_44444.4')])])
def test_protein_to_transcript_cache(accession, version, match_version,
                                     expected_forward, expected_reverse):
    """
    Get transcript for protein and check the resulting cache state.
    """
    ncbi.protein_to_transcript(accession, version, match_version)

    forward = [(key.split(':')[-1], redis.get(key) or None)
               for key in redis.keys('ncbi:transcript-to-protein:*')]
    assert sorted(forward) == sorted(expected_forward)

    reverse = [(redis.get(key) or None, key.split(':')[-1])
               for key in redis.keys('ncbi:protein-to-transcript:*')]
    assert sorted(reverse) == sorted(expected_reverse)
