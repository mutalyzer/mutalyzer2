"""
Tests for the mutalyzer.services.json module.
"""


from __future__ import unicode_literals

import json
import pytest
from spyne.model.fault import Fault
from spyne.server.null import NullServer

import mutalyzer
from mutalyzer import announce
from mutalyzer.config import settings
from mutalyzer import Scheduler
from mutalyzer.services.json import application


# Todo: We currently have no way of testing POST requests to the JSON API. We
#     had some tests for this, but they were removed with the new setup [1].
#     They depended on a patch to Spyne [2].
#
# [1] https://git.lumc.nl/mutalyzer/mutalyzer2/commit/fbf42d8c4a6c6ec27f9001686941b835867199ff
# [2] https://github.com/LUMC/spyne/commit/58660dec28d47b1c3bf1e46d20f55a913ad036cd


@pytest.fixture
def server():
    return NullServer(application, ostr=True)


@pytest.fixture
def api(server):
    def call(method, *args, **kwargs):
        r = getattr(server.service, method)(*args, **kwargs)
        return json.loads(''.join(r))
    return call


def test_checksyntax_valid(api):
    """
    Running checkSyntax with a valid variant name should return True.
    """
    r = api('checkSyntax', variant='AB026906.1:c.274G>T')
    assert r == {'valid': True, 'messages': []}


def test_checksyntax_invalid(api):
    """
    Running checkSyntax with an invalid variant name should return False
    and give at least one error message.
    """
    r = api('checkSyntax', variant='0:abcd')
    assert not r['valid']
    assert len(r['messages']) >= 1


def test_checksyntax_empty(api):
    """
    Running checkSyntax with no variant name should raise exception.
    """
    # The validator doesn't work with NullServer, so we cannot do this
    # test. See https://github.com/arskom/spyne/issues/318
    # r = api('checkSyntax')
    # assert r['faultcode'] == 'Client.ValidationError'
    pass


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_transcriptinfo_valid(api):
    """
    Running transcriptInfo with valid arguments should get us a Transcript
    object.
    """
    r = api('transcriptInfo', LOVD_ver='123', build='hg19',
            accNo='NM_002001.2')
    assert r['trans_start'] == -99
    assert r['trans_stop'] == 1066
    assert r['CDS_stop'] == 774


def test_info(api):
    """
    Running the info method should give us some version information.
    """
    r = api('info')
    assert isinstance(r['versionParts'], list)
    assert r['version'] == mutalyzer.__version__


def test_info_announcement(api):
    """
    Running the info method should show us the current announcement
    """
    announce.set_announcement('Test announcement')
    r = api('info')
    assert isinstance(r['announcement'], unicode)
    assert r['announcement'] == 'Test announcement'

    announce.set_announcement('New announcement')
    r = api('info')
    assert isinstance(r['announcement'], unicode)
    assert r['announcement'] == 'New announcement'

    announce.unset_announcement()
    r = api('info')
    assert not r.get('announcement')


def test_checksyntax_unicode(api):
    """
    Run checkSyntax with an invalid variant description containing
    non-ASCII unicode characters.
    """
    r = api('checkSyntax', 'La Pe\xf1a')
    assert not r['valid']
    assert len(r['messages']) == 1
    assert r['messages'][0]['errorcode'] == 'EPARSE'
    assert r['messages'][0]['message'] == 'Expected W:(0123...) (at char 2), (line:1, col:3)'


@pytest.mark.usefixtures('db')
def test_batchjob_unicode(api):
    """
    Submit a batch job with non-ASCII unicode characters in the input
    file.
    """
    variants = ['\u2026AB026906.1:c.274G>T',
                '\u2026AL449423.14(CDKN2A_v002):c.5_400del']
    expected = [['\u2026AB026906.1:c.274G>T',
                 '(grammar): Expected W:(0123...) (at char 0), (line:1, col:1)'],
                ['\u2026AL449423.14(CDKN2A_v002):c.5_400del',
                 '(grammar): Expected W:(0123...) (at char 0), (line:1, col:1)']]

    data = '\n'.join(variants) + '\n'  # .encode('base64')

    result = api('submitBatchJob', data.encode('utf-8'), 'SyntaxChecker')
    job_id = unicode(result)

    result = api('monitorBatchJob', job_id)
    assert int(result) == len(variants)

    scheduler = Scheduler.Scheduler()
    scheduler.process()

    result = api('monitorBatchJob', job_id)
    assert int(result) == 0

    result = api('getBatchJob', job_id)
    result = result.decode('base64').decode('utf-8').strip().split('\n')[1:]
    assert expected == [line.split('\t') for line in result]


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gene_location(api):
    """
    Get outer coordinates for gene.
    """
    r = api('getGeneLocation', 'SDHD', 'hg19')

    assert r == {'gene': 'SDHD',
                 'start': 111957571,
                 'stop': 111966518,
                 'orientation': 'forward',
                 'chromosome_name': 'chr11',
                 'chromosome_accession': 'NC_000011.9',
                 'assembly_name': 'GRCh37',
                 'assembly_alias': 'hg19'}


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gene_location_reverse(api):
    """
    Get outer coordinates for gene on the reverse strand.
    """
    r = api('getGeneLocation', 'DMD', 'hg19')

    assert r == {'gene': 'DMD',
                 'start': 31137345,
                 'stop': 33038317,
                 'orientation': 'reverse',
                 'chromosome_name': 'chrX',
                 'chromosome_accession': 'NC_000023.10',
                 'assembly_name': 'GRCh37',
                 'assembly_alias': 'hg19'}


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gene_location_default_build(api):
    """
    Get outer coordinates for gene without specifying the build.
    """
    r = api('getGeneLocation', 'SDHD')

    assert r == {'gene': 'SDHD',
                 'start': 111957571,
                 'stop': 111966518,
                 'orientation': 'forward',
                 'chromosome_name': 'chr11',
                 'chromosome_accession': 'NC_000011.9',
                 'assembly_name': 'GRCh37',
                 'assembly_alias': 'hg19'}


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gene_location_invalid_gene(api):
    """
    Get outer coordinates for gene that does not exist.
    """
    with pytest.raises(Fault):
        api('getGeneLocation', 'THISISNOTAGENE', 'hg19')


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_get_transcripts_mapping(api):
    """
    Test output of getTranscriptsMapping.
    """
    r = api('getTranscriptsMapping', 'hg19', 'chr11', 111955524, 111966518)
    assert r == [{'cds_start': 111957492,
                  'cds_stop': 111956019,
                  'transcript': 'NM_012459.2',
                  'name': 'NM_012459',
                  'stop': 111955524,
                  'start': 111957522,
                  'version': 2,
                  'gene': 'TIMM8B',
                  'orientation': '-'},
                 {'cds_start': None,
                  'cds_stop': None,
                  'transcript': 'NR_028383.1',
                  'name': 'NR_028383',
                  'stop': 111955524,
                  'start': 111957522,
                  'version': 1,
                  'gene': 'TIMM8B',
                  'orientation': '-'},
                 {'cds_start': 111957632,
                  'cds_stop': 111965694,
                  'transcript': 'NM_003002.2',
                  'name': 'NM_003002',
                  'stop': 111966518,
                  'start': 111957571,
                  'version': 2,
                  'gene': 'SDHD',
                  'orientation': '+'}]


def test_description_extract(api):
    """
    Test output of descriptionExtract.
    """
    r = api('descriptionExtract',
            'ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA',
            'ATGATTTGATCAGATACATGTGATACCGGTAGTTAGGACAA')
    assert r == {'allele': [{'end': 6,
                             'deleted': '',
                             'weight': 8,
                             'inserted': 'TT',
                             'start_offset': 0,
                             'start': 5,
                             'description': '5_6insTT',
                             'shift': 1,
                             'end_offset': 0,
                             'type': 'ins',
                             'sample_start': 6,
                             'sample_end': 7,
                             'sample_start_offset': 0,
                             'sample_end_offset': 0},
                            {'end': 17,
                             'deleted': 'G',
                             'weight': 7,
                             'inserted': '',
                             'start_offset': 0,
                             'start': 17,
                             'description': '17del',
                             'shift': 0,
                             'end_offset': 0,
                             'type': 'del',
                             'sample_start': 18,
                             'sample_end': 19,
                             'sample_start_offset': 0,
                             'sample_end_offset': 0},
                            {'end': 26,
                             'deleted': 'A',
                             'weight': 3,
                             'inserted': 'C',
                             'start_offset': 0,
                             'start': 26,
                             'description': '26A>C',
                             'shift': 0,
                             'end_offset': 0,
                             'type': 'subst',
                             'sample_start': 27,
                             'sample_end': 27,
                             'sample_start_offset': 0,
                             'sample_end_offset': 0},
                            {'end': 35,
                             'deleted': '',
                             'weight': 5,
                             'inserted': 'G',
                             'start_offset': 0,
                             'start': 35,
                             'description': '35dup',
                             'shift': 1,
                             'end_offset': 0,
                             'type': 'dup',
                             'sample_start': 37,
                             'sample_end': 37,
                             'sample_start_offset': 0,
                             'sample_end_offset': 0}],
                 'description': '[5_6insTT;17del;26A>C;35dup]'}


def test_description_extract_ref_too_long(api):
    """
    Test output of descriptionExtract with too long reference sequence.
    """
    with pytest.raises(Fault):
        api('descriptionExtract',
            'A' * (settings.EXTRACTOR_MAX_INPUT_LENGTH + 1),
            'A')


def test_description_extract_sample_too_long(api):
    """
    Test output of descriptionExtract with too long sample sequence.
    """
    with pytest.raises(Fault):
        api('descriptionExtract',
            'A' * (settings.EXTRACTOR_MAX_INPUT_LENGTH),
            'A' * (settings.EXTRACTOR_MAX_INPUT_LENGTH + 1))


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_transcript_order(api):
    """
    Test whether getGeneLocation and numberConversion have same strategy for
    selecting a transcript mapping.
    """
    result_GL = api('getGeneLocation', gene='VAMP7', build='hg19')
    result_NC = api('numberConversion', build='hg19',
                    variant='NM_001145149.2:c.100dup', gene='VAMP7')

    assert 'chromosome_accession' in result_GL
    assert len(result_NC) == 1

    acc = result_GL['chromosome_accession']
    assert result_NC[0][:len(acc)] == acc
