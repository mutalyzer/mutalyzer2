"""
Tests for the mutalyzer.services.soap module.
"""


from __future__ import unicode_literals

import bz2
import datetime
import os

from Bio import Entrez
from mock import patch
import pytest
from spyne.server.null import NullServer
from spyne.model.fault import Fault
from suds.client import Client

import mutalyzer
from mutalyzer.services.soap import application
from mutalyzer import Scheduler

from fixtures import with_references


@pytest.fixture
def server():
    return NullServer(application, ostr=True)


@pytest.fixture
def wsdl(tmpdir, server):
    wsdl_file = tmpdir.join('wsdl').ensure()
    server.doc.wsdl11.build_interface_document('/')
    wsdl_file.write(server.doc.wsdl11.get_interface_document())
    return unicode(wsdl_file)


@pytest.fixture
def client(wsdl):
    return Client('file://%s' % wsdl, cache=None)


@pytest.fixture
def api(server, client):
    def call(method, *args, **kwargs):
        r = getattr(server.service, method)(*args, **kwargs)
        # This seems to be the way to feed raw SOAP response strings to a
        # SUDS client, without having it talking to a real server.
        return getattr(client.service, method)(__inject={'reply': ''.join(r)})
    return call


def test_ping(api):
    """
    Running the ping method should return 'pong'.
    """
    r = api('ping')
    assert r == 'pong'


def test_checksyntax_valid(api):
    """
    Running checkSyntax with a valid variant name should return True.
    """
    r = api('checkSyntax', 'AB026906.1:c.274G>T')
    assert r.valid


def test_checksyntax_invalid(api):
    """
    Running checkSyntax with an invalid variant name should return False
    and give at least one error message.
    """
    r = api('checkSyntax', '0:abcd')
    assert not r.valid
    assert len(r.messages.SoapMessage) >= 1


def test_checksyntax_empty(api):
    """
    Running checkSyntax with no variant name should raise exception.
    """
    # The validator doesn't work with NullServer, so we cannot really do
    # these type of tests. However, in this case we implemented our own
    # check instead of relying on the validator.
    # See https://github.com/arskom/spyne/issues/318
    with pytest.raises(Fault):
        api('checkSyntax')


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_transcriptinfo_valid(api):
    """
    Running transcriptInfo with valid arguments should get us a Transcript
    object.
    """
    r = api('transcriptInfo',
            LOVD_ver='123', build='hg19', accNo='NM_002001.2')
    assert r.trans_start == -99
    assert r.trans_stop == 1066
    assert r.CDS_stop == 774


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_numberconversion_gtoc_valid(api):
    """
    Running numberConversion with valid g variant should give a list of
    c variant names.
    """
    r = api('numberConversion',
            build='hg19', variant='NC_000001.10:g.159272155del')
    assert type(r.string) == list
    assert 'NM_002001.2:c.1del' in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_numberconversion_ctog_valid(api):
    """
    Running numberConversion with valid c variant should give a list of
    g variant names.
    """
    r = api('numberConversion',
            build='hg19', variant='NM_002001.2:c.1del')
    assert type(r.string) == list
    assert 'NC_000001.10:g.159272155del' in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_numberconversion_gtoc_gene(api):
    """
    Running numberConversion with valid g variant and a gene name should
    give a list of c variant names on transcripts for the given gene.
    """
    r = api('numberConversion',
            build='hg19', variant='NC_000023.10:g.32827640G>A', gene='DMD')
    assert type(r.string) == list
    assert 'NM_004007.2:c.250C>T' in r.string
    assert 'NM_004011.3:c.-397314C>T' in r.string
    assert 'NM_004019.2:c.-1542694C>T' in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_numberconversion_gtoc_no_transcripts(api):
    """
    Running numberConversion with valid g variant but no transcripts
    close to it should give an empty list.
    """
    r = api('numberConversion',
            build='hg19', variant='chr7:g.345T>C')
    assert not r


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_numberconversion_gtoc_required_gene(api):
    """
    Running numberConversion with valid g variant but no transcripts
    close to it, but with a gene name, should give a list of c variant
    names on transcripts for the given gene.
    """
    r = api('numberConversion',
            build='hg19', variant='chr7:g.345T>C', gene='LOC100132858')
    assert type(r.string) == list
    # Fix for r536: disable the -u and +d convention.
    # assert 'XM_001715131.2:c.1155+d19483A>G' in r.string
    assert 'XM_001715131.2:c.*19483A>G' in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gettranscripts_lrg(api):
    """
    Running getTranscripts should give us overlapping transcripts.
    list of transcripts including LRG transcripts.
    """
    r = api('getTranscripts', build='hg19', chrom='chr1',
            pos=207646118)
    assert type(r.string) == list
    assert 'LRG_348t1' in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gettranscripts_mtdna(api):
    """
    Running getTranscripts should give us overlapping transcripts.
    list of transcripts, also on chrM.
    """
    r = api('getTranscripts', build='hg19', chrom='chrM',
            pos=10765)
    assert type(r.string) == list
    assert 'NC_012920(ND4_v001)' in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gettranscriptsrange_lrg(api):
    """
    Running getTranscriptsRange should give us overlapping transcripts.
    list of transcripts including LRG transcripts.
    """
    r = api('getTranscriptsRange', 'hg19', 'chr1', 207646118, 207646118, 1)
    assert type(r.string) == list
    assert 'LRG_348t1' in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gettranscriptsrange_mtdna(api):
    """
    Running getTranscripts should give us overlapping transcripts.
    list of transcripts, also on chrM.
    """
    r = api('getTranscriptsRange', 'hg19', 'chrM', 10765, 10765, 1)
    assert type(r.string) == list
    assert 'NC_012920(ND4_v001)' in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gettranscriptsbygenename_valid(api):
    """
    Running getTranscriptsByGeneName with valid gene name should give a
    list of transcripts.
    """
    r = api('getTranscriptsByGeneName',  build='hg19', name='DMD')
    assert type(r.string) == list
    for t in ['NM_004011.3',
              'NM_004019.2',
              'NM_004007.2']:
        assert t in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gettranscriptsbygenename_valid_lrg(api):
    """
    Running getTranscriptsByGeneName with valid gene name should give a
    list of transcripts including LRG transcripts.
    """
    r = api('getTranscriptsByGeneName', build='hg19', name='CR2')
    assert type(r.string) == list
    assert 'LRG_348t1' in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gettranscriptsbygenename_valid_mtdna(api):
    """
    Running getTranscriptsByGeneName with valid gene name should give a
    list of transcripts also on chrM.
    """
    r = api('getTranscriptsByGeneName', build='hg19', name='ND4')
    assert type(r.string) == list
    assert 'NC_012920.1(ND4_v001)' in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gettranscriptsbygenename_invalid(api):
    """
    Running getTranscriptsByGeneName with invalid gene name should not
    give a result.
    """
    r = api('getTranscriptsByGeneName', build='hg19', name='BOGUSGENE')
    assert not r


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gettranscriptsmapping_mtdna(api):
    """
    Running getTranscriptsMapping on mtDNA should give a list of transcripts
    including their complete identification in the transcript attribute.
    """
    r = api('getTranscriptsMapping', 'hg19', 'chrM', 10765, 10765, 1)
    assert type(r.TranscriptMappingInfo) == list
    assert len(r.TranscriptMappingInfo) == 1
    info = r.TranscriptMappingInfo[0]
    assert 'NC_012920' == info.name
    assert 1 == info.version
    assert 'ND4' == info.gene
    assert 'NC_012920.1(ND4_v001)' == info.transcript


@with_references('AF230870.1')
def test_gettranscriptsandinfo_valid(api):
    """
    Running getTranscriptsAndInfo with a valid genomic reference should
    give a list of TranscriptInfo objects.
    """
    r = api('getTranscriptsAndInfo', 'AF230870.1')
    assert type(r.TranscriptInfo) == list
    names = [t.name for t in r.TranscriptInfo]
    for t in ['mtmC2_v001',
              'mtmB2_v001']:
        assert t in names


@with_references('AL449423.14')
def test_gettranscriptsandinfo_restricted_valid(api):
    """
    Running getTranscriptsAndInfo with a valid genomic reference and a
    gene name should give a list of TranscriptInfo objects restricted
    to the gene.
    """
    r = api('getTranscriptsAndInfo', 'AL449423.14', 'CDKN2A')
    assert type(r.TranscriptInfo) == list
    names = [t.name for t in r.TranscriptInfo]
    for t in ['CDKN2A_v008',
              'CDKN2A_v007']:
        assert t in names
    for t in ['CDKN2B_v002',
              'CDKN2B_v001',
              'MTAP_v005',
              'C9orf53_v001']:
        assert t not in names


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gettranscriptsmapping(api):
    """
    Running getTranscriptsMapping should give a list of
    TranscriptMappingInfo objects.
    """
    r = api('getTranscriptsMapping',
            'hg19', 'chrX', 31200000, 31210000, 1)
    assert type(r.TranscriptMappingInfo) == list
    names = [t.name for t in r.TranscriptMappingInfo]
    for t in ('NM_004011',
              'NM_004019',
              'NM_004007'):
        assert t in names


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_mappinginfo(api):
    """
    Running mappingInfo should give a Mapping object.
    """
    r = api('mappingInfo',
            '3.0-beta-06', 'hg19', 'NM_001100.3', 'g.112037014G>T')
    assert r.endoffset == 117529978
    assert r.start_g == 112037014
    assert r.startoffset == 117529978
    assert r.mutationType == "subst"
    assert r.end_g == 112037014
    assert r.startmain == 1388
    assert r.endmain == 1388


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_mappinginfo(api):
    """
    Running mappingInfo should give a Mapping object.
    """
    r = api('mappingInfo',
            '3.0-beta-06', 'hg19', 'NM_002001.2', 'g.159272168G>T')
    assert r.endoffset == 0
    assert r.start_g == 159272168
    assert r.startoffset == 0
    assert r.mutationType == 'subst'
    assert r.end_g == 159272168
    assert r.startmain == 14
    assert r.endmain == 14


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_mappinginfo_compound(api):
    """
    Running mappingInfo with compound variant should give a Mapping
    object.
    """
    r = api('mappingInfo',
            '3.0-beta-06', 'hg19', 'NM_002001.2', 'g.[159272168G>T;159272174T>A]')
    assert r.endoffset == 0
    assert r.start_g == 159272168
    assert r.startoffset == 0
    assert r.mutationType == 'compound'
    assert r.end_g == 159272174
    assert r.startmain == 14
    assert r.endmain == 20


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_mappinginfo_reverse(api):
    """
    Running mappingInfo on a reverse transcript should give a Mapping
    object.
    """
    r = api('mappingInfo',
            '3.0-beta-06', 'hg19', 'NM_004011.3', 'g.31152229_31152239del')
    assert r.endoffset == 0
    assert r.start_g == 31152229
    assert r.startoffset == 0
    assert r.mutationType == 'del'
    assert r.end_g == 31152239
    assert r.startmain == 6981
    assert r.endmain == 6971


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_mappinginfo_compound_reverse(api):
    """
    Running mappingInfo with compound variant on a reverse transcript
    should give a Mapping object.
    """
    r = api('mappingInfo',
            '3.0-beta-06', 'hg19', 'NM_004011.3', 'g.[31152229_31152232del;31152235_31152239del]')
    assert r.endoffset == 0
    assert r.start_g == 31152229
    assert r.startoffset == 0
    assert r.mutationType == 'compound'
    assert r.end_g == 31152239
    assert r.startmain == 6981
    assert r.endmain == 6971


def test_info(api):
    """
    Running the info method should give us some version information.
    """
    r = api('info')
    assert type(r.versionParts.string) == list
    assert r.version == mutalyzer.__version__


@with_references('AB026906.1', 'AL449423.14', 'NM_003002.2')
def test_getcache(output, api):
    """
    Running the getCache method should give us the expected number of
    cache entries.
    """
    created_since = datetime.datetime.today() - datetime.timedelta(days=14)
    r = api('getCache', created_since)
    assert len(r.CacheEntry) == 3


def test_getdbsnpdescriptions(api):
    """
    Running getdbSNPDescriptions method should give us the expected HGVS
    descriptions for the given dbSNP id.
    """
    # Patch Retriever.snpConvert to return rs9919552.
    def mock_efetch(*args, **kwargs):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data',
                            'rs9919552.xml.bz2')
        return bz2.BZ2File(path)

    with patch.object(Entrez, 'efetch', mock_efetch):
        r = api('getdbSNPDescriptions', 'rs9919552')

    assert 'NC_000011.9:g.111959625C>T' in r.string
    assert 'NG_012337.2:g.7055C>T' in r.string
    assert 'NM_003002.3:c.204C>T' in r.string
    assert 'NP_002993.1:p.Ser68=' in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gettranscripts(api):
    """
    Running getTranscripts should give a list of transcripts.
    """
    r = api('getTranscripts',
            build='hg19', chrom='chrX', pos=32237295)
    assert type(r.string) == list
    for t in ['NM_004011',
              'NM_004007']:
        assert t in r.string


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_gettranscripts_with_versions(api):
    """
    Running getTranscripts with versions=True should give a list
    of transcripts with version numbers.
    """
    r = api('getTranscripts',
            build='hg19', chrom='chrX', pos=32237295, versions=True)
    assert type(r.string) == list
    for t in ['NM_004011.3',
              'NM_004007.2']:
        assert t in r.string


@with_references('NM_003002.2')
def test_runmutalyzer(api):
    """
    Just a runMutalyzer test.
    """
    r = api('runMutalyzer', 'NM_003002.2:c.274G>T')
    assert r.errors == 0
    assert r.genomicDescription == 'NM_003002.2:n.335G>T'
    assert 'NM_003002.2(SDHD_v001):c.274G>T' in r.transcriptDescriptions.string


@pytest.mark.usefixtures('db')
def test_runmutalyzer_reference_info_nm(api):
    """
    Get reference info for an NM variant without version.
    """
    # Patch GenBankRetriever.fetch to return the contents of NM_003002.2
    # for NM_003002.
    def mock_efetch(*args, **kwargs):
        if kwargs.get('id') != 'NM_003002':
            return Entrez.efetch(*args, **kwargs)
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data',
                            'NM_003002.2.gb.bz2')
        return bz2.BZ2File(path)

    with patch.object(Entrez, 'efetch', mock_efetch):
        r = api('runMutalyzer', 'NM_003002:c.274G>T')

    assert r.errors == 0
    assert r.referenceId == 'NM_003002.2'
    assert r.sourceId == 'NM_003002.2'
    assert r.sourceAccession == 'NM_003002'
    assert r.sourceVersion == '2'
    assert r.molecule == 'n'


@with_references('NM_003002.2')
def test_runmutalyzer_reference_info_nm_version(api):
    """
    Get reference info for an NM variant with version.
    """
    r = api('runMutalyzer', 'NM_003002.2:c.274G>T')
    assert r.errors == 0
    assert r.referenceId == 'NM_003002.2'
    assert r.sourceId == 'NM_003002.2'
    assert r.sourceAccession == 'NM_003002'
    assert r.sourceVersion == '2'
    assert r.molecule == 'n'


@with_references('LRG_1')
def test_runmutalyzer_reference_info_lrg(api):
    """
    Get reference info for an LRG variant.
    """
    r = api('runMutalyzer', 'LRG_1t1:c.266G>T')
    assert r.errors == 0
    assert r.referenceId == 'LRG_1'
    assert r.sourceId == 'LRG_1'
    assert r.molecule == 'g'


@with_references('NG_012772.1')
def test_runmutalyzer_reference_info_ng(api):
    """
    Get reference info for an NG variant without version.
    """
    # Patch GenBankRetriever.fetch to return the contents of NG_012772.1
    # for NG_012772.
    def mock_efetch(*args, **kwargs):
        if kwargs.get('id') != 'NG_012772':
            return Entrez.efetch(*args, **kwargs)
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data',
                            'NG_012772.1.gb.bz2')
        return bz2.BZ2File(path)

    with patch.object(Entrez, 'efetch', mock_efetch):
        r = api('runMutalyzer', 'NG_012772:g.18964del')

    assert r.errors == 0
    assert r.referenceId == 'NG_012772.1'
    assert r.sourceId == 'NG_012772.1'
    assert r.sourceAccession == 'NG_012772'
    assert r.sourceVersion == '1'
    assert r.molecule == 'g'


@with_references('NG_009105.1')
def test_runmutalyzer_reference_info_ng_version(api):
    """
    Get reference info for an NG variant with version.
    """
    r = api('runMutalyzer', 'NG_009105.1:g.18964del')
    assert r.errors == 0
    assert r.referenceId == 'NG_009105.1'
    assert r.sourceId == 'NG_009105.1'
    assert r.sourceAccession == 'NG_009105'
    assert r.sourceVersion == '1'
    assert r.molecule == 'g'


@with_references('NG_012772.1')
def test_runmutalyzer_reference_info_gi(api):
    """
    Get reference info for a GI variant.
    """
    r = api('runMutalyzer', 'gi256574794:g.18964del')
    assert r.errors == 1
    assert len(r.messages.SoapMessage) == 1
    assert r.messages.SoapMessage[0]['errorcode'] == 'EGINOTSUPPORTED'


@with_references('NM_000143.3')
def test_runmutalyzer_exons(api):
    """
    Exon table in runMutalyzer output.
    """
    r = api('runMutalyzer', 'NM_000143.3:c.630_636del')
    assert r.errors == 0
    expected_exons = [(1, 195, '-63', '132'),
                      (196, 330, '133', '267'),
                      (331, 441, '268', '378'),
                      (442, 618, '379', '555'),
                      (619, 801, '556', '738'),
                      (802, 967, '739', '904'),
                      (968, 1171, '905', '1108'),
                      (1172, 1299, '1109', '1236'),
                      (1300, 1453, '1237', '1390'),
                      (1454, 1867, '1391', '*271')]
    assert len(r.exons.ExonInfo) == len(expected_exons)
    for exon, expected_exon in zip(r.exons.ExonInfo, expected_exons):
        assert (exon.gStart, exon.gStop, exon.cStart, exon.cStop) == expected_exon


@with_references('AB026906.1', 'NM_003002.2', 'AL449423.14')
def test_batchjob(api):
    """
    Submit a batch job.
    """
    variants = ['AB026906.1(SDHD):g.7872G>T',
                'NM_003002.2:c.3_4insG',
                'AL449423.14(CDKN2A_v002):c.5_400del']
    data = '\n'.join(variants) + '\n' #.encode('base64')

    result = api('submitBatchJob', data.encode('utf-8'), 'NameChecker')
    job_id = unicode(result)

    result = api('monitorBatchJob', job_id)
    assert int(result) == len(variants)

    scheduler = Scheduler.Scheduler()
    scheduler.process()

    result = api('monitorBatchJob', job_id)
    assert int(result) == 0

    result = api('getBatchJob', job_id)
    assert len(result.decode('base64').strip().split('\n')) - 1 == len(variants)


@pytest.mark.usefixtures('db')
def test_batchjob_newlines_unix(api):
    """
    Submit a batch job with UNIX newlines.
    """
    variants = ['AB026906.1(SDHD):g.7872G>T',
                'NM_003002.2:c.3_4insG',
                'AL449423.14(CDKN2A_v002):c.5_400del']
    data = '\n'.join(variants) + '\n'

    result = api('submitBatchJob', data.encode('utf-8'), 'SyntaxChecker')
    job_id = unicode(result)

    result = api('monitorBatchJob', job_id)
    assert int(result) == len(variants)

    scheduler = Scheduler.Scheduler()
    scheduler.process()

    result = api('monitorBatchJob', job_id)
    assert int(result) == 0


@pytest.mark.usefixtures('db')
def test_batchjob_newlines_mac(api):
    """
    Submit a batch job with Mac newlines.
    """
    variants = ['AB026906.1(SDHD):g.7872G>T',
                'NM_003002.2:c.3_4insG',
                'AL449423.14(CDKN2A_v002):c.5_400del']
    data = '\r'.join(variants) + '\r'

    result = api('submitBatchJob', data.encode('utf-8'), 'SyntaxChecker')
    job_id = unicode(result)

    result = api('monitorBatchJob', job_id)
    assert int(result) == len(variants)

    scheduler = Scheduler.Scheduler()
    scheduler.process()

    result = api('monitorBatchJob', job_id)
    assert int(result) == 0


@pytest.mark.usefixtures('db')
def test_batchjob_newlines_windows(api):
    """
    Submit a batch job with Windows newlines.
    """
    variants = ['AB026906.1(SDHD):g.7872G>T',
                'NM_003002.2:c.3_4insG',
                'AL449423.14(CDKN2A_v002):c.5_400del']
    data = '\r\n'.join(variants) + '\r\n'

    result = api('submitBatchJob', data.encode('utf-8'), 'SyntaxChecker')
    job_id = unicode(result)

    result = api('monitorBatchJob', job_id)
    assert int(result) == len(variants)

    scheduler = Scheduler.Scheduler()
    scheduler.process()

    result = api('monitorBatchJob', job_id)
    assert int(result) == 0


@pytest.mark.usefixtures('db')
def test_batchjob_toobig(settings, api):
    """
    Submit the batch name checker with a too big input file.
    """
    seed = """
Lorem ipsum dolor sit amet, consectetuer adipiscing elit, sed diam nonummy
nibh euismod tincidunt ut laoreet dolore magna aliquam erat volutpat. Ut wisi
enim ad minim veniam, quis nostrud exerci tation ullamcorper suscipit lobortis
nisl ut aliquip ex ea commodo consequat. Duis autem vel eum iriure dolor in
hendrerit in vulputate velit esse molestie consequat, vel illum dolore eu
feugiat nulla facilisis at vero eros et accumsan et iusto odio dignissim qui
blandit praesent luptatum zzril delenit augue duis dolore te feugait nulla
facilisi."""
    data = seed
    # Very crude way of creating something big.
    while len(data) <= settings.MAX_FILE_SIZE:
        data += data

    try:
        api('submitBatchJob', data.encode('utf-8'), 'NameChecker')
        assert False
    except Fault as e:
        # - senv:Client.RequestTooLong: Raised by Spyne, depending on
        #     the max_content_length argument to the HttpBase constructor.
        # - EMAXSIZE: Raised by Mutalyzer, depending on the
        #     batchInputMaxSize configuration setting.
        assert e.faultcode in ('senv:Client.RequestTooLong', 'EMAXSIZE')


@pytest.mark.usefixtures('db')
def test_upload_local_genbank(api):
    """
    Upload local genbank file.
    """
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'data',
                        'AB026906.1.gb.bz2')
    with bz2.BZ2File(path) as f:
        data = f.read()

    result = api('uploadGenBankLocalFile', data)
    ud = unicode(result)

    r = api('runMutalyzer', ud + '(SDHD):g.7872G>T')
    assert r.errors == 0
    assert r.genomicDescription == ud + ':g.7872G>T'
    assert ud + '(SDHD_v001):c.274G>T' in r.transcriptDescriptions.string


def test_checksyntax_unicode(api):
    """
    Run checkSyntax with an invalid variant description containing
    non-ASCII unicode characters.
    """
    r = api('checkSyntax', 'La Pe\xf1a')
    assert not r.valid
    assert len(r.messages.SoapMessage) == 1
    assert r.messages.SoapMessage[0]['errorcode'] == 'EPARSE'
    assert r.messages.SoapMessage[0]['message'] ==  'Expected W:(0123...) (at char 2), (line:1, col:3)'


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
def test_get_transcripts_mapping(api):
    """
    Test output of getTranscriptsMapping.
    """
    r = api('getTranscriptsMapping', 'hg19', 'chr11',
            111955524, 111966518)
    assert len(r.TranscriptMappingInfo) == 3
    assert all(all(t_real[k] == t_expected[k] for k in t_expected)
               for t_real, t_expected in
               zip(r.TranscriptMappingInfo, [{'cds_start': 111957492,
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
                                              'orientation': '+'}]))


def test_description_extract(api):
    """
    Test output of descriptionExtract.
    """
    r = api('descriptionExtract',
            'ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA',
            'ATGATTTGATCAGATACATGTGATACCGGTAGTTAGGACAA')
    assert r['description'] == '[5_6insTT;17del;26A>C;35dup]'
    assert len(r['allele'].RawVar) == 4
    # For some reason, we get the empty string as `None` in SOAP.
    assert all(all((v_real[k] == v_expected[k]) or not(v_real[k] or v_expected[k])
                   for k in v_expected)
               for v_real, v_expected in
               zip(r['allele'].RawVar, [{'end': 6,
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
                                         'sample_end_offset': 0}]))
