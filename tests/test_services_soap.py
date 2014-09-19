"""
Tests for the SOAP interface to Mutalyzer.
"""


import bz2
import datetime
import logging
import os
import tempfile

from Bio import Entrez
from mock import patch
import pytest
from spyne.server.null import NullServer
from spyne.model.fault import Fault
from suds.client import Client

import mutalyzer
from mutalyzer.config import settings
from mutalyzer.output import Output
from mutalyzer.services.soap import application
from mutalyzer.sync import CacheSync
from mutalyzer import Scheduler

from fixtures import database, cache, hg19, hg19_transcript_mappings
from utils import MutalyzerTest
from utils import fix


# Suds logs an awful lot of things with level=DEBUG, including entire WSDL
# files and SOAP responses. On any error, this is all dumped to the console,
# which is very unconvenient. The following suppresses most of this.
logging.raiseExceptions = 0
logging.basicConfig(level=logging.INFO)
for logger in ('suds.metrics', 'suds.wsdl', 'suds.xsd.schema',
               'suds.xsd.sxbasic', 'suds.xsd.sxbase', 'suds.xsd.query',
               'suds.transport.http', 'suds.xsd.deplist', 'suds.mx.core',
               'suds.mx.literal', 'suds.resolver', 'suds.client',
               'suds.umx.typed'):
    logging.getLogger(logger).setLevel(logging.ERROR)


def _write_wsdl(server):
    server.doc.wsdl11.build_interface_document('/')
    wsdl = tempfile.NamedTemporaryFile(mode='w', delete=False)
    wsdl_filename = wsdl.name
    wsdl.write(server.doc.wsdl11.get_interface_document())
    wsdl.close()
    return wsdl_filename


class TestServicesSoap(MutalyzerTest):
    """
    Test the Mutalyzer SOAP interface.
    """
    def setup(self):
        super(TestServicesSoap, self).setup()
        self.server = NullServer(application, ostr=True)
        # Unfortunately there's no easy way to just give a SUDS client a
        # complete WSDL string, it only accepts a URL to it. So we create one.
        self.wsdl = _write_wsdl(self.server)
        self.client = Client('file://%s' % self.wsdl, cache=None)

    def teardown(self):
        super(TestServicesSoap, self).teardown()
        os.unlink(self.wsdl)

    def _call(self, method, *args, **kwargs):
        r = getattr(self.server.service, method)(*args, **kwargs)
        # This seems to be the way to feed raw SOAP response strings to a
        # SUDS client, without having it talking to a real server.
        return getattr(self.client.service, method)(__inject={'reply': ''.join(r)})

    def test_ping(self):
        """
        Running the ping method should return 'pong'.
        """
        r = self._call('ping')
        assert r == 'pong'

    def test_checksyntax_valid(self):
        """
        Running checkSyntax with a valid variant name should return True.
        """
        r = self._call('checkSyntax', 'AB026906.1:c.274G>T')
        assert r.valid == True

    def test_checksyntax_invalid(self):
        """
        Running checkSyntax with an invalid variant name should return False
        and give at least one error message.
        """
        r = self._call('checkSyntax', '0:abcd')
        assert r.valid == False
        assert len(r.messages.SoapMessage) >= 1

    def test_checksyntax_empty(self):
        """
        Running checkSyntax with no variant name should raise exception.
        """
        # The validator doesn't work with NullServer, so we cannot really do
        # these type of tests. However, in this case we implemented our own
        # check instead of relying on the validator.
        # See https://github.com/arskom/spyne/issues/318
        with pytest.raises(Fault):
            self._call('checkSyntax')

    @fix(database, hg19, hg19_transcript_mappings)
    def test_transcriptinfo_valid(self):
        """
        Running transcriptInfo with valid arguments should get us a Transcript
        object.
        """
        r = self._call('transcriptInfo',
                       LOVD_ver='123', build='hg19', accNo='NM_002001.2')
        assert r.trans_start == -99
        assert r.trans_stop == 1066
        assert r.CDS_stop == 774

    @fix(database, hg19, hg19_transcript_mappings)
    def test_numberconversion_gtoc_valid(self):
        """
        Running numberConversion with valid g variant should give a list of
        c variant names.
        """
        r = self._call('numberConversion',
                       build='hg19', variant='NC_000001.10:g.159272155del')
        assert type(r.string) == list
        assert 'NM_002001.2:c.1del' in r.string

    @fix(database, hg19, hg19_transcript_mappings)
    def test_numberconversion_ctog_valid(self):
        """
        Running numberConversion with valid c variant should give a list of
        g variant names.
        """
        r = self._call('numberConversion',
                       build='hg19', variant='NM_002001.2:c.1del')
        assert type(r.string) == list
        assert 'NC_000001.10:g.159272155del' in r.string

    @fix(database, hg19, hg19_transcript_mappings)
    def test_numberconversion_gtoc_gene(self):
        """
        Running numberConversion with valid g variant and a gene name should
        give a list of c variant names on transcripts for the given gene.
        """
        r = self._call('numberConversion',
                       build='hg19', variant='NC_000023.10:g.32827640G>A', gene='DMD')
        assert type(r.string) == list
        assert 'NM_004007.2:c.250C>T' in r.string
        assert 'NM_004011.3:c.-397314C>T' in r.string
        assert 'NM_004019.2:c.-1542694C>T' in r.string

    @fix(database, hg19, hg19_transcript_mappings)
    def test_numberconversion_gtoc_no_transcripts(self):
        """
        Running numberConversion with valid g variant but no transcripts
        close to it should give an empty list.
        """
        r = self._call('numberConversion',
                       build='hg19', variant='chr7:g.345T>C')
        assert not r

    @fix(database, hg19, hg19_transcript_mappings)
    def test_numberconversion_gtoc_required_gene(self):
        """
        Running numberConversion with valid g variant but no transcripts
        close to it, but with a gene name, should give a list of c variant
        names on transcripts for the given gene.
        """
        r = self._call('numberConversion',
                       build='hg19', variant='chr7:g.345T>C', gene='LOC100132858')
        assert type(r.string) == list
        # Fix for r536: disable the -u and +d convention.
        #assert 'XM_001715131.2:c.1155+d19483A>G' in r.string
        assert 'XM_001715131.2:c.*19483A>G' in r.string

    @fix(database, hg19, hg19_transcript_mappings)
    def test_gettranscriptsbygenename_valid(self):
        """
        Running getTranscriptsByGeneName with valid gene name should give a
        list of transcripts.
        """
        r = self._call('getTranscriptsByGeneName',
                       build='hg19', name='DMD')
        assert type(r.string) == list
        for t in ['NM_004011.3',
                  'NM_004019.2',
                  'NM_004007.2']:
            assert t in r.string

    @fix(database, hg19, hg19_transcript_mappings)
    def test_gettranscriptsbygenename_invalid(self):
        """
        Running getTranscriptsByGeneName with invalid gene name should not
        give a result.
        """
        r = self._call('getTranscriptsByGeneName',
                       build='hg19', name='BOGUSGENE')
        assert not r

    @fix(database, cache('AF230870.1'))
    def test_gettranscriptsandinfo_valid(self):
        """
        Running getTranscriptsAndInfo with a valid genomic reference should
        give a list of TranscriptInfo objects.
        """
        r = self._call('getTranscriptsAndInfo', 'AF230870.1')
        assert type(r.TranscriptInfo) == list
        names = [t.name for t in r.TranscriptInfo]
        for t in ['mtmC2_v001',
                  'mtmB2_v001']:
            assert t in names

    @fix(database, cache('AL449423.14'))
    def test_gettranscriptsandinfo_restricted_valid(self):
        """
        Running getTranscriptsAndInfo with a valid genomic reference and a
        gene name should give a list of TranscriptInfo objects restricted
        to the gene.
        """
        r = self._call('getTranscriptsAndInfo', 'AL449423.14', 'CDKN2A')
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

    @fix(database, hg19, hg19_transcript_mappings)
    def test_gettranscriptsmapping(self):
        """
        Running getTranscriptsMapping should give a list of
        TranscriptMappingInfo objects.
        """
        r = self._call('getTranscriptsMapping',
                       'hg19', 'chrX', 31200000, 31210000, 1)
        assert type(r.TranscriptMappingInfo) == list
        names = [t.name for t in r.TranscriptMappingInfo]
        for t in ('NM_004011',
                  'NM_004019',
                  'NM_004007'):
            assert t in names

    @fix(database, hg19, hg19_transcript_mappings)
    def test_mappinginfo(self):
        """
        Running mappingInfo should give a Mapping object.
        """
        r = self._call('mappingInfo',
                       '3.0-beta-06', 'hg19', 'NM_001100.3', 'g.112037014G>T')
        assert r.endoffset == 117529978
        assert r.start_g == 112037014
        assert r.startoffset == 117529978
        assert r.mutationType == "subst"
        assert r.end_g == 112037014
        assert r.startmain == 1388
        assert r.endmain == 1388

    @fix(database, hg19, hg19_transcript_mappings)
    def test_mappinginfo(self):
        """
        Running mappingInfo should give a Mapping object.
        """
        r = self._call('mappingInfo',
                       '3.0-beta-06', 'hg19', 'NM_002001.2', 'g.159272168G>T')
        assert r.endoffset == 0
        assert r.start_g == 159272168
        assert r.startoffset == 0
        assert r.mutationType == 'subst'
        assert r.end_g == 159272168
        assert r.startmain == 14
        assert r.endmain == 14

    @fix(database, hg19, hg19_transcript_mappings)
    def test_mappinginfo_compound(self):
        """
        Running mappingInfo with compound variant should give a Mapping
        object.
        """
        r = self._call('mappingInfo',
                       '3.0-beta-06', 'hg19', 'NM_002001.2', 'g.[159272168G>T;159272174T>A]')
        assert r.endoffset == 0
        assert r.start_g == 159272168
        assert r.startoffset == 0
        assert r.mutationType == 'compound'
        assert r.end_g == 159272174
        assert r.startmain == 14
        assert r.endmain == 20

    @fix(database, hg19, hg19_transcript_mappings)
    def test_mappinginfo_reverse(self):
        """
        Running mappingInfo on a reverse transcript should give a Mapping
        object.
        """
        r = self._call('mappingInfo',
                       '3.0-beta-06', 'hg19', 'NM_004011.3', 'g.31152229_31152239del')
        assert r.endoffset == 0
        assert r.start_g == 31152229
        assert r.startoffset == 0
        assert r.mutationType == 'del'
        assert r.end_g == 31152239
        assert r.startmain == 6981
        assert r.endmain == 6971

    @fix(database, hg19, hg19_transcript_mappings)
    def test_mappinginfo_compound_reverse(self):
        """
        Running mappingInfo with compound variant on a reverse transcript
        should give a Mapping object.
        """
        r = self._call('mappingInfo',
                       '3.0-beta-06', 'hg19', 'NM_004011.3', 'g.[31152229_31152232del;31152235_31152239del]')
        assert r.endoffset == 0
        assert r.start_g == 31152229
        assert r.startoffset == 0
        assert r.mutationType == 'compound'
        assert r.end_g == 31152239
        assert r.startmain == 6981
        assert r.endmain == 6971

    def test_info(self):
        """
        Running the info method should give us some version information.
        """
        r = self._call('info')
        assert type(r.versionParts.string) == list
        assert r.version == mutalyzer.__version__

    @fix(database, cache('AB026906.1', 'AL449423.14', 'NM_003002.2'))
    def test_getcache(self):
        """
        Running the getCache method should give us the expected number of
        cache entries.
        """
        created_since = datetime.datetime.today() - datetime.timedelta(days=14)

        output = Output(__file__)
        sync = CacheSync(output)

        r = self._call('getCache', created_since)
        assert len(r.CacheEntry) == 3

    def test_getdbsnpdescriptions(self):
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
            r = self._call('getdbSNPDescriptions', 'rs9919552')

        assert 'NC_000011.9:g.111959625C>T' in r.string
        assert 'NG_012337.2:g.7055C>T' in r.string
        assert 'NM_003002.3:c.204C>T' in r.string
        assert 'NP_002993.1:p.Ser68=' in r.string

    @fix(database, hg19, hg19_transcript_mappings)
    def test_gettranscripts(self):
        """
        Running getTranscripts should give a list of transcripts.
        """
        r = self._call('getTranscripts',
                       build='hg19', chrom='chrX', pos=32237295)
        assert type(r.string) == list
        for t in ['NM_004011',
                  'NM_004007']:
            assert t in r.string

    @fix(database, hg19, hg19_transcript_mappings)
    def test_gettranscripts_with_versions(self):
        """
        Running getTranscripts with versions=True should give a list
        of transcripts with version numbers.
        """
        r = self._call('getTranscripts',
                       build='hg19', chrom='chrX', pos=32237295, versions=True)
        assert type(r.string) == list
        for t in ['NM_004011.3',
                  'NM_004007.2']:
            assert t in r.string

    @fix(database, cache('NM_003002.2'))
    def test_runmutalyzer(self):
        """
        Just a runMutalyzer test.
        """
        r = self._call('runMutalyzer', 'NM_003002.2:c.274G>T')
        assert r.errors == 0
        assert r.genomicDescription == 'NM_003002.2:n.335G>T'
        assert 'NM_003002.2(SDHD_v001):c.274G>T' in r.transcriptDescriptions.string

    @fix(database)
    def test_runmutalyzer_reference_info_nm(self):
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
            r = self._call('runMutalyzer', 'NM_003002:c.274G>T')

        assert r.errors == 0
        assert r.referenceId == 'NM_003002.2'
        assert r.sourceId == 'NM_003002.2'
        assert r.sourceAccession == 'NM_003002'
        assert r.sourceVersion == '2'
        assert r.sourceGi == '222352156'
        assert r.molecule == 'n'

    @fix(database, cache('NM_003002.2'))
    def test_runmutalyzer_reference_info_nm_version(self):
        """
        Get reference info for an NM variant with version.
        """
        r = self._call('runMutalyzer', 'NM_003002.2:c.274G>T')
        assert r.errors == 0
        assert r.referenceId == 'NM_003002.2'
        assert r.sourceId == 'NM_003002.2'
        assert r.sourceAccession == 'NM_003002'
        assert r.sourceVersion == '2'
        assert r.sourceGi == '222352156'
        assert r.molecule == 'n'

    @fix(database, cache('LRG_1'))
    def test_runmutalyzer_reference_info_lrg(self):
        """
        Get reference info for an LRG variant.
        """
        r = self._call('runMutalyzer', 'LRG_1t1:c.266G>T')
        assert r.errors == 0
        assert r.referenceId == 'LRG_1'
        assert r.sourceId == 'LRG_1'
        assert r.molecule == 'g'

    @fix(database, cache('NG_012772.1'))
    def test_runmutalyzer_reference_info_ng(self):
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
            r = self._call('runMutalyzer', 'NG_012772:g.18964del')

        assert r.errors == 0
        assert r.referenceId == 'NG_012772.1'
        assert r.sourceId == 'NG_012772.1'
        assert r.sourceAccession == 'NG_012772'
        assert r.sourceVersion == '1'
        assert r.sourceGi == '256574794'
        assert r.molecule == 'g'

    @fix(database, cache('NG_009105.1'))
    def test_runmutalyzer_reference_info_ng_version(self):
        """
        Get reference info for an NG variant with version.
        """
        r = self._call('runMutalyzer', 'NG_009105.1:g.18964del')
        assert r.errors == 0
        assert r.referenceId == 'NG_009105.1'
        assert r.sourceId == 'NG_009105.1'
        assert r.sourceAccession == 'NG_009105'
        assert r.sourceVersion == '1'
        assert r.sourceGi == '216548283'
        assert r.molecule == 'g'

    @fix(database, cache('NG_012772.1'))
    def test_runmutalyzer_reference_info_gi(self):
        """
        Get reference info for a GI variant.
        """
        r = self._call('runMutalyzer', 'gi256574794:g.18964del')
        assert r.errors == 0
        assert r.referenceId == 'NG_012772.1'
        assert r.sourceId == 'NG_012772.1'
        assert r.sourceAccession == 'NG_012772'
        assert r.sourceVersion == '1'
        assert r.sourceGi == '256574794'
        assert r.molecule == 'g'

    @fix(database, cache('NM_000143.3'))
    def test_runmutalyzer_exons(self):
        """
        Exon table in runMutalyzer output.
        """
        r = self._call('runMutalyzer', 'NM_000143.3:c.630_636del')
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

    @fix(database, cache('AB026906.1', 'NM_003002.2', 'AL449423.14'))
    def test_batchjob(self):
        """
        Submit a batch job.
        """
        variants = ['AB026906.1(SDHD):g.7872G>T',
                    'NM_003002.2:c.3_4insG',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        data = '\n'.join(variants) + '\n' #.encode('base64')

        result = self._call('submitBatchJob', data, 'NameChecker')
        job_id = str(result)

        result = self._call('monitorBatchJob', job_id)
        assert int(result) == len(variants)

        scheduler = Scheduler.Scheduler()
        scheduler.process()

        result = self._call('monitorBatchJob', job_id)
        assert int(result) == 0

        result = self._call('getBatchJob', job_id)
        assert len(result.decode('base64').strip().split('\n')) - 1 == len(variants)

    @fix(database)
    def test_batchjob_newlines_unix(self):
        """
        Submit a batch job with UNIX newlines.
        """
        variants = ['AB026906.1(SDHD):g.7872G>T',
                    'NM_003002.2:c.3_4insG',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        data = '\n'.join(variants) + '\n'

        result = self._call('submitBatchJob', data, 'SyntaxChecker')
        job_id = str(result)

        result = self._call('monitorBatchJob', job_id)
        assert int(result) == len(variants)

        scheduler = Scheduler.Scheduler()
        scheduler.process()

        result = self._call('monitorBatchJob', job_id)
        assert int(result) == 0

    @fix(database)
    def test_batchjob_newlines_mac(self):
        """
        Submit a batch job with Mac newlines.
        """
        variants = ['AB026906.1(SDHD):g.7872G>T',
                    'NM_003002.2:c.3_4insG',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        data = '\r'.join(variants) + '\r'

        result = self._call('submitBatchJob', data, 'SyntaxChecker')
        job_id = str(result)

        result = self._call('monitorBatchJob', job_id)
        assert int(result) == len(variants)

        scheduler = Scheduler.Scheduler()
        scheduler.process()

        result = self._call('monitorBatchJob', job_id)
        assert int(result) == 0

    @fix(database)
    def test_batchjob_newlines_windows(self):
        """
        Submit a batch job with Windows newlines.
        """
        variants = ['AB026906.1(SDHD):g.7872G>T',
                    'NM_003002.2:c.3_4insG',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        data = '\r\n'.join(variants) + '\r\n'

        result = self._call('submitBatchJob', data, 'SyntaxChecker')
        job_id = str(result)

        result = self._call('monitorBatchJob', job_id)
        assert int(result) == len(variants)

        scheduler = Scheduler.Scheduler()
        scheduler.process()

        result = self._call('monitorBatchJob', job_id)
        assert int(result) == 0

    @fix(database)
    def test_batchjob_toobig(self):
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
            self._call('submitBatchJob', data.encode('base64'), 'NameChecker')
            assert False
        except Fault as e:
            # - senv:Client.RequestTooLong: Raised by Spyne, depending on
            #     the max_content_length argument to the HttpBase constructor.
            # - EMAXSIZE: Raised by Mutalyzer, depending on the
            #     batchInputMaxSize configuration setting.
            assert e.faultcode in ('senv:Client.RequestTooLong', 'EMAXSIZE')

    @fix(database)
    def test_upload_local_genbank(self):
        """
        Upload local genbank file.
        """
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data',
                            'AB026906.1.gb.bz2')
        with bz2.BZ2File(path) as f:
            data = f.read()

        result = self._call('upLoadGenBankLocalFile', data)
        ud = str(result)

        r = self._call('runMutalyzer', ud + '(SDHD):g.7872G>T')
        assert r.errors == 0
        assert r.genomicDescription == ud + ':g.7872G>T'
        assert ud + '(SDHD_v001):c.274G>T' in r.transcriptDescriptions.string
