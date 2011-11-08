"""
Tests for the SOAP interface to Mutalyzer.
"""


from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import os
from datetime import datetime, timedelta
import mutalyzer
from mutalyzer.output import Output
from mutalyzer.sync import CacheSync
from mutalyzer import Db
import logging
import urllib2
from suds.client import Client
from suds import WebFault
from nose.tools import *


# Suds logs an awful lot of things with level=DEBUG, including entire WSDL
# files and SOAP responses. On any error, this is all dumped to the console,
# which is very unconvenient. The following suppresses most of this.
logging.raiseExceptions = 0
logging.basicConfig(level=logging.INFO)
for logger in ('suds.metrics', 'suds.wsdl', 'suds.xsd.schema',
               'suds.xsd.sxbasic', 'suds.xsd.sxbase', 'suds.xsd.query',
               'suds.transport.http', 'suds.xsd.deplist', 'suds.mx.core',
               'suds.mx.literal', 'suds.resolver'):
    logging.getLogger(logger).setLevel(logging.ERROR)


WSDL_URL = 'http://localhost/mutalyzer/services/?wsdl'


class TestWSDL():
    """
    Test the Mutalyzer SOAP interface WSDL description.
    """
    def test_wsdl(self):
        """
        Test if the WSDL is available and looks somewhat sensible.
        """
        wsdl = urllib2.urlopen(WSDL_URL).read()
        assert wsdl.startswith("<?xml version='1.0' encoding='UTF-8'?>")
        assert 'name="Mutalyzer"' in wsdl


class TestWebservice():
    """
    Test the Mutalyzer SOAP interface.
    """

    def setUp(self):
        """
        Initialize webservice entrypoint.

        @todo: Start the standalone server and stop it in self.tearDown
        instead of depending on some running instance at a fixed address.
        """
        self.client = Client(WSDL_URL) #, cache=None)
        self.client.options.cache.setduration(seconds=120)

    def test_checksyntax_valid(self):
        """
        Running checkSyntax with a valid variant name should return True.
        """
        r = self.client.service.checkSyntax('AB026906.1:c.274G>T')
        assert_equal(r.valid, True)

    def test_checksyntax_invalid(self):
        """
        Running checkSyntax with an invalid variant name should return False
        and give at least one error message.
        """
        r = self.client.service.checkSyntax('0:abcd')
        assert_equal(r.valid, False)
        assert len(r.messages.SoapMessage) >= 1

    @raises(WebFault)
    def test_checksyntax_empty(self):
        """
        Running checkSyntax with no variant name should raise exception.
        """
        self.client.service.checkSyntax()

    def test_transcriptinfo_valid(self):
        """
        Running transcriptInfo with valid arguments should get us a Transcript
        object.
        """
        r = self.client.service.transcriptInfo(LOVD_ver='123', build='hg19',
                                               accNo='NM_002001.2')
        assert_equal(r.trans_start, -99)
        assert_equal(r.trans_stop, 1066)
        assert_equal(r.CDS_stop, 774)

    def test_numberconversion_gtoc_valid(self):
        """
        Running numberConversion with valid g variant should give a list of
        c variant names.
        """
        r = self.client.service.numberConversion(build='hg19',
                                                 variant='NC_000001.10:g.159272155del')
        assert_equal(type(r.string), list)
        assert 'NM_002001.2:c.1del' in r.string

    def test_numberconversion_ctog_valid(self):
        """
        Running numberConversion with valid c variant should give a list of
        g variant names.
        """
        r = self.client.service.numberConversion(build='hg19',
                                                 variant='NM_002001.2:c.1del')
        assert_equal(type(r.string), list)
        assert 'NC_000001.10:g.159272155del' in r.string

    def test_numberconversion_gtoc_gene(self):
        """
        Running numberConversion with valid g variant and a gene name should
        give a list of c variant names on transcripts for the given gene.
        """
        r = self.client.service.numberConversion(build='hg19',
                                                 variant='NC_000011.9:g.111959693G>T',
                                                 gene='C11orf57')
        assert_equal(type(r.string), list)
        assert 'NM_001082969.1:c.*2178+d3819G>T' in r.string
        assert 'NM_001082970.1:c.*2178+d3819G>T' in r.string
        assert 'NM_018195.3:c.*2178+d3819G>T' in r.string

    def test_numberconversion_gtoc_no_transcripts(self):
        """
        Running numberConversion with valid g variant but no transcripts
        close to it should give an empty list.
        """
        r = self.client.service.numberConversion(build='hg19',
                                                 variant='chr7:g.345T>C')
        assert_false(r)

    def test_numberconversion_gtoc_required_gene(self):
        """
        Running numberConversion with valid g variant but no transcripts
        close to it, but with a gene name, should give a list of c variant
        names on transcripts for the given gene.
        """
        r = self.client.service.numberConversion(build='hg19',
                                                 variant='chr7:g.345T>C',
                                                 gene='LOC100132858')
        assert_equal(type(r.string), list)
        assert 'XM_001715131.2:c.1155+d19483A>G' in r.string

    def test_gettranscriptsbygenename_valid(self):
        """
        Running getTranscriptsByGeneName with valid gene name should give a
        list of transcripts.
        """
        r = self.client.service.getTranscriptsByGeneName(build='hg19',
                                                         name='DMD')
        assert_equal(type(r.string), list)
        for t in ['NM_004006.2',
                  'NM_000109.3',
                  'NM_004021.2',
                  'NM_004009.3',
                  'NM_004007.2',
                  'NM_004018.2',
                  'NM_004022.2']:
            assert t in r.string

    def test_gettranscriptsbygenename_invalid(self):
        """
        Running getTranscriptsByGeneName with invalid gene name should not
        give a result.
        """
        r = self.client.service.getTranscriptsByGeneName(build='hg19',
                                                         name='BOGUSGENE')
        assert_false(r)

    def test_gettranscriptsandinfo_valid(self):
        """
        Running getTranscriptsAndInfo with a valid genomic reference should
        give a list of TranscriptInfo objects.
        """
        r = self.client.service.getTranscriptsAndInfo('AL449423.14')
        assert_equal(type(r.TranscriptInfo), list)
        names = [t.name for t in r.TranscriptInfo]
        for t in ['CDKN2B_v002',
                  'CDKN2B_v001',
                  'MTAP_v005',
                  'CDKN2A_v008',
                  'CDKN2A_v007',
                  'C9orf53_v001',
                  'CDKN2A_v001']:
            assert t in names

    def test_gettranscriptsandinfo_restricted_valid(self):
        """
        Running getTranscriptsAndInfo with a valid genomic reference and a
        gene name should give a list of TranscriptInfo objects restricted
        to the gene.
        """
        r = self.client.service.getTranscriptsAndInfo('AL449423.14', 'CDKN2A')
        assert_equal(type(r.TranscriptInfo), list)
        names = [t.name for t in r.TranscriptInfo]
        for t in ['CDKN2A_v008',
                  'CDKN2A_v007']:
            assert t in names
        for t in ['CDKN2B_v002',
                  'CDKN2B_v001',
                  'MTAP_v005',
                  'C9orf53_v001']:
            assert_false(t in names)

    def test_info(self):
        """
        Running the info method should give us some version information.
        """
        r = self.client.service.info()
        assert_equal(type(r.versionParts.string), list)
        assert_equal(r.version, mutalyzer.__version__)

    def test_getcache(self):
        """
        Running the getCache method should give us the expected number of
        cache entries.
        """
        created_since = datetime.today() - timedelta(days=14)

        database = Db.Cache()
        output = Output(__file__)
        sync = CacheSync(output, database)
        cache = sync.local_cache(created_since)

        r = self.client.service.getCache(created_since)
        if len(cache) > 0:
            assert_equal(len(r.CacheEntry), len(cache))

    def test_getdbsnpdescriptions(self):
        """
        Running getdbSNPDescriptions method should give us the expected HGVS
        descriptions for the given dbSNP id.
        """
        r = self.client.service.getdbSNPDescriptions('rs9919552')
        assert 'NC_000011.9:g.111959625C>T' in r.string
        assert 'NG_012337.1:g.7055C>T' in r.string
        assert 'NM_003002.2:c.204C>T' in r.string
        assert 'NP_002993.1:p.Ser68=' in r.string

    def test_gettranscripts(self):
        """
        Running getTranscripts should give a list of transcripts.
        """
        r = self.client.service.getTranscripts(build='hg19', chrom='chrX',
                                               pos=32237295)
        assert_equal(type(r.string), list)
        for t in ['NM_000109',
                  'NM_004006',
                  'NM_004007',
                  'NM_004009',
                  'NM_004010',
                  'NM_004011',
                  'NM_004012']:
            assert t in r.string

    def test_gettranscripts_with_versions(self):
        """
        Running getTranscripts with versions=True should give a list
        of transcripts with version numbers.
        """
        r = self.client.service.getTranscripts(build='hg19', chrom='chrX',
                                               pos=32237295, versions=True)
        assert_equal(type(r.string), list)
        for t in ['NM_000109.3',
                  'NM_004006.2',
                  'NM_004007.2',
                  'NM_004009.3',
                  'NM_004010.3',
                  'NM_004011.3',
                  'NM_004012.3']:
            assert t in r.string

    def test_ping(self):
        """
        Running the ping method should return 'pong'.
        """
        r = self.client.service.ping()
        assert_equal(r, 'pong')
