"""
Tests for the WSGI interface to Mutalyzer.

@todo: Tests for /upload.
"""


from __future__ import unicode_literals

#import logging; logging.basicConfig()
import bz2
from mock import patch
import os
from io import BytesIO

from Bio import Entrez
import lxml.html

from mutalyzer import announce, Scheduler
from mutalyzer.db import models
from mutalyzer.website import create_app

from fixtures import cache, database, hg19, hg19_transcript_mappings
from utils import MutalyzerTest
from utils import fix


BATCH_RESULT_URL = 'http://localhost/mutalyzer/Results_{id}.txt'


class TestWebsite(MutalyzerTest):
    """
    Test the Mutalyzer WSGI interface.
    """
    def setup(self):
        super(TestWebsite, self).setup()
        self.app = create_app().test_client()

    def test_homepage(self):
        """
        Expect the index HTML page.
        """
        r = self.app.get('/')
        assert r.status_code == 200
        assert 'Welcome to the Mutalyzer website' in r.data

    def test_about(self):
        """
        See if people get proper credit.
        """
        r = self.app.get('/about')
        assert r.status == '200 OK'
        assert 'Jonathan Vis' in r.data

    def test_non_existing(self):
        """
        Expect a 404 response.
        """
        r = self.app.get('/this/doesnotexist')
        assert r.status_code == 404

    @fix(database)
    def test_menu_links(self):
        """
        Test all links in the main menu.
        """
        ignore = []  # This could contain relative links we want to skip
        r = self.app.get('/')

        dom = lxml.html.fromstring(r.data)

        for link in dom.cssselect('#menu a'):
            href = link.get('href')
            if (href.startswith('http://') or
                href.startswith('https://') or
                href in ignore):
                continue
            if not href.startswith('/'):
                href = '/' + href

            r = self.app.get(href)
            assert r.status_code == 200

    def test_announcement(self):
        """
        We should always see the current announcement.
        """
        announce.set_announcement('Test announcement')
        r = self.app.get('/syntax-checker')
        assert r.status == '200 OK'
        assert 'Test announcement' in r.data

        announce.set_announcement('New announcement')
        r = self.app.get('/syntax-checker')
        assert r.status == '200 OK'
        assert 'New announcement' in r.data

        announce.unset_announcement()
        r = self.app.get('/syntax-checker')
        assert r.status == '200 OK'
        assert 'nnouncement' not in r.data

    def test_description_extractor(self):
        """
        Submit the variant description extractor.
        """
        r = self.app.get('/description-extractor', query_string={
                'reference_sequence': 'ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA',
                'variant_sequence': 'ATGATTTGATCAGATACATGTGATACCGGTAGTTAGGACAA'})
        assert 'g.[5_6insTT;17del;26A&gt;C;35dup]' in r.data

    def test_checksyntax_valid(self):
        """
        Submit the check syntax form with a valid variant.
        """
        r = self.app.get('/syntax-checker',
                         query_string={'description': 'AB026906.1:c.274G>T'})
        assert 'The syntax of this variant is OK!' in r.data

    def test_checksyntax_invalid(self):
        """
        Submit the check syntax form with an invalid variant.
        """
        r = self.app.get('/syntax-checker',
                         query_string={'description': 'AB026906.1:c.27'})
        assert 'Fatal' in r.data
        assert 'The &quot;^&quot; indicates the position where the error occurred' in r.data

    @fix(database, cache('NM_002001.2'))
    def test_check_valid(self):
        """
        Submit the name checker form with a valid variant.
        Should include form and main HTML layout.
        """
        r = self.app.get('/name-checker',
                         query_string={'description': 'NM_002001.2:g.1del'})
        assert '0 Errors' in r.data
        assert '0 Warnings' in r.data
        assert 'Raw variant 1: deletion of 1' in r.data
        assert '<input class="form-control form-pre example-target" type="text" name="description" id="hgvs" value="NM_002001.2:g.1del" placeholder="Mutation name using HGVS format">' in r.data

    def test_check_invalid(self):
        """
        Submit the name checker form with an invalid variant.
        """
        r = self.app.get('/name-checker',
                         query_string={'description': 'NM_002001.2'})
        assert '1 Error' in r.data
        assert '0 Warnings' in r.data
        assert 'The &quot;^&quot; indicates the position where the error occurred' in r.data

    @fix(database, cache('NP_064445.1'))
    def test_check_protein_reference(self):
        """
        Submit the name checker form with a protein reference sequence (not
        supported).
        """
        r = self.app.get('/name-checker',
                         query_string={'description': 'NP_064445.1:c.274G>T'})
        assert '1 Error' in r.data
        assert '0 Warnings' in r.data
        assert 'Protein reference sequences are not supported' in r.data

    @fix(database, cache('NM_002001.2'))
    def test_check_noninteractive(self):
        """
        Submit the name checker form non-interactively.
        Should not include form and main layout HTML.
        """
        r = self.app.get('/name-checker',
                         query_string={'description': 'NM_002001.2:g.1del',
                                       'standalone': '1'})
        assert '<a href="#bottom" class="hornav">go to bottom</a>' not in r.data
        assert '<input type="text" name="description" value="NM_002001.2:g.1del" style="width:100%">' not in r.data
        assert '0 Errors' in r.data
        assert '0 Warnings' in r.data
        assert 'Raw variant 1: deletion of 1' in r.data

    @fix(database, cache('NG_012772.1'))
    def test_check_interactive_links(self):
        """
        Submitting interactively should have links to transcripts also
        interactive.
        """
        r = self.app.get('/name-checker',
                         query_string={'description': 'NG_012772.1:g.128del'})
        assert '0 Errors' in r.data
        assert 'href="/name-checker?description=NG_012772.1%3Ag.128del"' in r.data
        assert 'href="/name-checker?description=NG_012772.1%28BRCA2_v001%29%3Ac.-5100del"' in r.data

    def test_snp_converter_valid(self):
        """
        Submit the SNP converter form with a valid SNP.
        """
        # Patch Retriever.snpConvert to return rs9919552.
        def mock_efetch(*args, **kwargs):
            path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'data',
                                'rs9919552.xml.bz2')
            return bz2.BZ2File(path)

        with patch.object(Entrez, 'efetch', mock_efetch):
            r = self.app.get('/snp-converter',
                             query_string={'rs_id': 'rs9919552'})
        assert '0 Errors' in r.data
        assert '0 Warnings' in r.data
        assert 'NC_000011.9:g.111959625C&gt;T' in r.data
        assert 'NG_012337.2:g.7055C&gt;T' in r.data
        assert 'NM_003002.3:c.204C&gt;T' in r.data
        assert 'NP_002993.1:p.Ser68=' in r.data

    def test_snp_converter_invalid(self):
        """
        Submit the SNP converter form with an invalid SNP.
        """
        r = self.app.get('/snp-converter',
                         query_string={'rs_id': 'r9919552'})

        assert '1 Error' in r.data
        assert '0 Warnings' in r.data
        assert 'Fatal' in r.data
        assert 'This is not a valid dbSNP id' in r.data

    @fix(database, hg19, hg19_transcript_mappings)
    def test_position_converter_c2g(self):
        """
        Submit the position converter form with a valid variant.
        """
        r = self.app.get('/position-converter',
                         query_string={'assembly_name_or_alias': 'hg19',
                                       'description': 'NM_003002.2:c.204C>T'})
        assert 'NC_000011.9:g.111959625C&gt;T' in r.data

    @fix(database, hg19, hg19_transcript_mappings)
    def test_position_converter_g2c(self):
        """
        Submit the position converter form with a valid variant.
        """
        r = self.app.get('/position-converter',
                         query_string={'assembly_name_or_alias': 'hg19',
                                       'description': 'NC_000011.9:g.111959625C>T'})
        assert 'NM_003002.2:c.204C&gt;T' in r.data

    def _batch(self, job_type='name-checker', assembly_name_or_alias=None,
               file="", size=0, header='', lines=None):
        """
        Submit a batch form.

        @kwarg batch_type: Type of batch job to test. One of name-checker,
                           syntax-checker, position-converter.
        @kwarg argument: Optional extra argument for the batch job.
        @kwarg file: String with variants to use as input for the batch job.
        @kwarg size: Number of variants in input.
        @kwarg header: Message that must be found in the batch job result.
        @kwarg lines: Number of result rows expected.

        @return: The batch result document.
        @rtype: string
        """
        data = {'job_type': job_type,
                'email': 'test@test.test',
                'file': (BytesIO(file.encode('utf-8')), 'test.txt')}
        if assembly_name_or_alias is not None:
            data['assembly_name_or_alias'] = assembly_name_or_alias

        r = self.app.post('/batch-jobs',
                          data=data)
        progress_url = '/' + r.location.split('/')[-1]

        r = self.app.get(progress_url)
        assert '<div id="if_items_left">' in r.data
        assert '<div id="ifnot_items_left" style="display:none">' in r.data
        assert ('<span id="items_left">%d</span>' % size) in r.data

        scheduler = Scheduler.Scheduler()
        scheduler.process()

        r = self.app.get(progress_url)
        assert '<div id="if_items_left" style="display:none">' in r.data
        assert '<div id="ifnot_items_left">' in r.data

        dom = lxml.html.fromstring(r.data)
        result_url = dom.cssselect('#ifnot_items_left a')[0].attrib['href']

        if not lines:
            lines = size

        r = self.app.get(result_url)
        assert 'text/plain' in r.headers['Content-Type']
        assert header in r.data
        assert len(r.data.strip().split('\n')) - 1 == lines

        return r.data

    @fix(database, cache('AB026906.1', 'NM_003002.2', 'AL449423.14'))
    def test_batch_namechecker(self):
        """
        Submit the batch name checker form.
        """
        variants=['AB026906.1(SDHD):g.7872G>T',
                  'NM_003002.2:c.3_4insG',
                  'AL449423.14(CDKN2A_v002):c.5_400del']
        self._batch('name-checker',
                    file='\n'.join(variants),
                    size=len(variants),
                    header='Input\tErrors and warnings')

    @fix(database)
    def test_batch_namechecker_extra_tab(self):
        """
        Submit the batch syntax checker form with lines ending with tab
        characters.
        """
        variants=['AB026906.1(SDHD):g.7872G>T\t',
                  'AB026906.1(SDHD):g.7872G>T\t',
                  'AB026906.1(SDHD):g.7872G>T\t']
        self._batch('syntax-checker',
                    file='\n'.join(variants),
                    size=len(variants) * 2,
                    lines=len(variants),
                    header='Input\tStatus')

    @fix(database)
    def test_batch_syntaxchecker(self):
        """
        Submit the batch syntax checker form.
        """
        variants = ['AB026906.1(SDHD):g.7872G>T',
                    'NM_003002.1:c.3_4insG',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        self._batch('syntax-checker',
                    file='\n'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    @fix(database, hg19, hg19_transcript_mappings)
    def test_batch_positionconverter(self):
        """
        Submit the batch position converter form.
        """
        variants = ['NM_003002.2:c.204C>T',
                    'NC_000011.9:g.111959625C>T']
        self._batch('position-converter',
                    assembly_name_or_alias='hg19',
                    file='\n'.join(variants),
                    size=len(variants),
                    header='Input Variant')

    @fix(database)
    def test_batch_syntaxchecker_newlines_unix(self):
        """
        Submit batch syntax checker job with Unix line endings.
        """
        variants = ['AB026906.1(SDHD):g.7872G>T',
                    'NM_003002.1:c.3_4insG',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        self._batch('syntax-checker',
                    file='\n'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    @fix(database)
    def test_batch_syntaxchecker_newlines_mac(self):
        """
        Submit batch syntax checker job with Mac line endings.
        """
        variants = ['AB026906.1(SDHD):g.7872G>T',
                    'NM_003002.1:c.3_4insG',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        self._batch('syntax-checker',
                    file='\r'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    @fix(database)
    def test_batch_syntaxchecker_newlines_windows(self):
        """
        Submit batch syntax checker job with Windows line endings.
        """
        variants = ['AB026906.1(SDHD):g.7872G>T',
                    'NM_003002.1:c.3_4insG',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        self._batch('syntax-checker',
                    file='\r\n'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    @fix(database)
    def test_batch_syntaxchecker_newlines_big_unix(self):
        """
        Submit big batch syntax checker job with Unix line endings.
        """
        samples = ['AB026906.1(SDHD):g.7872G>T',
                   'NM_003002.1:c.3_4insG',
                   'AL449423.14(CDKN2A_v002):c.5_400del']
        variants = []
        # Create 240 variants out of 3 samples
        for i in range(80):
            variants.extend(samples)
        self._batch('syntax-checker',
                    file='\n'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    @fix(database)
    def test_batch_syntaxchecker_newlines_big_mac(self):
        """
        Submit big batch syntax checker job with Mac line endings.
        """
        samples = ['AB026906.1(SDHD):g.7872G>T',
                   'NM_003002.1:c.3_4insG',
                   'AL449423.14(CDKN2A_v002):c.5_400del']
        variants = []
        # Create 240 variants out of 3 samples
        for i in range(80):
            variants.extend(samples)
        self._batch('syntax-checker',
                    file='\r'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    @fix(database)
    def test_batch_syntaxchecker_newlines_big_windows(self):
        """
        Submit big batch syntax checker job with Windows line endings.
        """
        samples = ['AB026906.1(SDHD):g.7872G>T',
                   'NM_003002.1:c.3_4insG',
                   'AL449423.14(CDKN2A_v002):c.5_400del']
        variants = []
        # Create 240 variants out of 3 samples
        for i in range(80):
            variants.extend(samples)
        self._batch('syntax-checker',
                    file='\r\n'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    @fix(database)
    def test_batch_syntaxchecker_oldstyle(self):
        """
        Submit the batch syntax checker form with old style input file.
        """
        variants = ['AccNo\tGenesymbol\tMutation',
                    'AB026906.1\tSDHD\tg.7872G>T',
                    'NM_003002.1\t\tc.3_4insG',
                    'AL449423.14\tCDKN2A_v002\tc.5_400del']
        self._batch('syntax-checker',
                    file='\n'.join(variants),
                    size=len(variants)-1,
                    header='Input\tStatus')

    @fix(database, cache('AB026906.1'))
    def test_batch_namechecker_restriction_sites(self):
        """
        Submit the batch name checker form and see if restriction site effects
        are added.
        """
        variants=['AB026906.1:c.274G>T',
                  'AB026906.1:c.[274G>T;143A>G;15G>T]']
        results = self._batch('name-checker',
                              file='\n'.join(variants),
                              size=len(variants),
                              header='Input\tErrors and warnings').strip().split('\n')
        assert 'Restriction Sites Created\tRestriction Sites Deleted' in results[0]
        assert 'CviQI,RsaI\tBccI' in results[1]
        assert 'CviQI,RsaI;HhaI,HinP1I;SfcI\tBccI;;BpmI,BsaXI (2),LpnPI,MnlI' in results[2]

    @fix(database)
    def test_batch_multicolumn(self):
        """
        Submit the batch syntax checker with a multiple-colums input file.

        This by the way also tests for the correct order of batch results.
        """
        variants = [('AB026906.1(SDHD):g.7872G>T', 'NM_003002.1:c.3_4insG'),
                    ('NM_003002.1:c.3_4insG', 'AB026906.1(SDHD):g.7872G>T'),
                    ('AL449423.14(CDKN2A_v002):c.5_400del', 'AL449423.14(CDKN2A_v002):c.5_400del')]
        result = self._batch('syntax-checker',
                             file='\n'.join(['\t'.join(r) for r in variants]),
                             size=len(variants) * 2,
                             header='Input\tStatus',
                             lines=len(variants))
        for line in result.splitlines()[1:]:
            assert len(line.split('\t')) == len(variants[0]) * 2

    def test_download_py(self):
        """
        Download a Python example client for the web service.
        """
        r = self.app.get('/downloads/client-suds.py')
        assert 'text/plain' in r.headers['Content-Type']
        assert '#!/usr/bin/env python' in r.data

    def test_download_rb(self):
        """
        Download a Ruby example client for the web service.
        """
        r = self.app.get('/downloads/client-savon.rb')
        assert 'text/plain' in r.headers['Content-Type']
        assert '#!/usr/bin/env ruby' in r.data

    def test_download_cs(self):
        """
        Download a C# example client for the web service.
        """
        r = self.app.get('/downloads/client-mono.cs')
        assert 'text/plain' in r.headers['Content-Type']
        assert 'public static void Main(String [] args) {' in r.data

    def test_download_php(self):
        """
        Download a PHP example client for the web service.
        """
        r = self.app.get('/downloads/client-php.php')
        assert 'text/plain' in r.headers['Content-Type']
        assert '<?php' in r.data

    def test_downloads_batchtest(self):
        """
        Download the batch test example file.
        """
        r = self.app.get('/downloads/batchtestnew.txt')
        assert 'text/plain' in r.headers['Content-Type']
        assert 'NM_003002.1:c.3_4insG' in r.data

    def test_annotated_soap_api(self):
        """
        Test the SOAP documentation generated from the WSDL.
        """
        r = self.app.get('/soap-api')
        assert 'text/html' in r.headers['Content-Type']
        assert 'Web Service: Mutalyzer' in r.data

    @fix(database, cache('NG_012337.1'))
    def test_getgs(self):
        """
        Test the /getGS interface used by LOVD2.
        """
        r = self.app.get('/getGS',
                         query_string={'variantRecord': 'NM_003002.2',
                                       'forward': '1',
                                       'mutationName': 'NG_012337.1:g.7055C>T'},
                         follow_redirects=True)
        assert '0 Errors' in r.data
        assert '0 Warnings' in r.data
        assert 'Raw variant 1: substitution at 7055' in r.data
        assert 'go to bottom' not in r.data
        assert '<input' not in r.data

    @fix(database, cache('NG_012337.1'))
    def test_getgs_coding_multiple_transcripts(self):
        """
        Test the /getGS interface on a coding description and genomic
        reference with multiple transcripts.
        """
        r = self.app.get('/getGS',
                         query_string={'variantRecord': 'NM_003002.2',
                                       'forward': '1',
                                       'mutationName': 'NG_012337.1:c.45A>T'},
                         follow_redirects=False)
        assert '/name-checker?' in r.location
        assert 'description=NG_012337.1' in r.location

    @fix(database, cache('NG_008939.1'))
    def test_getgs_variant_error(self):
        """
        Test the /getGS interface on a variant description with an error.
        """
        # The error is that position c.45 is a C, not an A.
        r = self.app.get('/getGS',
                         query_string={'variantRecord': 'NM_000532.4',
                                       'forward': '1',
                                       'mutationName': 'NG_008939.1:c.45A>T'},
                         follow_redirects=False)
        assert '/name-checker?' in r.location
        assert 'description=NG_008939.1' in r.location

    @fix(database, hg19, hg19_transcript_mappings)
    def test_variantinfo_g2c(self):
        """
        Test the /Variant_info interface used by LOVD2 (g to c).
        """
        r = self.app.get('/Variant_info',
                         query_string={'LOVD_ver': '2.0-29',
                                       'build': 'hg19',
                                       'acc': 'NM_203473.1',
                                       'var': 'g.48374289_48374389del'})
        assert 'text/plain' in r.headers['Content-Type']
        expected = '\n'.join(['1020', '0', '1072', '48', '48374289', '48374389', 'del'])
        assert r.data == expected

    @fix(database, hg19, hg19_transcript_mappings)
    def test_variantinfo_c2g(self):
        """
        Test the /Variant_info interface used by LOVD2 (c to g).
        """
        r = self.app.get('/Variant_info',
                         query_string={'LOVD_ver': '2.0-29',
                                       'build': 'hg19',
                                       'acc': 'NM_203473.1',
                                       'var': 'c.1020_1072+48del'})
        assert 'text/plain' in r.headers['Content-Type']
        expected = '\n'.join(['1020', '0', '1072', '48', '48374289', '48374389', 'del'])
        assert r.data == expected

    @fix(database, hg19, hg19_transcript_mappings)
    def test_variantinfo_c2g_downstream(self):
        """
        Test the /Variant_info interface used by LOVD2 (c variant downstream
        notation to g).
        """
        r = self.app.get('/Variant_info',
                         query_string={'LOVD_ver': '2.0-29',
                                       'build': 'hg19',
                                       'acc': 'NM_203473.1',
                                       'var': 'c.1709+d187del'})
        assert 'text/plain' in r.headers['Content-Type']
        expected = '\n'.join(['1709', '187', '1709', '187', '48379389', '48379389', 'del'])
        assert r.data == expected

    @fix(database, hg19, hg19_transcript_mappings)
    def test_variantinfo_no_variant(self):
        """
        Test the /Variant_info interface used by LOVD2 (without variant).
        """
        r = self.app.get('/Variant_info',
                         query_string={'LOVD_ver': '2.0-29',
                                       'build': 'hg19',
                                       'acc': 'NM_203473.1'})
        assert 'text/plain' in r.headers['Content-Type']
        assert 'text/plain' in r.content_type
        expected = '\n'.join(['-158', '1709', '1371'])
        assert r.data == expected

    @fix(database, hg19, hg19_transcript_mappings)
    def test_variantinfo_ivs(self):
        """
        Test the /Variant_info interface used by LOVD2 (with IVS positioning).
        """
        r = self.app.get('/Variant_info',
                         query_string={'LOVD_ver': '2.0-33',
                                       'build': 'hg19',
                                       'acc': 'NM_000249.3',
                                       'var': 'c.IVS10+3A>G'})
        assert 'text/plain' in r.headers['Content-Type']
        expected = '\n'.join(['884', '3', '884', '3', '37059093', '37059093', 'subst'])
        assert r.data == expected

    @fix(database)
    def test_upload_local_file(self):
        """
        Test the genbank uploader.
        """
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data',
                            'AB026906.1.gb.bz2')
        r = self.app.post('/reference-loader',
                          data={'method': 'upload',
                                'file': (bz2.BZ2File(path), 'AB026906.1.gb')})
        assert 'Your reference sequence was loaded successfully.' in r.data

        dom = lxml.html.fromstring(r.data)
        reference_url = dom.cssselect('#reference_download')[0].attrib['href']

        r = self.app.get(reference_url)
        assert r.data == bz2.BZ2File(path).read()

    @fix(database)
    def test_upload_local_file_invalid(self):
        """
        Test the genbank uploader with a non-genbank file.
        """
        r = self.app.post('/reference-loader',
                          data={'method': 'upload',
                                'file': (BytesIO('this is not a genbank file'.encode('utf-8')), 'AB026906.1.gb')})
        assert 'Your reference sequence was loaded successfully.' not in r.data
        assert 'The file could not be parsed.' in r.data

    @fix(database, cache('NM_002001.2'))
    def test_reference(self):
        """
        Test if reference files are cached.
        """
        r = self.app.get('/name-checker',
                         query_string={'description': 'NM_002001.2:g.1del'})
        assert '0 Errors' in r.data

        r = self.app.get('/reference/NM_002001.2.gb')
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data',
                            'NM_002001.2.gb.bz2')
        assert r.data == bz2.BZ2File(path).read()

    @fix(database, cache('NM_002001.2'))
    def test_reference_head(self):
        """
        Test if reference files are cached, by issuing a HEAD request.
        """
        r = self.app.get('/name-checker',
                         query_string={'description': 'NM_002001.2:g.1del'})
        assert '0 Errors' in r.data

        r = self.app.head('/reference/NM_002001.2.gb')
        assert r.status_code == 200

    @fix(database)
    def test_reference_head_none(self):
        """
        Test if non-existing reference files gives a 404 on a HEAD request.
        """
        r = self.app.head('/reference/NM_002001.2.gb')
        assert r.status_code == 404

    @fix(database, hg19, hg19_transcript_mappings, cache('NM_003002.2'))
    def test_bed(self):
        """
        BED track for variant.
        """
        r = self.app.get('/bed',
                         query_string={'description': 'NM_003002.2:c.274G>T'})
        assert 'text/plain' in r.headers['Content-Type']
        assert '\t'.join(['chr11', '111959694', '111959695', '274G>T', '0', '+']) in r.data

    @fix(database, hg19, hg19_transcript_mappings, cache('NM_000132.3'))
    def test_bed_reverse(self):
        """
        BED track for variant on reverse strand.
        """
        r = self.app.get('/bed',
                         query_string={'description': 'NM_000132.3:c.[4374A>T;4380_4381del]'})
        assert 'text/plain' in r.headers['Content-Type']
        assert '\t'.join(['chrX', '154157690', '154157691', '4374A>T', '0', '-']) in r.data
        assert '\t'.join(['chrX', '154157683', '154157685', '4380_4381del', '0', '-']) in r.data

    def test_checksyntax_unicode(self):
        """
        Run check syntax form with an invalid variant description containing
        non-ASCII unicode characters.
        """
        r = self.app.get('/syntax-checker',
                         query_string={'description': 'La Pe\xf1a'})
        body = r.get_data(as_text=True)
        assert 'Fatal' in body
        assert 'The &quot;^&quot; indicates the position where the error occurred' in body
        assert 'Expected W:(0123...) (at char 2), (line:1, col:3)' in body

    @fix(database)
    def test_batch_unicode(self):
        """
        Submit a batch form with non-ASCII unicode characters in the input
        file.
        """
        file = '\n'.join(['\u2026AB026906.1:c.274G>T',
                          '\u2026AL449423.14(CDKN2A_v002):c.5_400del'])
        expected = [['\u2026AB026906.1:c.274G>T',
                     '(grammar): Expected W:(0123...) (at char 0), (line:1, col:1)'],
                    ['\u2026AL449423.14(CDKN2A_v002):c.5_400del',
                     '(grammar): Expected W:(0123...) (at char 0), (line:1, col:1)']]

        data = {'job_type': 'syntax-checker',
                'email': 'test@test.test',
                'file': (BytesIO(file.encode('utf-8')), 'test.txt')}

        r = self.app.post('/batch-jobs',
                          data=data)
        progress_url = '/' + r.location.split('/')[-1]

        assert models.BatchJob.query.first().email == 'test@test.test'

        scheduler = Scheduler.Scheduler()
        scheduler.process()

        r = self.app.get(progress_url)

        dom = lxml.html.fromstring(r.data)
        result_url = dom.cssselect('#ifnot_items_left a')[0].attrib['href']

        r = self.app.get(result_url)
        assert 'text/plain' in r.headers['Content-Type']

        result = r.get_data(as_text=True).strip().split('\n')[1:]
        assert expected == [line.split('\t') for line in result]

    @fix(database)
    def test_batch_unicode_email(self):
        """
        Submit a batch form with non-ASCII unicode characters in the email
        address.
        """
        file = '\n'.join(['AB026906.1:c.274G>T',
                          'AL449423.14(CDKN2A_v002):c.5_400del'])
        expected = [['AB026906.1:c.274G>T',
                     'OK'],
                    ['AL449423.14(CDKN2A_v002):c.5_400del',
                     'OK']]

        data = {'job_type': 'syntax-checker',
                'email': 'pe\xf1a@test.test',
                'file': (BytesIO(file.encode('utf-8')), 'test.txt')}

        r = self.app.post('/batch-jobs',
                          data=data)
        progress_url = '/' + r.location.split('/')[-1]

        assert models.BatchJob.query.first().email == 'pe\xf1a@test.test'

        scheduler = Scheduler.Scheduler()
        scheduler.process()

        r = self.app.get(progress_url)

        dom = lxml.html.fromstring(r.data)
        result_url = dom.cssselect('#ifnot_items_left a')[0].attrib['href']

        r = self.app.get(result_url)
        assert 'text/plain' in r.headers['Content-Type']

        result = r.get_data(as_text=True).strip().split('\n')[1:]
        assert expected == [line.split('\t') for line in result]
