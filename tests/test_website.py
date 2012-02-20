"""
Tests for the WSGI interface to Mutalyzer.

Uses WebTest, see:
  http://pythonpaste.org/webtest/
  http://blog.ianbicking.org/2010/04/02/webtest-http-testing/

I just installed webtest by 'easy_install webtest'.

@todo: Tests for /upload.
"""


#import logging; logging.basicConfig()
import os
import re
import urllib2
import time
import web
from nose.tools import *
from webtest import TestApp
import logging
import urllib


import mutalyzer
from mutalyzer import website
from mutalyzer.util import slow


# TAL logs an awful lot of things with level=DEBUG. On any error, this is all
# dumped to the console, which is very unconvenient. The following suppresses
# most of this.
logging.raiseExceptions = 0
logging.basicConfig(level=logging.INFO)
logging.getLogger('simpleTAL.HTMLTemplateCompiler').setLevel(logging.ERROR)


BATCH_RESULT_URL = 'http://localhost/mutalyzer/Results_{id}.txt'


class TestWSGI():
    """
    Test the Mutalyzer WSGI interface.
    """

    def setUp(self):
        """
        Initialize test application.
        """
        web.config.debug = False
        application = website.app.wsgifunc()
        self.app = TestApp(application)

    def test_root(self):
        """
        Expect the index HTML page.
        """
        r = self.app.get('')
        assert_equal(r.status, '301 Moved Permanently')
        assert r.location.endswith('/')
        r = r.follow()
        assert_equal(r.status, '200 OK')
        # We check for <html> to make sure the menu template is included
        r.mustcontain('<html>',
                      'Welcome to the Mutalyzer web site',
                      '</html>')

    def test_index(self):
        """
        Expect the index HTML page.
        """
        r = self.app.get('/')
        assert_equal(r.status, '200 OK')
        # We check for <html> to make sure the menu template is included
        r.mustcontain('<html>',
                      'Welcome to the Mutalyzer web site',
                      '</html>')

    def test_index_explicit(self):
        """
        Expect the index HTML page.
        """
        r = self.app.get('/index')
        assert_equal(r.status, '200 OK')
        # We check for <html> to make sure the menu template is included
        r.mustcontain('<html>',
                      'Welcome to the Mutalyzer web site',
                      '</html>')

    def test_about(self):
        """
        See if my name is on the About page ;)
        """
        r = self.app.get('/about')
        assert_equal(r.status, '200 OK')
        r.mustcontain('Martijn Vermaat')

    def test_non_existing(self):
        """
        Expect a 404 response.
        """
        r = self.app.get('/this/doesnotexist', status=404)

    def test_menu_links(self):
        """
        Test all links in the main menu.
        """
        ignore = []  # This could contain relative links we want to skip
        r = self.app.get('/')
        for link in r.lxml.cssselect('#menu a'):
            href = link.get('href')
            if href.startswith('http://') or href.startswith('https://') \
            or href in ignore:
                continue
            if not href.startswith('/'):
                href = '/' + href
            self.app.get(href)

    def test_description_extractor(self):
        """
        Submit the variant description extractor.
        """
        r = self.app.get('/descriptionExtract')
        form = r.forms[0]
        form['referenceSeq'] = 'ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA'
        form['variantSeq'] = 'ATGATTTGATCAGATACATGTGATACCGGTAGTTAGGACAA'
        r = form.submit()
        r.mustcontain('g.[5_6insTT;17del;26A&gt;C;35dup]')

    def test_checksyntax_valid(self):
        """
        Submit the check syntax form with a valid variant.
        """
        r = self.app.get('/syntaxCheck')
        form = r.forms[0]
        form['variant'] = 'AB026906.1:c.274G>T'
        r = form.submit()
        r.mustcontain('The syntax of this variant is OK!')

    def test_checksyntax_invalid(self):
        """
        Submit the check syntax form with an invalid variant.
        """
        r = self.app.get('/syntaxCheck')
        form = r.forms[0]
        form['variant'] = 'AB026906.1:c.27'
        r = form.submit()
        r.mustcontain('Fatal',
                      'Details of the parse error')

    def test_check_valid(self):
        """
        Submit the name checker form with a valid variant.
        Should include form and main HTML layout.
        """
        r = self.app.get('/check')
        form = r.forms[0]
        form['mutationName'] = 'NM_002001.2:g.1del'
        r = form.submit()
        r.mustcontain('0 Errors',
                      '0 Warnings',
                      'Raw variant 1: deletion of 1',
                      '<a href="#bottom" class="hornav">go to bottom</a>',
                      '<input value="NM_002001.2:g.1del" type="text" name="mutationName" style="width:100%">')

    def test_check_more_valid(self):
        """
        Test the name checker for some more variants.
        """
        def check_name(name):
            r = self.app.post('/check', {'mutationName': name})
            r.mustcontain('0 Errors')
        names = ['NG_012337.1:g.7055C>T']
        for name in names:
            check_name(name)

    def test_check_invalid(self):
        """
        Submit the name checker form with an invalid variant.
        """
        r = self.app.get('/check')
        form = r.forms[0]
        form['mutationName'] = 'NM_002001.2'
        r = form.submit()
        r.mustcontain('1 Error',
                      '0 Warnings',
                      'Details of the parse error')

    def test_check_protein_reference(self):
        """
        Submit the name checker form with a protein reference sequence (not
        supported).
        """
        r = self.app.get('/check')
        form = r.forms[0]
        form['mutationName'] = 'BAA81889.1:c.274G>T'
        r = form.submit()
        r.mustcontain('1 Error',
                      '0 Warnings',
                      'Protein reference sequences are not supported')

    def test_check_noninteractive(self):
        """
        Submit the name checker form non-interactively.
        Should not include form and main layout HTML.
        """
        r = self.app.get('/check?mutationName=NM_002001.2:g.1del')
        assert_false('<a href="#bottom" class="hornav">go to bottom</a>' in r)
        assert_false('<input value="NM_002001.2:g.1del" type="text" name="mutationName" style="width:100%">' in r)
        r.mustcontain('0 Errors',
                      '0 Warnings',
                      'Raw variant 1: deletion of 1',
                      '<html>',
                      '</html>')

    def test_check_browser_link(self):
        """
        Submit the name checker form with a coding variant on a transcript.
        Should include link to UCSC Genome Browser.
        """
        r = self.app.get('/check')
        form = r.forms[0]
        form['mutationName'] = 'NM_003002.2:c.274G>T'
        r = form.submit()
        bed_track = urllib.quote(r.environ['wsgi.url_scheme'] + '://' + r.environ['HTTP_HOST'] + '/bed?variant=' + urllib.quote('NM_003002.2:c.274G>T'))
        r.mustcontain('<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&amp;position=chr11:111959685-111959705&amp;hgt.customText=%s">View original variant in UCSC Genome Browser</a>' % bed_track)

    def test_checkforward(self):
        """
        A checkForward request should set the given variant in the session and
        redirect to the name checker.
        """
        r = self.app.get('/checkForward?mutationName=NM_002001.2:g.1del')
        assert_equal(r.status, '303 See Other')
        assert r.location.endswith('/check')
        r = r.follow()
        r.mustcontain('0 Errors',
                      '0 Warnings',
                      'Raw variant 1: deletion of 1',
                      '<a href="#bottom" class="hornav">go to bottom</a>',
                      '<input value="NM_002001.2:g.1del" type="text" name="mutationName" style="width:100%">')

    def test_snp_converter_valid(self):
        """
        Submit the SNP converter form with a valid SNP.
        """
        r = self.app.get('/snp')
        form = r.forms[0]
        form['rsId'] = 'rs9919552'
        r = form.submit()
        r.mustcontain('0 Errors',
                      '0 Warnings',
                      'NC_000011.9:g.111959625C>T',
                      'NG_012337.1:g.7055C>T',
                      'NM_003002.2:c.204C>T',
                      'NP_002993.1:p.Ser68=')

    def test_snp_converter_invalid(self):
        """
        Submit the SNP converter form with an invalid SNP.
        """
        r = self.app.get('/snp')
        form = r.forms[0]
        form['rsId'] = 'r9919552'
        r = form.submit()
        r.mustcontain('1 Error',
                      '0 Warnings',
                      'Fatal',
                      'This is not a valid dbSNP id')

    def test_position_converter_c2g(self):
        """
        Submit the position converter form with a valid variant.
        """
        r = self.app.get('/positionConverter')
        form = r.forms[0]
        form['build'] = 'hg19'
        form['variant'] = 'NM_003002.2:c.204C>T'
        r = form.submit()
        r.mustcontain('NC_000011.9:g.111959625C>T')

    def test_position_converter_g2c(self):
        """
        Submit the position converter form with a valid variant.
        """
        r = self.app.get('/positionConverter')
        form = r.forms[0]
        form['build'] = 'hg19'
        form['variant'] = 'NC_000011.9:g.111959625C>T'
        r = form.submit()
        r.mustcontain('NM_003002.2:c.204C>T')

    @slow
    def _batch(self, batch_type='NameChecker', arg1=None, file="", size=0,
               header='', lines=None):
        """
        Submit a batch form.

        @kwarg batch_type: Type of batch job to test. One of NameChecker,
                           SyntaxChecker, PositionConverter.
        @kwarg arg1: Optional extra argument for the batch job.
        @kwarg file: String with variants to use as input for the batch job.
        @kwarg size: Number of variants in input.
        @kwarg header: Message that must be found in the batch job result.
        @kwarg lines: Number of result rows expected.

        @return: The batch result document.
        @rtype: string

        @note: Since the batch files are processed by a running batch daemon
            process, the result gets written to the directory defined by the
            system-wide configuration (e.g. /var/mutalyzer/cache), thus
            inaccessible for the TestApp instance under our current user.
            The 'solution' for this is to download the results via a running
            webserver that should be using the same configuration as the batch
            daemon. Yes, this is a hack.
        """
        r = self.app.get('/batch')
        form = r.forms[0]
        if arg1:
            form['arg1'] = arg1
        form['batchType'] = batch_type
        form['batchEmail'] = 'test@test.test'
        form.set('batchFile', ('test_%s.txt' % batch_type,
                               file))
        r = form.submit()
        id = r.lxml.cssselect('#jobID')[0].get('value')
        max_tries = 60
        for i in range(max_tries):
            r = self.app.get('/progress?jobID=' + id + '&totalJobs=' + str(size) + '&ajax=1')
            assert_equal(r.content_type, 'text/plain')
            #print '%s: %s' % (batch_type, r.body)
            if r.body == 'OK': break
            assert re.match('[0-9]+', r.body)
            time.sleep(2)
        assert_equal(r.body, 'OK')
        # Actually, this only means the last entry was taken from the database
        # queue. It might still be processing, in which case we miss some
        # expected output. So let's wait a few seconds.
        time.sleep(2)
        # This is a hack to get to the batch results (see @note above).
        response = urllib2.urlopen(BATCH_RESULT_URL.format(id=id))
        assert_equal(response.info().getheader('Content-Type'), 'text/plain')
        result = response.read()
        assert header in result
        if not lines:
            lines = size
        assert_equal(len(result.strip().split('\n')) - 1, lines)
        return result

    def test_batch_namechecker(self):
        """
        Submit the batch name checker form.
        """
        variants=['AB026906.1(SDHD):g.7872G>T',
                  'NM_003002.1:c.3_4insG',
                  'AL449423.14(CDKN2A_v002):c.5_400del']
        self._batch('NameChecker',
                    file='\n'.join(variants),
                    size=len(variants),
                    header='Input\tErrors | Messages')

    def test_batch_syntaxchecker(self):
        """
        Submit the batch syntax checker form.
        """
        variants = ['AB026906.1(SDHD):g.7872G>T',
                    'NM_003002.1:c.3_4insG',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        self._batch('SyntaxChecker',
                    file='\n'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    def test_batch_positionconverter(self):
        """
        Submit the batch position converter form.
        """
        variants = ['NM_003002.2:c.204C>T',
                    'NC_000011.9:g.111959625C>T']
        self._batch('PositionConverter',
                    arg1='hg19',
                    file='\n'.join(variants),
                    size=len(variants),
                    header='Input Variant')

    def test_batch_syntaxchecker_newlines_unix(self):
        """
        Submit the batch syntax checker form with unix line endings.
        """
        variants = ['AB026906.1(SDHD):g.7872G>T',
                    'NM_003002.1:c.3_4insG',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        self._batch('SyntaxChecker',
                    file='\n'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    def test_batch_syntaxchecker_newlines_mac(self):
        """
        Submit the batch syntax checker form with mac line endings.
        """
        variants = ['AB026906.1(SDHD):g.7872G>T',
                    'NM_003002.1:c.3_4insG',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        self._batch('SyntaxChecker',
                    file='\r'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    def test_batch_syntaxchecker_newlines_windows(self):
        """
        Submit the batch syntax checker form with windows line endings.
        """
        variants = ['AB026906.1(SDHD):g.7872G>T',
                    'NM_003002.1:c.3_4insG',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        self._batch('SyntaxChecker',
                    file='\r\n'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    def test_batch_syntaxchecker_newlines_big_unix(self):
        """
        Submit the batch syntax checker form with unix line ending
        styles and a big input file.
        """
        samples = ['AB026906.1(SDHD):g.7872G>T',
                   'NM_003002.1:c.3_4insG',
                   'AL449423.14(CDKN2A_v002):c.5_400del']
        variants = []
        # Create 240 variants out of 3 samples
        for i in range(80):
            variants.extend(samples)
        self._batch('SyntaxChecker',
                    file='\n'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    def test_batch_syntaxchecker_newlines_big_mac(self):
        """
        Submit the batch syntax checker form with mac line ending
        styles and a big input file.
        """
        samples = ['AB026906.1(SDHD):g.7872G>T',
                   'NM_003002.1:c.3_4insG',
                   'AL449423.14(CDKN2A_v002):c.5_400del']
        variants = []
        # Create 240 variants out of 3 samples
        for i in range(80):
            variants.extend(samples)
        self._batch('SyntaxChecker',
                    file='\r'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    def test_batch_syntaxchecker_newlines_big_windows(self):
        """
        Submit the batch syntax checker form with windows line ending
        styles and a big input file.
        """
        samples = ['AB026906.1(SDHD):g.7872G>T',
                   'NM_003002.1:c.3_4insG',
                   'AL449423.14(CDKN2A_v002):c.5_400del']
        variants = []
        # Create 240 variants out of 3 samples
        for i in range(80):
            variants.extend(samples)
        self._batch('SyntaxChecker',
                    file='\r\n'.join(variants),
                    size=len(variants),
                    header='Input\tStatus')

    def test_batch_syntaxchecker_oldstyle(self):
        """
        Submit the batch syntax checker form with old style input file.
        """
        variants = ['AccNo\tGenesymbol\tMutation',
                    'AB026906.1\tSDHD\tg.7872G>T',
                    'NM_003002.1\t\tc.3_4insG',
                    'AL449423.14\tCDKN2A_v002\tc.5_400del']
        self._batch('SyntaxChecker',
                    file='\n'.join(variants),
                    size=len(variants)-1,
                    header='Input\tStatus')

    @slow
    def test_batch_namechecker_restriction_sites(self):
        """
        Submit the batch name checker form and see if restriction site effects
        are added.

        Note that we use the @slow decorator here even though it is already
        applied to self._batch. The reason is that we use the result from
        self._batch, which does not exist if @slow checks are disabled.
        """
        variants=['AB026906.1:c.274G>T',
                  'AB026906.1:c.[274G>T;143A>G;15G>T]']
        results = self._batch('NameChecker',
                              file='\n'.join(variants),
                              size=len(variants),
                              header='Input\tErrors | Messages')
        assert 'Restriction Sites Created\tRestriction Sites Deleted' in results[0]
        assert 'CviQI,RsaI\tBccI' in results[1]
        assert 'CviQI,RsaI;HinP1I,HhaI;SfcI\tBccI;;BpmI,MnlI,BsaXI (2),Hpy188III' in results[2]

    @slow
    def test_batch_syntaxchecker_toobig(self):
        """
        Submit the batch syntax checker with a too big input file.
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
        file = seed
        # Very crude way of creating something at least 6MB in size
        while len(file) < 6000000:
            file += file
        r = self.app.get('/batch')
        form = r.forms[0]
        form['batchType'] = 'SyntaxChecker'
        form['batchEmail'] = 'm.vermaat.hg@lumc.nl'
        form.set('batchFile', ('test_batch_toobig.txt',
                               file))
        r = form.submit(status=413)
        assert_equal(r.content_type, 'text/plain')

    @slow
    def test_batch_multicolumn(self):
        """
        Submit the batch syntax checker with a multiple-colums input file.

        This by the way also tests for the correct order of batch results.
        """
        variants = [('AB026906.1(SDHD):g.7872G>T', 'NM_003002.1:c.3_4insG'),
                    ('NM_003002.1:c.3_4insG', 'AB026906.1(SDHD):g.7872G>T'),
                    ('AL449423.14(CDKN2A_v002):c.5_400del', 'AL449423.14(CDKN2A_v002):c.5_400del')]
        result = self._batch('SyntaxChecker',
                             file='\n'.join(['\t'.join(r) for r in variants]),
                             size=len(variants) * 2,
                             header='Input\tStatus',
                             lines=len(variants))
        for line in result.splitlines()[1:]:
            assert_equal(len(line.split('\t')), len(variants[0]) * 2)

    def test_download_py(self):
        """
        Download a Python example client for the webservice.
        """
        r = self.app.get('/download/client-suds.py')
        assert_equal(r.content_type, 'text/plain')
        r.mustcontain('#!/usr/bin/env python')

    def test_download_rb(self):
        """
        Download a Ruby example client for the webservice.
        """
        r = self.app.get('/download/client-savon.rb')
        assert_equal(r.content_type, 'text/plain')
        r.mustcontain('#!/usr/bin/env ruby')

    def test_download_cs(self):
        """
        Download a C# example client for the webservice.
        """
        r = self.app.get('/download/client-mono.cs')
        assert_equal(r.content_type, 'text/plain')
        r.mustcontain('public static void Main(String [] args) {')

    def test_download_php(self):
        """
        Download a PHP example client for the webservice.
        """
        r = self.app.get('/download/client-php.php')
        assert_equal(r.content_type, 'text/plain')
        r.mustcontain('<?php')

    def test_downloads_batchtest(self):
        """
        Download the batch test example file.
        """
        r = self.app.get('/downloads/batchtestnew.txt')
        assert_equal(r.content_type, 'text/plain')
        r.mustcontain('NM_003002.1:c.3_4insG')

    def test_soap_documentation(self):
        """
        Test the SOAP documentation generated from the WSDL.
        """
        r = self.app.get('/documentation')
        assert_equal(r.content_type, 'text/html')
        r.mustcontain('Web Service: Mutalyzer')

    def test_getgs(self):
        """
        Test the /getGS interface used by LOVD2.
        """
        r = self.app.get('/getGS?variantRecord=NM_003002.2&forward=1&mutationName=NG_012337.1:g.7055C%3ET')
        r.mustcontain('0 Errors',
                      '0 Warnings',
                      'Raw variant 1: substitution at 7055')
        assert_equal(r.body.find('go to bottom'), -1)
        assert_equal(r.body.find('<input'), -1)

    def test_variantinfo_g2c(self):
        """
        Test the /Variant_info interface used by LOVD2 (g to c).
        """
        r = self.app.get('/Variant_info?LOVD_ver=2.0-29&build=hg19&acc=NM_203473.1&var=g.48374289_48374389del')
        assert_equal(r.content_type, 'text/plain')
        expected = '\n'.join(['1020', '0', '1072', '48', '48374289', '48374389', 'del'])
        assert_equal(r.body, expected)

    def test_variantinfo_c2g(self):
        """
        Test the /Variant_info interface used by LOVD2 (c to g).
        """
        r = self.app.get('/Variant_info?LOVD_ver=2.0-29&build=hg19&acc=NM_203473.1&var=c.1020_1072%2B48del')
        assert_equal(r.content_type, 'text/plain')
        expected = '\n'.join(['1020', '0', '1072', '48', '48374289', '48374389', 'del'])
        assert_equal(r.body, expected)

    def test_variantinfo_c2g_downstream(self):
        """
        Test the /Variant_info interface used by LOVD2 (c variant downstream
        notation to g).
        """
        r = self.app.get('/Variant_info?LOVD_ver=2.0-29&build=hg19&acc=NM_203473.1&var=c.1709%2Bd187del')
        assert_equal(r.content_type, 'text/plain')
        expected = '\n'.join(['1709', '187', '1709', '187', '48379389', '48379389', 'del'])
        assert_equal(r.body, expected)

    def test_variantinfo_no_variant(self):
        """
        Test the /Variant_info interface used by LOVD2 (without variant).
        """
        r = self.app.get('/Variant_info?LOVD_ver=2.0-32&build=hg19&acc=NM_001083962.1')
        assert_equal(r.content_type, 'text/plain')
        expected = '\n'.join(['-612', '7720', '2016'])
        assert_equal(r.body, expected)

    def test_variantinfo_ivs(self):
        """
        Test the /Variant_info interface used by LOVD2 (with IVS positioning).
        """
        r = self.app.get('/Variant_info?LOVD_ver=2.0-33&build=hg19&acc=NM_000249.3&var=c.IVS10%2B3A%3EG')
        assert_equal(r.content_type, 'text/plain')
        expected = '\n'.join(['884', '3', '884', '3', '37059093', '37059093', 'subst'])
        assert_equal(r.body, expected)

    def test_upload_local_file(self):
        """
        Test the genbank uploader.

        @todo: Use another genbank file to get a UD number and check that
        we can then check variants using that UD number.
        @todo: This genbank file location is bogus. The tests directory is not
            included with the package installation.
        """
        test_genbank_file = os.path.join(os.path.split(mutalyzer.package_root())[0], 'tests/data/AB026906.1.gb')
        r = self.app.get('/upload')
        form = r.forms[0]
        form['invoermethode'] = 'file'
        form.set('bestandsveld', ('test_upload.gb',
                                  open(test_genbank_file, 'r').read()))
        r = form.submit()
        r.mustcontain('Your reference sequence was loaded successfully.')

    def test_upload_local_file_invalid(self):
        """
        Test the genbank uploader with a non-genbank file.

        @note: We add the current time to the file contents to make sure it is
            not recognized by its hash.
        """
        r = self.app.get('/upload')
        form = r.forms[0]
        form['invoermethode'] = 'file'
        form.set('bestandsveld', ('test_upload.gb',
                                  'this is not a genbank file (%s)\n' % time.ctime()))
        r = form.submit()
        r.mustcontain('The file could not be parsed.')
        print r.body

    def test_reference(self):
        """
        Test if reference files are cached.
        """
        r = self.app.get('/check')
        form = r.forms[0]
        form['mutationName'] = 'AB026906.1:c.274G>T'
        r = form.submit()
        r.mustcontain('0 Errors',
                      '1 Warning',
                      'Raw variant 1: substitution at 7872',
                      '<a href="#bottom" class="hornav">go to bottom</a>',
                      '<input value="AB026906.1:c.274G&gt;T" type="text" name="mutationName" style="width:100%">')
        r = self.app.get('/Reference/AB026906.1.gb')
        assert_equal(r.content_type, 'text/plain')
        assert_equal(r.content_length, 26427)
        r.mustcontain('ggaaaaagtc tctcaaaaaa cctgctttat')

    def test_reference_head(self):
        """
        Test if reference files are cached, by issuing a HEAD request.

        Note: The WebTest module also checks that the response to a HEAD
            request is empty, as it should be.
        """
        r = self.app.get('/check')
        form = r.forms[0]
        form['mutationName'] = 'AB026906.1:c.274G>T'
        r = form.submit()
        r.mustcontain('0 Errors',
                      '1 Warning',
                      'Raw variant 1: substitution at 7872',
                      '<a href="#bottom" class="hornav">go to bottom</a>',
                      '<input value="AB026906.1:c.274G&gt;T" type="text" name="mutationName" style="width:100%">')
        r = self.app.head('/Reference/AB026906.1.gb')
        assert_equal(r.content_type, 'text/plain')

    def test_reference_head_none(self):
        """
        Test if non-existing reference files gives a 404 on a HEAD request.
        """
        r = self.app.head('/Reference/AB026906.78.gb', status=404)

    def test_bed(self):
        """
        BED track for variant.
        """
        r = self.app.get('/bed?variant=NM_003002.2%3Ac.274G%3ET')
        assert_equal(r.content_type, 'text/plain')
        r.mustcontain('\t'.join(['chr11', '111959694', '111959695', '274G>T', '0', '+']))

    def test_bed_reverse(self):
        """
        BED track for variant on reverse strand.
        """
        r = self.app.get('/bed?variant=NM_000132.3%3Ac.%5B4374A%3ET%3B4380_4381del%5D')
        assert_equal(r.content_type, 'text/plain')
        r.mustcontain('\t'.join(['chrX', '154157690', '154157691', '4374A>T', '0', '-']))
        r.mustcontain('\t'.join(['chrX', '154157683', '154157685', '4380_4381del', '0', '-']))
