#!/usr/bin/env python

"""
Tests for the WSGI interface to Mutalyzer.

See:
  http://pythonpaste.org/webtest/
  http://blog.ianbicking.org/2010/04/02/webtest-http-testing/

I just installed webtest by 'easy_install webtest'.
"""

import unittest
from webtest import TestApp

# Todo: Can this be done in a more elegant way?
import site
site.addsitedir('..')
from wsgi import application

class TestWSGI(unittest.TestCase):
    """
    Test the Mutalyzer WSGI interface.
    """

    def setUp(self):
        """
        Initialize test application.
        """
        self.app = TestApp(application)

    def test_index(self):
        """
        Expect the index HTML page.
        """
        r = self.app.get('/')
        self.assertEqual(r.status, '200 OK')
        # We check for <html> to make sure the menu template is included
        r.mustcontain('<html>',
                      'Welcome to the Mutalyzer web site',
                      '</html>')

    def test_non_existing(self):
        """
        Expect a 404 response.
        """
        r = self.app.get('/this/doesnotexist', status=404)

    def test_menu_links(self):
        """
        Test all links in the main menu.
        """
        ignore = ['external',   # Todo: should not be a link
                  'disclaimer', # Todo: add this page
                  'bugtracker']
        r = self.app.get('/')
        for link in r.lxml.cssselect('#menu a'):
            href = link.get('href')
            if href.startswith('http://') or href in ignore:
                continue
            if not href.startswith('/'):
                href = '/' + href
            self.app.get(href)

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

    def test_check_noninteractive(self):
        """
        Submit the name checker form non-interactively.
        Should not include form and main layout HTML.
        """
        r = self.app.get('/check?mutationName=NM_002001.2:g.1del')
        self.assertFalse('<a href="#bottom" class="hornav">go to bottom</a>' in r)
        self.assertFalse('<input value="NM_002001.2:g.1del" type="text" name="mutationName" style="width:100%">' in r)
        r.mustcontain('0 Errors',
                      '0 Warnings',
                      'Raw variant 1: deletion of 1',
                      '<html>',
                      '</html>')

    def test_checkforward(self):
        """
        A checkForward request should set the given variant in the session and
        redirect to the name checker.
        """
        r = self.app.get('/checkForward?mutationName=NM_002001.2:g.1del')
        self.assertEqual(r.status, '303 See Other')
        self.assertTrue(r.location.endswith('/check'))

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
                      'NG_012337.1:g.7055C>T',
                      'NM_003002.2:c.204C>T',
                      'NT_033899.8:g.15522041C>T')

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

    def test_download_py(self):
        """
        Download a Python example client for the webservice.
        """
        r = self.app.get('/download/client-suds.py')
        self.assertEqual(r.content_type, 'text/plain')
        r.mustcontain('#!/usr/bin/env python')

    def test_download_cs(self):
        """
        Download a C# example client for the webservice.
        """
        r = self.app.get('/download/client-mono.cs')
        self.assertEqual(r.content_type, 'text/plain')
        r.mustcontain('public static void Main(string [] args) {')

    def test_downloads_batchtest(self):
        """
        Download the batch test example file.
        """
        r = self.app.get('/downloads/batchtestnew.txt')
        self.assertEqual(r.content_type, 'text/plain')
        r.mustcontain('NM_003002.1:c.3_4insG')

    def test_reference(self):
        """
        Download a reference file.
        """
        r = self.app.get('/Reference/AB026906.1.gb')
        self.assertEqual(r.content_type, 'text/plain')
        self.assertEqual(r.content_length, 26427)
        r.mustcontain('ggaaaaagtc tctcaaaaaa cctgctttat')

    def test_soap_documentation(self):
        """
        Test the SOAP documentation generated from the WSDL.

        @todo: The way this is implemented (the app makes a request to the
               webservice), we cannot test it directly.
       """
        #r = self.app.get('/documentation')
        #self.assertEqual(r.content_type, 'text/html')
        #r.mustcontain('Web Service: MutalyzerService')
        pass

if __name__ == '__main__':
    unittest.main()
