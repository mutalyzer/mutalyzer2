#!/usr/bin/env python

"""
Tests for the SOAP interface to Mutalyzer.
"""

import logging; logging.raiseExceptions = 0
import urllib2
from suds.client import Client
from suds import WebFault
import unittest

WSDL_URL = 'http://mutalyzer.martijn/services/?wsdl'

class TestWSDL(unittest.TestCase):
    """
    Test the Mutalyzer SOAP interface WSDL description.
    """
    def test_wsdl(self):
        """
        Test if the WSDL is available and looks somewhat sensible.
        """
        wsdl = urllib2.urlopen(WSDL_URL).read()
        self.assertTrue(wsdl.startswith("<?xml version='1.0' encoding='UTF-8'?>"))
        self.assertTrue('name="Mutalyzer"' in wsdl)

class TestWebservice(unittest.TestCase):
    """
    Test the Mutalyzer SOAP interface.
    """

    def setUp(self):
        """
        Initialize webservice entrypoint.

        @todo: Start the standalone server and stop it in self.tearDown
        instead of depending on some running instance at a fixed address.
        """
        self.client = Client(WSDL_URL)

    def test_checksyntax_valid(self):
        """
        Running checkSyntax with a valid variant name should say it is valid.
        """
        r = self.client.service.checkSyntax('AB026906.1:c.274G>T')
        self.assertEqual(r, 'The syntax of this variant is OK!')

    def test_checksyntax_invalid(self):
        """
        Running checkSyntax with an invalid variant name should say it is
        invalid.
        """
        r = self.client.service.checkSyntax('0:abcd')
        self.assertEqual(r, 'This variant does not have the right syntax. Please try again.')

    def test_checksyntax_empty(self):
        """
        Running checkSyntax with no variant name should raise exception.
        """
        try:
            self.client.service.checkSyntax()
            self.fail('Expected WebFault exception')
        except WebFault, e:
            self.assertEqual(e.fault.faultstring,
                             'The variant argument is not provided.')

    def test_transcriptinfo_valid(self):
        """
        Running transcriptInfo with valid arguments should get us a Transcript
        object.
        """
        r = self.client.service.transcriptInfo(LOVD_ver='123', build='hg19',
                                               accNo='NM_002001.2')
        self.assertEqual(r.trans_start, -99)
        self.assertEqual(r.trans_stop, 1066)
        self.assertEqual(r.CDS_stop, 774)

    def test_numberconversion_gtoc_valid(self):
        """
        Running numberConversion with valid g variant should give a list of
        c variant names.
        """
        r = self.client.service.numberConversion(build='hg19',
                                                 variant='NC_000001.10:g.159272155del')
        self.assertEqual(type(r), list)
        self.assertTrue('NM_002001.2:c.1del' in r)

    def test_numberconversion_ctog_valid(self):
        """
        Running numberConversion with valid c variant should give a list of
        g variant names.
        """
        r = self.client.service.numberConversion(build='hg19',
                                                 variant='NM_002001.2:c.1del')
        self.assertEqual(type(r), list)
        self.assertTrue('NC_000001.10:g.159272155del' in r)

    def test_gettranscriptsbygenename_valid(self):
        """
        Running getTranscriptsByGeneName with valid gene name should give a
        list of transcripts.
        """
        r = self.client.service.getTranscriptsByGeneName(build='hg19',
                                                         name='DMD')
        self.assertEqual(type(r), list)
        for t in ['NM_004006.2',
                  'NM_000109.3',
                  'NM_004021.2',
                  'NM_004009.3',
                  'NM_004007.2',
                  'NM_004018.2',
                  'NM_004022.2']:
            self.assertTrue(t in r)

if __name__ == '__main__':
    unittest.main()
