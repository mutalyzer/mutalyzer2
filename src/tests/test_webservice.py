#!/usr/bin/env python
import logging; logging.raiseExceptions = 0
from suds.client import Client, WebFault
import unittest

class TestWebservice(unittest.TestCase):
    """
    Test the Mutalyzer SOAP interface.
    """

    def setUp(self):
        """
        Initialize webservice entrypoint.
        """
        self.client = Client('http://mutalyzer2service/mutalyzer2service?wsdl')

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
        self.assertRaises(WebFault, self.client.service.checkSyntax, '')

if __name__ == '__main__':
    unittest.main()
