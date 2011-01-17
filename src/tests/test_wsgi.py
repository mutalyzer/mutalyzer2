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

    def test_checksyntax(self):
        """
        Submit the check syntax form with a valid variant.
        """
        r = self.app.get('/syntaxCheck')
        form = r.forms[0]
        form['variant'] = 'AB026906.1:c.274G>T'
        r = form.submit()
        r.mustcontain('The syntax of this variant is OK!')

if __name__ == '__main__':
    unittest.main()
