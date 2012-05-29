"""
Tests for the SOAP interface to Mutalyzer.
"""


from urllib import urlencode
import simplejson as json
import requests
from nose.tools import *
import mutalyzer


SERVICE_ROOT = 'http://localhost/mutalyzer/json/'


def call(method, **kwargs):
    """
    Do a HTTP/RPC request and decode the JSON response.
    """
    r = requests.get(SERVICE_ROOT + method + '?' + urlencode(kwargs))
    return json.loads(r.content)


class TestServicesJson():
    """
    Test the Mutalyzer HTTP/RPC+JSON interface.
    """
    def test_checksyntax_valid(self):
        """
        Running checkSyntax with a valid variant name should return True.
        """
        r = call('checkSyntax', variant='AB026906.1:c.274G>T')
        assert_equal(r['CheckSyntaxOutput']['valid'], True)

    def test_checksyntax_invalid(self):
        """
        Running checkSyntax with an invalid variant name should return False
        and give at least one error message.
        """
        r = call('checkSyntax', variant='0:abcd')
        assert_equal(r['CheckSyntaxOutput']['valid'], False)
        assert len(r['CheckSyntaxOutput']['messages']['SoapMessage']) >= 1

    def test_checksyntax_empty(self):
        """
        Running checkSyntax with no variant name should raise exception.
        """
        r = call('checkSyntax')
        assert r['Fault']

    def test_transcriptinfo_valid(self):
        """
        Running transcriptInfo with valid arguments should get us a Transcript
        object.
        """
        r = call('transcriptInfo', LOVD_ver='123', build='hg19',
                 accNo='NM_002001.2')
        assert_equal(r['Transcript']['trans_start'], -99)
        assert_equal(r['Transcript']['trans_stop'], 1066)
        assert_equal(r['Transcript']['CDS_stop'], 774)

    def test_info(self):
        """
        Running the info method should give us some version information.
        """
        r = call('info')
        assert_equal(type(r['InfoOutput']['versionParts']['string']), list)
        assert_equal(r['InfoOutput']['version'], mutalyzer.__version__)
