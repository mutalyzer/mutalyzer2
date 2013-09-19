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
        assert_equal(r['checkSyntaxResponse']['checkSyntaxResult']['valid'], True)

    def test_checksyntax_invalid(self):
        """
        Running checkSyntax with an invalid variant name should return False
        and give at least one error message.
        """
        r = call('checkSyntax', variant='0:abcd')
        assert_equal(r['checkSyntaxResponse']['checkSyntaxResult']['valid'], False)
        assert len(r['checkSyntaxResponse']['checkSyntaxResult']['messages']['SoapMessage']) >= 1

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
        assert_equal(r['transcriptInfoResponse']['transcriptInfoResult']['trans_start'], -99)
        assert_equal(r['transcriptInfoResponse']['transcriptInfoResult']['trans_stop'], 1066)
        assert_equal(r['transcriptInfoResponse']['transcriptInfoResult']['CDS_stop'], 774)

    def test_info(self):
        """
        Running the info method should give us some version information.
        """
        r = call('info')
        assert_equal(type(r['infoResponse']['infoResult']['versionParts']['string']), list)
        assert_equal(r['infoResponse']['infoResult']['version'], mutalyzer.__version__)

    def test_method_get(self):
        """
        Using a GET request.
        """
        r = requests.get(SERVICE_ROOT + 'checkSyntax', params={'variant': 'AB026906.1:c.274G>T'})
        assert_equal(json.loads(r.content)['checkSyntaxResponse']['checkSyntaxResult']['valid'], True)

    def test_method_get_large(self):
        """
        Using a GET request with too large data.
        """
        r = requests.get(SERVICE_ROOT + 'checkSyntax', params={'variant': 'AB026906.1:c.274G>T', 'dummy': 'a' * 1000000})
        assert '414 Request-URI Too Large' in r.content

    def test_method_post_query_string(self):
        """
        Using a POST request with query string parameters.
        """
        r = requests.post(SERVICE_ROOT + 'checkSyntax', params={'variant': 'AB026906.1:c.274G>T'})
        assert_equal(json.loads(r.content)['checkSyntaxResponse']['checkSyntaxResult']['valid'], True)

    def test_method_post_query_string_large(self):
        """
        Using a POST request with too large data in body.
        """
        r = requests.post(SERVICE_ROOT + 'checkSyntax', params={'variant': 'AB026906.1:c.274G>T', 'dummy': 'a' * 1000000})
        assert '414 Request-URI Too Large' in r.content

    def test_method_post_body(self):
        """
        Using a POST request with data in body.
        """
        r = requests.post(SERVICE_ROOT + 'checkSyntax', data={'variant': 'AB026906.1:c.274G>T'})
        assert_equal(json.loads(r.content)['checkSyntaxResponse']['checkSyntaxResult']['valid'], True)

    def test_method_post_body_large(self):
        """
        Using a POST request with data in body.
        """
        r = requests.post(SERVICE_ROOT + 'checkSyntax', data={'variant': 'AB026906.1:c.274G>T', 'dummy': 'a' * 1000000})
        assert_equal(json.loads(r.content)['checkSyntaxResponse']['checkSyntaxResult']['valid'], True)
