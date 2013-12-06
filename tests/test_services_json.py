"""
Tests for the JSON interface to Mutalyzer.
"""


from nose.tools import *
import simplejson as json
from spyne.server.null import NullServer
import mutalyzer
from mutalyzer.services.json import application


class TestServicesJson():
    """
    Test the Mutalyzer HTTP/RPC+JSON interface.
    """
    def setUp(self):
        """
        Initialize test server.
        """
        self.server = NullServer(application, ostr=True)

    def _call(self, method, *args, **kwargs):
        r = getattr(self.server.service, method)(*args, **kwargs)
        return json.loads(''.join(r))

    def test_checksyntax_valid(self):
        """
        Running checkSyntax with a valid variant name should return True.
        """
        r = self._call('checkSyntax', variant='AB026906.1:c.274G>T')
        assert_equal(r, {'valid': True, 'messages': []})

    def test_checksyntax_invalid(self):
        """
        Running checkSyntax with an invalid variant name should return False
        and give at least one error message.
        """
        r = self._call('checkSyntax', variant='0:abcd')
        assert_equal(r['valid'], False)
        assert len(r['messages']) >= 1

    def test_checksyntax_empty(self):
        """
        Running checkSyntax with no variant name should raise exception.
        """
        # The validator doesn't work with NullServer, so we cannot do this
        # test. See https://github.com/arskom/spyne/issues/318
        #r = self._call('checkSyntax')
        #assert_equal(r['faultcode'], 'Client.ValidationError')
        pass

    def test_transcriptinfo_valid(self):
        """
        Running transcriptInfo with valid arguments should get us a Transcript
        object.
        """
        r = self._call('transcriptInfo', LOVD_ver='123', build='hg19',
                       accNo='NM_002001.2')
        assert_equal(r['trans_start'], -99)
        assert_equal(r['trans_stop'], 1066)
        assert_equal(r['CDS_stop'], 774)

    def test_info(self):
        """
        Running the info method should give us some version information.
        """
        r = self._call('info')
        assert_equal(type(r['versionParts']), list)
        assert_equal(r['version'], mutalyzer.__version__)
