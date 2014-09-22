"""
Tests for the JSON interface to Mutalyzer.
"""


import simplejson as json
from spyne.server.null import NullServer
import mutalyzer
from mutalyzer import announce
from mutalyzer.services.json import application

from fixtures import database, hg19, hg19_transcript_mappings
from utils import MutalyzerTest
from utils import fix


# Todo: We currently have no way of testing POST requests to the JSON API. We
#     had some tests for this, but they were removed with the new setup [1].
#     They depended on a patch to Spyne [2].
#
# [1] https://git.lumc.nl/mutalyzer/mutalyzer/commit/fbf42d8c4a6c6ec27f9001686941b835867199ff
# [2] https://github.com/LUMC/spyne/commit/58660dec28d47b1c3bf1e46d20f55a913ad036cd


class TestServicesJson(MutalyzerTest):
    """
    Test the Mutalyzer HTTP/RPC+JSON interface.
    """
    def setup(self):
        super(TestServicesJson, self).setup()
        self.server = NullServer(application, ostr=True)

    def _call(self, method, *args, **kwargs):
        r = getattr(self.server.service, method)(*args, **kwargs)
        return json.loads(''.join(r))

    def test_checksyntax_valid(self):
        """
        Running checkSyntax with a valid variant name should return True.
        """
        r = self._call('checkSyntax', variant='AB026906.1:c.274G>T')
        assert r == {'valid': True, 'messages': []}

    def test_checksyntax_invalid(self):
        """
        Running checkSyntax with an invalid variant name should return False
        and give at least one error message.
        """
        r = self._call('checkSyntax', variant='0:abcd')
        assert r['valid'] == False
        assert len(r['messages']) >= 1

    def test_checksyntax_empty(self):
        """
        Running checkSyntax with no variant name should raise exception.
        """
        # The validator doesn't work with NullServer, so we cannot do this
        # test. See https://github.com/arskom/spyne/issues/318
        #r = self._call('checkSyntax')
        #assert r['faultcode'] == 'Client.ValidationError'
        pass

    @fix(database, hg19, hg19_transcript_mappings)
    def test_transcriptinfo_valid(self):
        """
        Running transcriptInfo with valid arguments should get us a Transcript
        object.
        """
        r = self._call('transcriptInfo', LOVD_ver='123', build='hg19',
                       accNo='NM_002001.2')
        assert r['trans_start'] == -99
        assert r['trans_stop'] == 1066
        assert r['CDS_stop'] == 774

    def test_info(self):
        """
        Running the info method should give us some version information.
        """
        r = self._call('info')
        assert type(r['versionParts']) == list
        assert r['version'] == mutalyzer.__version__

    def test_info_announcement(self):
        """
        Running the info method should show us the current announcement
        """
        announce.set_announcement('Test announcement')
        r = self._call('info')
        assert type(r['announcement']) == str
        assert r['announcement'] == 'Test announcement'

        announce.set_announcement('New announcement')
        r = self._call('info')
        assert type(r['announcement']) == str
        assert r['announcement'] == 'New announcement'

        announce.unset_announcement()
        r = self._call('info')
        assert not r.get('announcement')
