"""
Base tests for Mutalyzer.
"""


from __future__ import unicode_literals

import urllib2

import pytest


@pytest.mark.skipif(pytest.config.getvalue('allow_network'),
                    reason='--allow-network was specified')
def test_outgoing_fails():
    """
    Internet should not be reachable.
    """
    with pytest.raises(AssertionError):
        urllib2.urlopen('http://www.lumc.nl')
