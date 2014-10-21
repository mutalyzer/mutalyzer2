#!/usr/bin/env python
"""
Get static version information from a Mutalyzer installation.

Usage:
  {command}

The version information is retrieved from the Mutalyzer SOAP web service and
printed to standard output.
"""


from __future__ import unicode_literals

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client


from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main():
    """
    Get static version information and print this to standard output.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.info()

    if result:
        print 'Version: %s' % result.version
        print 'Version parts: %s' % ', '.join(
            p for p in result.versionParts.string)
        print 'Release date: %s' % result.releaseDate
        print 'Nomenclature version: %s' % result.nomenclatureVersion
        print 'Nomenclature version parts: %s' % ', '.join(
            p for p in result.nomenclatureVersionParts.string)
        print 'Server name: %s' % result.serverName
        print 'Contact e-mail: %s' % result.contactEmail


if __name__ == '__main__':
    if len(sys.argv) != 1:
        print format_usage()
        sys.exit(1)
    main()
