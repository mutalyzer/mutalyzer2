#!/usr/bin/env python
"""
Get cache entries from a Mutalyzer installation.

Usage:
  {command} days

  days: Retrieve entries of at most this number of days old.

The cache entries are retrieved from the Mutalyzer SOAP webservice and printed
to standard output.
"""


from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from datetime import datetime, timedelta
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main(days):
    """
    Get cache entries and print them to standard output.
    """
    created_since = datetime.today() - timedelta(days=days)
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.getCache(created_since)

    if result:
        for entry in result.CacheEntry:
            print 'Entry %s created at %s' % (entry.name, entry.created)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print format_usage()
        sys.exit(1)
    try:
        main(int(sys.argv[1]))
    except ValueError:
        print 'First argument must be an integer'
        sys.exit(1)
