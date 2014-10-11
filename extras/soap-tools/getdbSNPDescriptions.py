#!/usr/bin/env python
"""
Get HGVS descriptions for a dbSNP rs number.

Usage:
  {command} rs_number

  rs_number: A valid dbSNP rs number, e.g. 'rs9919552'.

The HGVS descriptions are retrieved from the Mutalyzer SOAP web service and
printed to standard output.
"""


from __future__ import unicode_literals

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main(rs_number):
    """
    Get HGVS descriptions and print them to standard output.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.getdbSNPDescriptions(rs_number)

    if result:
        for description in result.string:
            print description


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print format_usage()
        sys.exit(1)
    main(sys.argv[1])
