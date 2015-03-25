#!/usr/bin/env python
"""
Convert a variant description from c. to g. notation or vice versa.

Usage:
  {command} build description

  build: Human genome reference build to use, i.e. 'hg18' or 'hg19'.
  description: Variant description to convert.

The converted HGVS description(s) is (are) retrieved from the Mutalyzer SOAP
web service and printed to standard output.
"""


from __future__ import unicode_literals

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://127.0.0.1:8081/?wsdl'


def main(build, description):
    """
    Convert variant description from c. to g. notation or vice versa.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.numberConversion(build, description)

    if result:
        for description in result.string:
            print description
    else:
        print 'No descriptions returned.'


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print format_usage()
        sys.exit(1)
    main(*sys.argv[1:])
