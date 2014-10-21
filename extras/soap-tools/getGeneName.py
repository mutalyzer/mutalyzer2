#!/usr/bin/env python
"""
Get gene name from a transcript accession number.

Usage:
  {command} transcript [build]

  transcript: Transcript accession number, e.g. 'NM_002001.2'.
  build: Human genome build, 'hg18' or 'hg19' (default: 'hg19').

The transcript gene name is retrieved from the Mutalyzer SOAP web service and
printed to standard output.
"""


from __future__ import unicode_literals

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main(transcript, build='hg19'):
    """
    Get transcript gene name and print it to standard output.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.getGeneName(build, transcript)

    if result:
        print str(result)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print format_usage()
        sys.exit(1)
    main(*sys.argv[1:])
