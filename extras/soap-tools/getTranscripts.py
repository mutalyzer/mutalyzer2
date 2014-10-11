#!/usr/bin/env python
"""
Get transcript accession numbers that overlap with a chromosomal position.

Usage:
  {command} chromosome position

  chromosome: Chromosome to lookup transcripts for (e.g. 'chrX').
  position: Position to lookup overlapping transcripts for.

The transcript accession numbers are retrieved from the Mutalyzer SOAP
web service and printed to standard output.
"""


from __future__ import unicode_literals

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main(chromosome, position):
    """
    Get transcript accession numbers and print them to standard output.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.getTranscripts('hg19', chromosome, position, True)

    if result:
        for transcript in result.string:
            print transcript


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print format_usage()
        sys.exit(1)
    main(*sys.argv[1:])
