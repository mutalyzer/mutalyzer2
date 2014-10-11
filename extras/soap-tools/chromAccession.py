#!/usr/bin/env python
"""
Get chromosome accession number from a Mutalyzer installation.

Usage:
  {command} chromosome [build]

  chromosome: Chromosome to get accession number for, e.g. 'chr2'.
  build: Human genome build, 'hg18' or 'hg19' (default: 'hg19').

The accession number is retrieved from the Mutalyzer SOAP web service and
printed to standard output.
"""


from __future__ import unicode_literals

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main(chromosome, build='hg19'):
    """
    Get accession number and print it to standard output.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.chromAccession(build=build, name=chromosome)

    if result:
        print 'Chromosome accession number: %s' % str(result)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print format_usage()
        sys.exit(1)
    main(*sys.argv[1:])
