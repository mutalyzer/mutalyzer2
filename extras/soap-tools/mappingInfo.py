#!/usr/bin/env python
"""
Get mapping information from a Mutalyzer installation.

Usage:
  {command} transcript variant [build]

  transcript: Transcript accession number, e.g. 'NM_001008541.1'.
  variant: Variant to map, e.g. 'g.112039014G>T'.
  build: Human genome build, 'hg18' or 'hg19' (default: 'hg19').

The mapping information is retrieved from the Mutalyzer SOAP web service and
printed to standard output.
"""


from __future__ import unicode_literals

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main(transcript, variant, build='hg19'):
    """
    Get mapping information and print it to standard output.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.mappingInfo(LOVD_ver='3.0-beta-06', build=build,
                                 accNo=transcript, variant=variant)

    if result:
        print 'Genomic start: %s' % result.start_g
        print 'Genomic end: %s' % result.end_g
        print 'Transcript start main: %s' % result.startmain
        print 'Transcript start offset: %s' % result.startoffset
        print 'Transcript end main: %s' % result.endmain
        print 'Transcript end offset: %s' % result.endoffset
        print 'Mutation type: %s' % result.mutationType


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print format_usage()
        sys.exit(1)
    main(*sys.argv[1:])
