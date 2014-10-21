#!/usr/bin/env python
"""
Get transcript and product name from a Mutalyzer installation.

Usage:
  {command} genomic_reference transcript_reference

  genomic_reference: Genomic reference in which to lookup the transcript.
  transcript_reference: Reference of the transcript to lookup.

The transcript and product name are retrieved from the Mutalyzer SOAP
web service and printed to standard output.
"""


from __future__ import unicode_literals

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main(genomic_reference, transcript_reference):
    """
    Get cache entries and print them to standard output.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.getGeneAndTranscript(genomic_reference,
                                          transcript_reference)

    if result:
        print 'Transcript: %s' % result.transcriptName
        print 'Product: %s' % result.productName


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print format_usage()
        sys.exit(1)
    main(*sys.argv[1:])
