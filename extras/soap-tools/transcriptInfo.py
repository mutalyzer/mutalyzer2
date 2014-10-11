#!/usr/bin/env python
"""
Get transcript information from a Mutalyzer installation.

Usage:
  {command} transcript

  transcript: Transcript accession number, e.g. 'NM_002001.2'.

The transcript information is retrieved from the Mutalyzer SOAP web service
and printed to standard output.
"""


from __future__ import unicode_literals

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main(transcript):
    """
    Get transcript information and print it to standard output.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.transcriptInfo(LOVD_ver='123', build='hg19',
                                    accNo=transcript)

    if result:
        print 'Transcript start: %s' % result.trans_start
        print 'Transcript stop: %s' % result.trans_stop
        print 'CDS stop: %s' % result.CDS_stop


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print format_usage()
        sys.exit(1)
    main(sys.argv[1])
