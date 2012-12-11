#!/usr/bin/env python
"""
Extract the HGVS variant description from a reference sequence and an observed
sequence.

Usage:
  {command} reference observed

  reference: Reference DNA sequence.
  observed: Observed DNA sequence.

The extracted HGVS description is retrieved from the Mutalyzer SOAP web
service and printed to standard output.
"""


from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main(reference, observed):
    """
    Extract the HGVS variant description from a reference sequence and an
    observed sequence.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.descriptionExtract(reference, observed)

    if result:
        print result.description
    else:
        print 'No description returned.'


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print format_usage()
        sys.exit(1)
    main(*sys.argv[1:])
