#!/usr/bin/env python
"""
Run the Mutalyzer syntaxchecker on a variant description.

Usage:
  {command} description

  description: Variant description to check.

The syntaxchecker results are retrieved from the Mutalyzer SOAP web service
and printed to standard output.
"""


from __future__ import unicode_literals

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://127.0.0.1:8081/?wsdl'


def main(description):
    """
    Run the Mutalyzer syntaxchecker and print results to standard output.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.checkSyntax(description)

    if result.valid:
        print 'Syntax OK!'
    else:
        print 'Syntax NOT OK:'
        for message in result.messages.SoapMessage:
            print message.message


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print format_usage()
        sys.exit(1)
    main(sys.argv[1])
