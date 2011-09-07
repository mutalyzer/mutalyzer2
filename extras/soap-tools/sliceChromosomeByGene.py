#!/usr/bin/env python
"""
Create a slice of a chromosome by gene.

Usage:
  {command} gene

  gene: Gene symbol of the gene to slice.

A slice containing the gene with 5000 upstream bases and 2000 downstream bases
is created with the Mutalyzer SOAP webservice. The resulting UD number is
printed to standard output.
"""


from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client
from suds import WebFault

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main(gene):
    """
    Slice chromosome by gene and print UD number to standard output.
    """
    service = Client(WSDL_LOCATION, cache=None).service

    try:
        print service.sliceChromosomeByGene(gene, 'Human', 5000, 2000)
    except WebFault as message:
        print message


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print format_usage()
        sys.exit(1)
    main(sys.argv[1])
