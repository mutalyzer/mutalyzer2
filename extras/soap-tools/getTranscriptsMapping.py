#!/usr/bin/env python
"""
Get extended information on transcripts contained in a chromosomal region.

Usage:
  {command} build chromosome pos1 pos2 [tolerant]

  build: Human Genome reference build, e.g. 'hg19'.
  chromosome: Chromosome name, e.g. 'chr12'.
  pos1: Beginning of chromosomal region, e.g. 230838269.
  pos2: End of chromosomal region, e.g. 230938269.
  tolerant: Optionally set region matching to be tolerant, e.g. 1.

The transcript information is retrieved from the Mutalyzer SOAP webservice and
printed to standard output.
"""


from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main(build, chromosome, pos1, pos2, tolerant=0):
    """
    Get extended transcript information and print this to standard output.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.getTranscriptsMapping(build, chromosome, int(pos1),
                                           int(pos2), int(tolerant))

    if result:
        for t in result.TranscriptMappingInfo:
            print """Transcript: %s
  Version: %s
  Gene: %s
  Protein: %s
  Orientation: %s
  Start: %s
  Stop: %s
  CDS start: %s
  CDS stop: %s""" % \
            (t.name, t.version, t.gene, t.protein, t.orientation,
             t.start, t.stop, t.cds_start, t.cds_stop)


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print format_usage()
        sys.exit(1)
    main(*sys.argv[1:])
