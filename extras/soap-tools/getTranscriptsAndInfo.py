#!/usr/bin/env python
"""
Get extended information on transcripts contained in a genomic reference.

Usage:
  {command} genomic_reference [gene]

  genomic_reference: Genomic reference to look for transcripts in, for
      example 'AL449423.14'.
  gene: Optionally restrict results to transcripts for this gene.

The transcript information is retrieved from the Mutalyzer SOAP web service
and printed to standard output.
"""


from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'


def main(genomic_reference, gene=None):
    """
    Get extended transcript information and print this to standard output.
    """
    service = Client(WSDL_LOCATION, cache=None).service
    result = service.getTranscriptsAndInfo(genomic_reference, gene)

    if result:
        for t in result.TranscriptInfo:
            print """Transcript: %s
  ID: %s
  Product: %s
  Locus tag: %s
  Link method: %s
  Translation:
    Start: %s (c), %s (g), %s (chrom)
    End: %s (c), %s (g), %s (chrom)
    Sortable end: %s
  CDS:
    Start: %s (c), %s (g), %s (chrom)
    End: %s (c), %s (g), %s (chrom)""" % \
            (t.name, t.id, t.product, t.locusTag, t.linkMethod, t.cTransStart,
             t.gTransStart, t.chromTransStart if 'chromTransStart' in t else '-',
             t.cTransEnd, t.gTransEnd, t.chromTransEnd if 'chromTransEnd' in t else '-',
             t.sortableTransEnd, t.cCDSStart, t.gCDSStart,
             t.chromCDSStart if 'chromCDSStart' in t else '-', t.cCDSStop, t.gCDSStop,
             t.chromCDSStop if 'chromCDSStop' in t else '-')

            if 'proteinTranscript' in t:
                print """  Protein:
    Name: %s
    ID: %s
    Product: %s""" % \
                (t.proteinTranscript.name, t.proteinTranscript.id,
                 t.proteinTranscript.product)

            if 'exons' in t:
                print '  Exons:'
                for e in t.exons.ExonInfo:
                    print '    %s - %s (c), %s - %s (g), %s - %s (chrom)' % \
                          (e.cStart, e.cStop, e.gStart, e.gStop,
                           e.chromStart if 'chromStart' in e else '-',
                           e.chromStop if 'chromStop' in e else '-')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print format_usage()
        sys.exit(1)
    main(*sys.argv[1:])
