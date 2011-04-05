#!/usr/bin/env python

import sys
from suds.client import Client

URL = 'http://www.mutalyzer.nl/2.0/services/?wsdl'

if len(sys.argv) < 2:
    print 'Please provide a genomic reference'
    sys.exit(1)

c = Client(URL, cache=None)
o = c.service

# Example: AL449423.14
print 'Getting transcript and info ' + sys.argv[1] + ' ...'

gene = None

if len(sys.argv) > 2:
    print 'Restricting transcripts to ' + sys.argv[2]
    gene = sys.argv[2]

r = o.getTranscriptsAndInfo(sys.argv[1], gene)

if r:
    for t in r.TranscriptInfo:
        print '%s' % t.name
        print
        print 'ID:                 %s' % t.id
        print 'Product:            %s' % t.product
        print 'Locus tag:          %s' % t.locusTag
        print 'Link method:        %s' % t.linkMethod
        print
        print 'Translation start:  %s (c), %s (g)' % (t.cTransStart, t.gTransStart)
        print 'Translation end:    %s (c), %s (g)' % (t.cTransEnd, t.gTransEnd)
        print 'Sortable end:       %s' % t.sortableTransEnd
        print
        print 'Coding start:       %s (c), %s (g)' % (t.cCDSStart, t.gCDSStart)
        print 'Coding stop:        %s (c), %s (g)' % (t.cCDSStop, t.gCDSStop)
        print
        if t.proteinTranscript:
            print 'Protein name:       %s' % t.proteinTranscript.name
            print 'Protein ID:         %s' % t.proteinTranscript.id
            print 'Protein product:    %s' % t.proteinTranscript.product
        else:
            print 'No protein transcript'
        print
        if t.exons:
            for e in t.exons.ExonInfo:
                print 'Exon:               %s - %s (c), %s - %s (g)' % (e.cStart, e.cStop, e.gStart, e.gStop)
        else:
            print 'Huh, no exons?!'
        print
        print
