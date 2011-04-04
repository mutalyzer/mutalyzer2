#!/usr/bin/env python

import sys
from suds.client import Client

URL = 'http://www.mutalyzer.nl/2.0/services/?wsdl'

if len(sys.argv) < 2:
    print 'Please provide a gene name'
    sys.exit(1)

c = Client(URL, cache=None)
o = c.service

print 'Getting transcripts ' + sys.argv[1] + ' ...'

r = o.getTranscriptsByGeneName('hg19', sys.argv[1])

if r:
    for t in r.string:
        print t
