#!/usr/bin/env python

import sys
from suds.client import Client

URL = 'http://localhost/mutalyzer/services/?wsdl'

if len(sys.argv) < 3:
    print 'Please provide a genomic reference and a transcript reference'
    sys.exit(1)

c = Client(URL, cache=None)
o = c.service

print 'Getting gene and transcript ' + sys.argv[1] + ' / ' + sys.argv[2] + ' ...'

r = o.getGeneAndTranscript(sys.argv[1], sys.argv[2])

print r
