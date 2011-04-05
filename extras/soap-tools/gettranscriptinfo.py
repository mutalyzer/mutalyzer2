#!/usr/bin/env python

import sys
from suds.client import Client

URL = 'http://www.mutalyzer.nl/2.0/services/?wsdl'

if len(sys.argv) < 2:
    print 'Please provide a transcript'
    sys.exit(1)

c = Client(URL, cache=None)
o = c.service

print 'Getting transcript info ' + sys.argv[1] + ' ...'

r = o.transcriptInfo(LOVD_ver="123", build="hg19", accNo=sys.argv[1])
print 'CDS stop: %s' % r.CDS_stop
print 'Trans start: %s' % r.trans_start
print 'Trans stop: %s' % r.trans_stop
