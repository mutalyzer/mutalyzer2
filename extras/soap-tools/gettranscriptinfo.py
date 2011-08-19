#!/usr/bin/env python

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds.client import Client

URL = 'http://localhost/mutalyzer/services/?wsdl'

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
