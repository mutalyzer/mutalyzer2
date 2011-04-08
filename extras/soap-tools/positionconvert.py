#!/usr/bin/env python

# Usage:
#   ./positionconvert.py hg18 'NC_000011.9:g.111959695G>T'

import sys
from suds.client import Client  # https://fedorahosted.org/suds/

URL = 'http://localhost/mutalyzer/services/?wsdl'

if len(sys.argv) < 3:
    print 'Please provide a human genome build and a variant'
    sys.exit(1)

c = Client(URL)
o = c.service

print 'Running position converter for variant ' + sys.argv[1] + ' ...\n'

r = o.numberConversion(sys.argv[1], sys.argv[2])

if r:
    for v in r.string:
        print v
else:
    print 'No variants returned.'
