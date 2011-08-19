#!/usr/bin/env python

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from datetime import datetime, timedelta
from suds.client import Client

URL = 'http://localhost/mutalyzer/services/?wsdl'

c = Client(URL, cache=None)
o = c.service

days = 1
if len(sys.argv) > 1:
    days = int(sys.argv[1])

created_since = datetime.today() - timedelta(days=days)

print 'Getting cache...'

cache = o.getCache(created_since)

if cache:
    for r in cache.CacheEntry:
        print r.name
        if 'gi' in r:
            print 'GI: %s' % r.gi
        print 'Hash: %s' % r.hash
        if 'chromosomeName' in r:
            print r.chromosomeName
        if 'chromosomeStart' in r:
            print r.chromosomeStart
        if 'chromosomeStop' in r:
            print r.chromosomeStop
        if 'chromosomeOrientation' in r:
            print r.chromosomeOrientation
        if 'url' in r:
            print r.url
        print 'Created: %s' % r.created
        if 'cached' in r:
            print 'Cached as %s' % r.cached
        print
