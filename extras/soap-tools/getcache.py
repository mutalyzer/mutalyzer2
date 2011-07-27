#!/usr/bin/env python

# Monkey patch suds, because for some weird reason the location
# http://www.w3.org/2001/xml.xsd is used for the XML namespace, but the W3C
# seems to respond too slow on that url. We use therefore use
# http://www.w3.org/2009/01/xml.xsd which fixes this.
from suds.xsd.sxbasic import Import
_import_open = Import.open
def _import_open_patched(self, *args, **kwargs):
    if self.location == 'http://www.w3.org/2001/xml.xsd':
        self.location = 'http://www.w3.org/2009/01/xml.xsd'
    return _import_open(self, *args, **kwargs)
Import.open = _import_open_patched

import sys
from datetime import datetime, timedelta
from suds.client import Client

URL = 'http://zwarterita/mutalyzer/services/?wsdl'

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
