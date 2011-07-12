#!/usr/bin/env python

# Get Mutalyzer static version information from SOAP service.
#
# See http://www.mutalyzer.nl/2.0/webservices
#
# Usage:
#   python info.py
#
# This code is in the public domain; it can be used for whatever purpose
# with absolutely no restrictions.

import sys
from suds.client import Client  # https://fedorahosted.org/suds/

URL = 'http://localhost/mutalyzer/services/?wsdl'

c = Client(URL, cache=None)
o = c.service

print 'Getting version information...'
print

r = o.info()

print 'Version: %s' % r.version
print 'Version parts: %s' % ', '.join(p for p in r.versionParts.string)
print 'Release date: %s' % r.releaseDate
print 'Nomenclature version: %s' % r.nomenclatureVersion
print 'Nomenclature version parts: %s' % ', '.join(p for p in r.nomenclatureVersionParts.string)
print 'Server name: %s' % r.serverName
print 'Contact e-mail: %s' % r.contactEmail
