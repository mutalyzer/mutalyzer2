#!/usr/bin/env python

# Example SOAP client for the Mutalyzer web service in Python using the
# suds library.
#
# See {{ url_for('.webservices', _external=True) }}
#
# Usage:
#   python client-suds.py 'NM_002001.2:c.1del'
#
# This code is in the public domain; it can be used for whatever purpose
# with absolutely no restrictions.

import sys
from suds.client import Client

URL = '{{ soap_wsdl_url }}'

if len(sys.argv) < 2:
    print 'Please provide a variant'
    sys.exit(1)

c = Client(URL, cache=None)
o = c.service

print 'Checking ' + sys.argv[1] + ' ...'

r = o.checkSyntax(sys.argv[1])

if r.valid:
    print 'Valid!'
else:
    print 'Not valid!'

if r.messages:
    for m in r.messages.SoapMessage:
        print 'Message (%s): %s' % (m.errorcode, m.message)
