#!/usr/bin/env python

# Example SOAP client for the Mutalyzer web service in Python using the
# SOAPpy library.
#
# See {{ url_for('.webservices', _external=True) }}
#
# Usage:
#   python client-soappy.py < textfile
#
# This code is in the public domain; it can be used for whatever purpose
# with absolutely no restrictions.

import sys
from SOAPpy import WSDL

URL = '{{ soap_wsdl_url }}'

o = WSDL.Proxy(URL)

for line in sys.stdin.readlines() :
    for word in line.split() :
        if o.checkSyntax(variant = word).valid == 'true' :
            print word
