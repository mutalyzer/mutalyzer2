#!/usr/bin/env python

# Given a variant in HGVS notation, print raw variants with ASCII-art
# visualisation and alternative descriptions. If VERBOSE is True, also print
# original and mutated sequences.
#
# See http://www.mutalyzer.nl/2.0/webservices
#
# Usage:
#   python namecheck.py 'AB026906.1:c.274delG'
#
# This code is in the public domain; it can be used for whatever purpose
# with absolutely no restrictions.

import sys
from suds.client import Client  # https://fedorahosted.org/suds/

VERBOSE = True

URL = 'http://www.mutalyzer.nl/2.0/services/?wsdl'

if len(sys.argv) < 2:
    print 'Please provide a variant'
    sys.exit(1)

c = Client(URL)
o = c.service

print 'Running name checker for variant ' + sys.argv[1] + ' ...\n'

r = o.runMutalyzer(sys.argv[1])

if r.rawVariants:
    for v in r.rawVariants.RawVariant:
        print 'Raw variant: %s' % v.description
        print '%s\n' % v.visualisation

if VERBOSE:
    print 'Original:\n%s\n' % r.original
    print 'Mutated:\n%s\n' % r.mutated
    print 'origMRNA:\n%s\n' % r.origMRNA
    print 'mutatedMRNA:\n%s\n' % r.mutatedMRNA
    print 'origCDS:\n%s\n' % r.origCDS
    print 'newCDS:\n%s\n' % r.newCDS
    print 'origProtein:\n%s\n' % r.origProtein
    print 'newProtein:\n%s\n' % r.newProtein
    print 'altProtein:\n%s\n' % r.altProtein

print 'Errors: %s' % r.errors
print 'Warnings: %s' % r.warnings
print 'Summary: %s\n' % r.summary

if r.messages:
    for m in r.messages.SoapMessage:
        print 'Error %s: %s\n' % (m.errorcode, m.message)

print 'Chromosomal description: %s' % r.chromDescription
print 'Genomic description: %s' % r.genomicDescription

if r.transcriptDescriptions:
    print 'Affected transcripts:'
    print '\n'.join(r.transcriptDescriptions.string)
if r.proteinDescriptions:
    print 'Affected proteins:'
    print '\n'.join(r.proteinDescriptions.string)
