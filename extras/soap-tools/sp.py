#!/usr/bin/env python

# Example SOAP client for the Mutalyzer web service in Python using the
# SOAPpy library.
#
# See https://www.mutalyzer.nl/webservices
#
# Usage:
#   python sp.py
#
# This code is in the public domain; it can be used for whatever purpose
# with absolutely no restrictions.

from __future__ import unicode_literals

import sys
from SOAPpy import WSDL

o = WSDL.Proxy("http://localhost/mutalyzer/service.wsdl")

# Get all transcripts that are hit when we look at position 159272155 on
# chromosome 1.
print "hg19", "chr1", 159272155
r = o.getTranscripts(build="hg19", chrom="chr1", pos=159272155)
if r:
    # This seems to be a bug in SOAPpy. Arrays of length 1 are
    # flattened, so we cannot iterate over them.
    if not isinstance(r.string, list):
        r.string = [r.string]
    for i in r.string:
        print i, o.getGeneName(build="hg19", accno=i)

# Get all transcripts and genes that have (part of) a transcript in the range
# 159272155-159372155 on chromosome 1.
print "\n", "hg19", "chr1", 159272155, 159372155, 1
r = o.getTranscriptsRange(build="hg19", chrom="chr1", pos1=159272155,
                          pos2=159372155, method=1)
if r:
    # This seems to be a bug in SOAPpy. Arrays of length 1 are
    # flattened, so we cannot iterate over them.
    if not isinstance(r.string, list):
        r.string = [r.string]
    for i in r.string:
        print i, o.getGeneName(build="hg19", accno=i)

# Get all transcripts and genes that have the entire transcript in the range
# 159272155-159372155 on chromosome 1.
print "\n", "hg19", "chr1", 159272155, 159372155, 0
r = o.getTranscriptsRange(build="hg19", chrom="chr1", pos1=159272155,
                          pos2=159372155, method=0)
if r:
    # This seems to be a bug in SOAPpy. Arrays of length 1 are
    # flattened, so we cannot iterate over them.
    if not isinstance(r.string, list):
        r.string = [r.string]
    for i in r.string:
        print i, o.getGeneName(build="hg19", accno=i)

print "\n", "hg19", "NM_002001.2", "c.2del"
r = o.mappingInfo(LOVD_ver="123", build="hg19", accNo="NM_002001.2",
                  variant="c.1del")
print r.mutationType
print r.start_g
print r.end_g

print "\n", "hg19", "NM_002002.2"
r = o.transcriptInfo(LOVD_ver="123", build="hg19", accNo="NM_002001.2")
print r.CDS_stop
print r.trans_start
print r.trans_stop

print "\n", "hg19", "NM_002001.2:c.1del"
r = o.numberConversion(build="hg19", variant="NM_002001.2:c.1del")
if r:
    # This seems to be a bug in SOAPpy. Arrays of length 1 are
    # flattened, so we cannot iterate over them.
    if not isinstance(r.string, list):
        r.string = [r.string]
    for i in r.string:
        print i

print "\n", "hg19", "DMD"
r = o.getTranscriptsByGeneName(build="hg19", name="DMD")
if r:
    # This seems to be a bug in SOAPpy. Arrays of length 1 are
    # flattened, so we cannot iterate over them.
    if not isinstance(r.string, list):
        r.string = [r.string]
    for i in r.string:
        print i

print "\n", "NM_002001.2:g.1del"
r = o.runMutalyzer(variant="NM_002001.2:g.1del")
print r.original
print r.mutated
print r.origMRNA
print r.mutatedMRNA
print r.origCDS
print r.newCDS
print r.origProtein
print r.newProtein
print r.altProtein
print r.errors
print r.warnings
print r.summary
