#!/usr/bin/python

import logging; logging.basicConfig()
from SOAPpy import WSDL

url = '<tal tal:replace = "path"></tal>/service/?wsdl'
o = WSDL.Proxy(url)

# Get all transcripts that are hit when we look at position 159272155 on
# chromosome 1.
print "hg19", "chr1", 159272155
transcripts = o.getTranscripts(build = "hg19", chrom = "chr1",
                               pos = 159272155)
if not type(transcripts) == list: transcripts = [transcripts]
for i in transcripts :
    print i, o.getGeneName(build = "hg19", accno = i)

# Get all transcripts and genes that have (part of) a transcript in the range
# 159272155-159372155 on chromosome 1
print "\n", "hg19", "chr1", 159272155, 159372155, 1
transcripts = o.getTranscriptsRange(build = "hg19", chrom = "chr1",
                                    pos1 = 159272155, pos2 = 159372155, method = 1)
if not type(transcripts) == list: transcripts = [transcripts]
for i in transcripts :
    print i, o.getGeneName(build = "hg19", accno = i)

# Get all transcripts and genes that have the entire transcript in the range
# 159272155-159372155 on chromosome 1
print "\n", "hg19", "chr1", 159272155, 159372155, 0
transcripts = o.getTranscriptsRange(build = "hg19", chrom = "chr1",
                                    pos1 = 159272155, pos2 = 159372155, method = 0)
if not type(transcripts) == list: transcripts = [transcripts]
for i in transcripts :
    print i, o.getGeneName(build = "hg19", accno = i)

print "\n"
print o.mappingInfo(LOVD_ver = "123", build = "hg19", accNo = "NM_002001.2",
                    variant = "c.1del")

print "\n"
print o.transcriptInfo(LOVD_ver = "123", build = "hg19", accNo = "NM_002001.2")

print "\n"
print o.numberConversion(build = "hg19", variant = "NM_002001.2:c.1del")

m = o.transcriptInfo(LOVD_ver = "123", build = "hg19", accNo = "NM_002001.2")
print m.CDS_stop

transcripts = o.getTranscriptsByGeneName(build = "hg19", name = "DMD")
if not type(transcripts) == list: transcripts = [transcripts]
for i in transcripts :
    print i

mutalyzerOutput = o.runMutalyzer(variant = "NM_002001.2:g.1del")
print mutalyzerOutput.original
print mutalyzerOutput.mutated
print mutalyzerOutput.origMRNA
print mutalyzerOutput.mutatedMRNA
print mutalyzerOutput.origCDS
print mutalyzerOutput.newCDS
print mutalyzerOutput.origProtein
print mutalyzerOutput.newProtein
print mutalyzerOutput.altProtein
print mutalyzerOutput.errors
print mutalyzerOutput.warnings
print mutalyzerOutput.summary
