#!/usr/bin/python

from SOAPpy import WSDL

o = WSDL.Proxy("http://<tal tal:replace = "path"></tal>/service.wsdl")

# Get all transcripts that are hit when we look at position 159272155 on 
# chromosome 1.
print "chr1", 159272155
for i in eval(o.getTranscripts("chr1", 159272155)) :
    print i, o.getGeneName(i)

# Get all transcripts and genes that have (part of) a transcript in the range 
# 159272155-159372155 on chromosome 1
print
print "chr1", 159272155, 159372155, 1
for i in eval(o.getTranscriptsRange("chr1", 159272155, 159372155, 1)) :
    print i, o.getGeneName(i)

# Get all transcripts and genes that have the entire transcript in the range 
# 159272155-159372155 on chromosome 1
print
print "chr1", 159272155, 159372155, 0
for i in eval(o.getTranscriptsRange("chr1", 159272155, 159372155, 0)) :
    print i, o.getGeneName(i)

# Get variant information.
print 
print "123", "123", "NM_002001.2", "c.1_12del"
for i in eval(o.varInfo("123", "123", "NM_002001.2", "c.1_12del")) :
    print i,

# Get variant information.
print
print
print "123", "123", "NM_002001.2", ""
for i in eval(o.varInfo("123", "123", "NM_002001.2", "")) :
    print i,
