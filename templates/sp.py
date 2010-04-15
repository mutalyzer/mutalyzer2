#!/usr/bin/python

from SOAPpy import WSDL

o = WSDL.Proxy("http://<tal tal:replace = "path"></tal>/service.wsdl")

# Get all transcripts that are hit when we look at position 159272155 on 
# chromosome 1.
print "hg19", "chr1", 159272155
for i in eval(o.getTranscripts(build = "hg19", chrom = "chr1", 
              pos = 159272155)) :
    print i, o.getGeneName(build = "hg19", accno = i)

# Get all transcripts and genes that have (part of) a transcript in the range 
# 159272155-159372155 on chromosome 1
print
print "hg19", "chr1", 159272155, 159372155, 1
for i in eval(o.getTranscriptsRange(build = "hg19", chrom = "chr1", 
              pos1 = 159272155, pos2 = 159372155, method = 1)) :
    print i, o.getGeneName(build = "hg19", accno = i)

# Get all transcripts and genes that have the entire transcript in the range 
# 159272155-159372155 on chromosome 1
print
print "hg19", "chr1", 159272155, 159372155, 0
for i in eval(o.getTranscriptsRange(build = "hg19", chrom = "chr1", 
              pos1 = 159272155, pos2 = 159372155, method = 0)) :
    print i, o.getGeneName(build = "hg19", accno = i)

# Get variant information.
print 
print "123", "hg19", "NM_002001.2", "c.1_12del"
for i in eval(o.varInfo(LOVD_ver = "123", build = "hg19", 
              accno = "NM_002001.2", var = "c.1_12del")) :
    print i,

# Get variant information.
print
print
print "123", "hg19", "NM_002001.2", "\"\""
for i in eval(o.varInfo(LOVD_ver = "123", build = "hg19", 
              accno = "NM_002001.2", var = "")) :
    print i,

# Get variant information.
print
print 
print "123", "hg18", "NM_002001.2", "c.1_12del"
for i in eval(o.varInfo(LOVD_ver = "123", build = "hg18", 
              accno = "NM_002001.2", var = "c.1_12del")) :
    print i,

# Get variant information.
print
print
print "123", "hg18", "NM_002001.2", "\"\""
for i in eval(o.varInfo(LOVD_ver = "123", build = "hg18", 
              accno = "NM_002001.2", var = "")) :
    print i,
