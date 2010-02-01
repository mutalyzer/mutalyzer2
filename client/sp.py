#!/usr/bin/python

from SOAPpy import WSDL

o = WSDL.Proxy("http://localhost/mutalyzer2/service.wsdl")

for i in eval(o.getTranscripts("chr1", 159272155)) :
    print i, o.getGeneName(i)
