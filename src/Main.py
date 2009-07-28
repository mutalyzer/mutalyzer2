#!/usr/bin/python

import Retriever
import Crossmap

retriever = Retriever.Retriever()
record = retriever.loadrecord("GI16944057")
print record

print "---"

RNAf = [5002, 5125, 27745, 27939, 58661, 58762, 74680, 74767, 103409, 103528, 119465, 119537, 144687, 144810, 148418, 149215]
CDSf = [27925, 27939, 58661, 58762, 74680, 74736]

map = Crossmap.Crossmap(RNAf, CDSf, 1)
print map.test("*31+d100")
print map.g2x(map.test("*31+d100"))
print map.test("*31-u100")
print map.g2x(map.test("*31-u100"))
print map.test("30+100")
print map.g2x(map.test("31+100"))
print map.test("-31-100")
print map.g2x(map.test("-31-100"))
