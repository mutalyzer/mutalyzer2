#!/usr/bin/env python

import sys

import json

import describe

class MyEncoder(json.JSONEncoder):
    def default(self, o):
        return o.__dict__

def main():
    if len(sys.argv) < 3:
        print "usage: " + sys.argv[0] + " reference sample"
        exit()
    #if

    f = open(sys.argv[1], "r")
    ref = f.read()
    f.close()
    f = open(sys.argv[2], "r")
    alt = f.read()
    f.close()

    extracted_allele = describe.describe_dna(ref, alt)

    print "Description Extractor Version " + describe.extractor.VERSION
    #print "HGVS: " + describe.allele_description(extracted_allele)
    print "JSON: " + json.dumps({"reference_sequence": ref, "sample_sequence": alt, "allele_description": extracted_allele}, cls=MyEncoder)

#main

if __name__ == "__main__":
    main()
