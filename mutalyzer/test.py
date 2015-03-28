#!/usr/bin/env python

from __future__ import unicode_literals

import json

import describe

class MyEncoder(json.JSONEncoder):
    def default(self, o):
        return o.__dict__

def main():
    ref = "ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGATATTTTT"
    alt = "ACGTCGGTTCGCTAGCTTCGGGGGATAGATAGATATATAGAGATATTTTT"

    extracted_allele = describe.describe_dna(ref, alt)

    print extracted_allele
    print json.dumps({"reference_sequence": ref, "sample_sequence": alt,
        "allele_description": extracted_allele}, cls=MyEncoder)
#main

if __name__ == "__main__":
    main()
