"""
Mutalyzer command-line name checker.

.. todo: Refactor this file.
"""


from __future__ import unicode_literals

import argparse
import json

import extractor

from . import _cli_string
from .. import output
from .. import variantchecker


# TODO: This seems like a bit of a weird trick. In any case, we should
#   follow the pattern from the `json` documentation to call
#   `JSONEncoder.default(self, o)` as a fallback.
#
# https://docs.python.org/2/library/json.html#json.JSONEncoder.default
class AlleleEncoder(json.JSONEncoder):
    def default(self, o):
        json_object = o.__dict__
        json_object.update({"hgvs": unicode(o), "weight": o.weight()})

        return json_object
    #default
#MyEncoder


def check_name(description):
    """
    Run the name checker.
    """
    O = output.Output(__file__)

    O.addMessage(__file__, -1, "INFO", "Received variant " + description)

    RD = variantchecker.check_variant(description, O)

    O.addMessage(__file__, -1, "INFO", "Finished processing variant " + description)

    ### OUTPUT BLOCK ###
    gn = O.getOutput("genename")
    if gn :
        print "Gene Name: " + gn[0]
    tv = O.getOutput("transcriptvariant")
    if tv :
        print "Transcript variant: " + tv[0]
        print
    #if

    for i in O.getMessages() :
        print i
    errors, warnings, summary = O.Summary()
    print summary
    print

    if not errors:
        print "Overview of the raw variants:"
        for i in O.getOutput("visualisation"):
            for j in range(len(i)):
                print i[j]
            print
        #for

        print "Genomic description:"
        print O.getIndexedOutput('genomicDescription', 0, '')

        print "\nChromosomal description:"
        print O.getOutput("genomicChromDescription")

        print "\nAffected transcripts:"
        for i in O.getOutput('descriptions'):
            print i
        print "\nAffected proteins:"
        for i in O.getOutput('protDescriptions'):
            print i

        print "\nOld protein:"
        for i in O.getOutput("oldProteinFancyText"):
          print i

        print "\nNew protein:"
        for i in O.getOutput("newProteinFancyText"):
          print i

        print "\nAlternative protein:"
        for i in O.getOutput("altProteinFancyText"):
          print i

        print "\nExon information:"
        for i in O.getOutput("exonInfo") :
            print i

        print "\nCDS  information:"
        print O.getOutput("cdsStart_c"), O.getOutput("cdsStop_c")
        print O.getOutput("cdsStart_g"), O.getOutput("cdsStop_g")

        print "\nEffect on Restriction sites:"
        for i in O.getOutput("restrictionSites") :
            print i

        print "\nLegend:"
        for i in O.getOutput("legends") :
            print i

        reference_sequence = O.getIndexedOutput("original", 0)
        sample_sequence = O.getIndexedOutput("mutated", 0)

        described_allele = extractor.describe_dna(reference_sequence,
                                                  sample_sequence)
        #described_protein_allele = describe.describe(
        #    O.getIndexedOutput("oldProtein", 0),
        #    O.getIndexedOutput("newProtein", 0, default=""),
        #    DNA=False)
        described_protein_allele = ""

        described = described_protein = '(skipped)'

        if described_allele:
            described = described_allele
        if described_protein_allele:
            described_protein = described_protein_allele

        print "\nExperimental services:"
        print described
        print described_protein
        #print "+++ %s" % O.getOutput("myTranscriptDescription")
        print json.dumps({
            #"reference_sequence": reference_sequence,
            #"sample_sequence": sample_sequence,
            "allele_description": described_allele}, cls=AlleleEncoder)


def main():
    """
    Command-line interface to the name checker.
    """
    parser = argparse.ArgumentParser(
        description='Mutalyzer command-line name checker.')
    parser.add_argument(
        'description', metavar='DESCRIPTION', type=_cli_string,
        help='variant description to run the name checker on')

    args = parser.parse_args()
    check_name(args.description)


if __name__ == '__main__':
    main()
