"""
Mutalyzer command-line name checker.

.. todo: Refactor this file.
"""


from __future__ import unicode_literals

import argparse
import sys

from . import _cli_string
from .. import describe
from .. import output
from .. import variantchecker


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

        allele = describe.describe(O.getIndexedOutput("original", 0),
            O.getIndexedOutput("mutated", 0))
        prot_allele = describe.describe(O.getIndexedOutput("oldprotein", 0),
            O.getIndexedOutput("newprotein", 0, default=""), DNA=False)

        extracted = extractedProt = '(skipped)'

        if allele:
            extracted = describe.alleleDescription(allele)
        if prot_allele:
            extractedProt = describe.alleleDescription(prot_allele)

        print "\nExperimental services:"
        print extracted
        print extractedProt
        #print "+++ %s" % O.getOutput("myTranscriptDescription")


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
