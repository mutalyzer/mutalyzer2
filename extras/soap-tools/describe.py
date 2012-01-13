#!/usr/bin/python

"""
@requires: sys
@requires: argparse
@requires: Bio.Seq
@requires: suds.client.Client
"""


import sys
import argparse
import Bio.Seq
from suds.client import Client

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()
from mutalyzer.util import longest_common_prefix, longest_common_suffix

WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'

def LongestCommonSubstring(S1, S2) :
    """
    """

    M = [[0] * (len(S2) + 1) for i in xrange(len(S1) + 1)]
    longest, x_longest = 0, 0

    for x in xrange(1, len(S1) + 1) :
        for y in xrange(1, len(S2) + 1) :
            if S1[x - 1] == S2[y - 1] :
                M[x][y] = M[x - 1][y - 1] + 1

                if M[x][y] > longest :
                    longest = M[x][y]
                    x_longest  = x
                #if
            #if
            else :
                M[x][y] = 0
        #for
    #for

    return S1[x_longest - longest:x_longest]
#LongestCommonSubstring

def DNA_description(s1, s2, s1_start, s1_end, s2_start, s2_end) :
    """
    """
    # TODO: Roll the variants to the 3' end.

    # Nothing happened.
    if s1 == s2:
        return '='

    # Insertion / Duplication.
    if s1_start == s1_end :
        ins_length = s2_end - s2_start

        if s2_start - ins_length >= 0 and \
            s1[s1_start - ins_length:s1_start] == s2[s2_start:s2_end] :

            if ins_length == 1 :
                return '%idup' % (s1_start)
            return '%i_%idup' % (s1_start - ins_length + 1, s1_end)
        #if
        return '%i_%iins%s' % (s1_start, s1_start + 1, s2[s2_start:s2_end])
    #if

    # Deletion.
    if s2_start == s2_end :
        if s1_start + 1 == s1_end :
            return '%idel' % (s1_start + 1)
        return '%i_%idel' % (s1_start + 1, s1_end) 
    #if

    # Substitution.
    if s1_start + 1 == s1_end and s2_start + 1 == s2_end :
        return '%i%s>%s' % (s1_start + 1, s1[s1_start], s2[s2_start])

    # Simple InDel.
    if s1_start + 1 == s1_end :
        return '%idelins%s' % (s1_start + 1, s2[s2_start:s2_end])

    # At this stage, we either have an inversion, an indel or a Compound
    # variant.
    lcs_f = LongestCommonSubstring(s1[s1_start:s1_end], s2[s2_start:s2_end])
    lcs_r = LongestCommonSubstring(s1[s1_start:s1_end],
        Bio.Seq.reverse_complement(s2[s2_start:s2_end]))
    lcs_f_len = len(lcs_f)
    lcs_r_len = len(lcs_r)

    # Inversion or Compound variant.
    if max(lcs_f_len, lcs_r_len) > 3 : # TODO: This is not a good criterium.

        # Inversion.
        if lcs_f_len <= lcs_r_len :

            # Simple Inversion.
            if s1_end - s1_start == lcs_r_len :
                return '%i_%iinv' % (s1_start + 1, s1_end)

            # Determine the position of the inversion in both the reference and
            #   the mutated sequence.
            # The extra reverse_complement() is needed because searching for
            #   lcs_f in s2 might give a different result.
            r1_len, r2_len = map(len, s1[s1_start:s1_end].split(lcs_r))
            m1_len, m2_len = map(len, 
                Bio.Seq.reverse_complement(s2[s2_start:s2_end]).split(lcs_r))

            # The flanks of the inversion (but not both) can be empty, so we
            # generate descriptions conditionally.
            leftVariant = ''
            rightVariant = ''
            # TODO: remove the lcs from the leftVariant and the lcp from the
            #   rightVariant.
            if r1_len and m2_len :
                leftVariant = '%s;' % DNA_description(s1, s2,
                    s1_start, s1_start + r1_len, s2_start, s2_start + m2_len)
            if r2_len and m1_len :
                rightVariant = ';%s' % DNA_description(s1, s2,
                    s1_end - r2_len, s1_end, s2_end - m1_len, s2_end)

            return '%s%i_%iinv%s' % (leftVariant,
                s1_start + r1_len + 1, s1_end - r2_len, rightVariant)
        #if

        # Compound variant.
        else :
            # Determine the position of the longest common substring in both
            # the reference and the mutated sequence.
            r1_len, r2_len = map(len, s1[s1_start:s1_end].split(lcs_f))
            m1_len, m2_len = map(len, s2[s2_start:s2_end].split(lcs_f))

            return '%s;%s' % (
                DNA_description(s1, s2,
                    s1_start, s1_start + r1_len, s2_start, s2_start + m1_len),
                DNA_description(s1, s2,
                    s1_end - r2_len, s1_end, s2_end - m2_len, s2_end))
        #else
    #if

    # Default InDel.
    return '%i_%idelins%s' % (s1_start + 1, s1_end, s2[s2_start:s2_end])
#DNA_description

def describe(description) :
    """
    """
    service = Client(WSDL_LOCATION, cache = None).service
    result = service.runMutalyzer(description)

    s1 = str(result.original)
    s2 = str(result.mutated)
    lcp = len(longest_common_prefix(s1, s2))
    lcs = len(longest_common_suffix(s1[lcp:], s2[lcp:]))
    s1_end = len(s1) - lcs
    s2_end = len(s2) - lcs

    print(result.genomicDescription)
    print(DNA_description(s1, s2, lcp, s1_end, lcp, s2_end))
#describe

def main() :
    """
    Main entry point.
    """

    parser = argparse.ArgumentParser(
        prog = 'describe',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = "",
        epilog = "")

    parser.add_argument('-d', dest = 'description', type = str,
        required = True, help = 'HGVS description of a variant.')

    arguments = parser.parse_args()

    describe(arguments.description)
#main

if __name__ == "__main__" :
    main()
