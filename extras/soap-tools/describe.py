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
from mutalyzer.util import palinsnoop, roll

WSDL_LOCATION = "http://localhost/mutalyzer/services/?wsdl"

class RawVar() :
    def __init__(self, start = 0, start_offset = 0, end = 0, end_offset = 0,
        type = "", deleted = "", inserted = "", hgvs = "=") :
        """
        """

        self.start = start
        self.start_offset = start_offset
        self.end = end
        self.end_offset= end_offset
        self.type = type
        self.deleted = deleted
        self.inserted = inserted
        self.hgvs = hgvs
    #__init__
#RawVar

def LongestCommonSubstring(s1, s2) :
    """
    Find the longest common substring between {s1} and {s2}.

    Mainly copied from:
    http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/
        Longest_common_substring#Python

    @arg s1: String 1.
    @type s1: str
    @arg s2: String 2.
    @type s2: str

    @returns: The end locations and the length of the longest common substring.
    @rtype: tuple(int, int, int)
    """

    len_s1 = len(s1)
    len_s2 = len(s2)
    M = [[0] * (len_s2 + 1) for i in xrange(len_s1 + 1)]
    longest, x_longest, y_longest = 0, 0, 0

    for x in xrange(1, len_s1 + 1) :
        for y in xrange(1, len_s2 + 1) :
            if s1[x - 1] == s2[y - 1] :
                M[x][y] = M[x - 1][y - 1] + 1

                if M[x][y] > longest :
                    longest = M[x][y]
                    x_longest = x
                    y_longest = y
                #if
            #if
            else :
                M[x][y] = 0
        #for
    #for

    return x_longest, y_longest, longest
#LongestCommonSubstring

def DNA_description(s1, s2, s1_start, s1_end, s2_start, s2_end, allele) :
    """
    Give the HGVS description of the change from {s1} to {s2} in the range
    {s1_start}..{s1_end} on {s1} and {s2_start}..{s2_end} on {s2}.

    arg s1: Sequence 1.
    type s1: str
    arg s2: Sequence 2.
    type s2: str
    arg s1_start: Start of the range on {s1}.
    type s1_start: int
    arg s1_end: End of the range on {s1}.
    type s1_end: int
    arg s2_start: Start of the range on {s2}.
    type s2_start: int
    arg s2_end: End of the range on {s2}.
    type s2_end: int

    @returns: HGVS description.
    @rval: str
    """

    # Nothing happened.
    if s1 == s2:
        return "="

    # Insertion / Duplication.
    if s1_start == s1_end :
        ins_length = s2_end - s2_start
        dummy, shift = roll(s2, s2_start + 1, s2_end + 1)

        s1_start += shift + 1
        s1_end += shift + 1
        s2_start += shift + 1
        s2_end += shift + 1

        if s2_start - ins_length >= 0 and \
            s1[s1_start - ins_length:s1_start] == s2[s2_start:s2_end] :

            if ins_length == 1 : # FIXME we forgot the ins of length 1
                hgvs = "%idup" % (s1_start)
                allele.append(RawVar(start = s1_start,
                    type = "dup", hgvs = hgvs))
                return hgvs
            #if

            hgvs = "%i_%idup" % (s1_start - ins_length + 1, s1_end)
            allele.append(RawVar(start = s1_start - ins_length + 1,
                end = s1_end, type = "dup", hgvs = hgvs))
            return hgvs
        #if
        hgvs = "%i_%iins%s" % (s1_start, s1_start + 1, s2[s2_start:s2_end])
        allele.append(RawVar(start = s1_start, end = s1_start + 1,
            inserted = s2[s2_start:s2_end], type = "ins", hgvs = hgvs))
        return hgvs
    #if

    # Deletion.
    if s2_start == s2_end :
        dummy, shift = roll(s1, s1_start + 1, s1_end)

        if s1_start + 1 == s1_end :
            hgvs = "%idel" % (s1_start + shift + 1)
            allele.append(RawVar(start = s1_start + shift + 1, type = "del",
                hgvs = hgvs))
            return hgvs
        #if
        hgvs = "%i_%idel" % (s1_start + shift + 1, s1_end + shift)
        allele.append(RawVar(start = s1_start + shift + 1,
            s1_end = s1_end + shift, type = "del", hgvs = hgvs))
        return hgvs 
    #if

    # Substitution.
    if s1_start + 1 == s1_end and s2_start + 1 == s2_end :
        hgvs = "%i%s>%s" % (s1_start + 1, s1[s1_start], s2[s2_start])
        allele.append(RawVar(start = s1_start + 1, deleted = s1[s1_start],
            inserted = s2[s2_start], type = "subst", hgvs = hgvs))
        return hgvs
    #if

    # Simple InDel.
    if s1_start + 1 == s1_end :
        hgvs = "%idelins%s" % (s1_start + 1, s2[s2_start:s2_end])
        allele.append(RawVar(start = s1_start + 1,
            inserted = s2[s2_start:s2_end], type = "delins", hgvs = hgvs))
        return hgvs
    #if

    # At this stage, we either have an inversion, an indel or a Compound
    # variant.
    s1_end_f, s2_end_f, lcs_f_len = LongestCommonSubstring(s1[s1_start:s1_end],
        s2[s2_start:s2_end])
    s1_end_r, s2_end_r, lcs_r_len = LongestCommonSubstring(s1[s1_start:s1_end],
        Bio.Seq.reverse_complement(s2[s2_start:s2_end]))

    # Palindrome snooping.
    trim = palinsnoop(s1[s1_start + s1_end_r - lcs_r_len:s1_start + s1_end_r])
    if trim < 0 :     # Full palindrome.
        lcs_r_len = 0 # s1_end_r and s2_end_r should not be used after this.

    # Inversion or Compound variant.
    default = "%i_%idelins%s" % (s1_start + 1, s1_end, s2[s2_start:s2_end])

    if not max(lcs_f_len, lcs_r_len) : # Optimisation, not really needed.
        allele.append(RawVar(start = s1_start + 1, end = s1_end,
            inserted = s2[s2_start:s2_end], type = "delins", hgvs = hgvs))
        return default

    # Inversion.
    if lcs_f_len <= lcs_r_len :

        if trim > 0 : # Partial palindrome.
            s1_end_r -= trim
            s2_end_r -= trim
            lcs_r_len -= 2 * trim
        #if

        # Simple Inversion.
        if s2_end - s2_start == lcs_r_len and s1_end - s1_start == lcs_r_len :
            hgvs = "%i_%iinv" % (s1_start + 1, s1_end)
            allele.append(RawVar(start = s1_start + 1, end = s1_end,
                type = "inv", hgvs = hgvs))
            return hgvs
        #if

        r1_len = s1_end_r - lcs_r_len
        r2_len = s1_end - s1_start - s1_end_r
        m1_len = s2_end_r - lcs_r_len
        m2_len = s2_end - s2_start - s2_end_r

        # The flanks of the inversion (but not both) can be empty, so we
        # generate descriptions conditionally.
        leftVariant = ""
        rightVariant = ""
        if r1_len or m2_len :
            lcs = len(longest_common_suffix(s1[s1_start:s1_start + r1_len],
                s2[s2_start:s2_start + m2_len]))
            leftVariant = "%s;" % DNA_description(s1, s2,
                s1_start, s1_start + r1_len - lcs,
                s2_start, s2_start + m2_len - lcs, allele)
        #if
        if r2_len or m1_len :
            lcp = len(longest_common_prefix(s1[s1_end - r2_len:s1_end],
                s2[s2_end - m1_len:s2_end]))
            rightVariant = ";%s" % DNA_description(s1, s2,
                s1_end - r2_len + lcp, s1_end, s2_end - m1_len + lcp, s2_end,
                allele)
        #if

        partial = "%s%i_%iinv%s" % (leftVariant, s1_start + r1_len + 1,
            s1_end - r2_len, rightVariant)
    #if

    # Compound variant.
    else :
        r1_len = s1_end_f - lcs_f_len
        r2_len = s1_end - s1_start - s1_end_f
        m1_len = s2_end_f - lcs_f_len
        m2_len = s2_end - s2_start - s2_end_f

        partial = "%s;%s" % (
            DNA_description(s1, s2,
                s1_start, s1_start + r1_len, s2_start, s2_start + m1_len, 
                allele),
            DNA_description(s1, s2,
                s1_end - r2_len, s1_end, s2_end - m2_len, s2_end, allele))
    #else

    if len(partial) <= len(default) :
        return partial
    return default
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

    if result.rawVariants :
        for i in result.rawVariants.RawVariant :
            print i.description
            print i.visualisation
            print
        #for
    print(result.genomicDescription)

    allele = []
    print(DNA_description(s1, s2, lcp, s1_end, lcp, s2_end, allele))
    for i in allele :
        print i.hgvs
#describe

def main() :
    """
    Main entry point.
    """

    parser = argparse.ArgumentParser(
        prog = "describe",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = "",
        epilog = "")

    parser.add_argument("-d", dest = "description", type = str,
        required = True, help = "HGVS description of a variant.")

    arguments = parser.parse_args()

    describe(arguments.description)
#main

if __name__ == "__main__" :
    main()
