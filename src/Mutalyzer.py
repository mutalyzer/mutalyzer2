#!/usr/bin/env python

"""
The nomenclature checker.

@todo: Use exceptions for failure handling.
@todo: End vs stop. I guess we should use start/stop (end goes with beginning).
       Or first/last, or acceptor/donor. Anyway, CDS is always denoted with
       start/stop.
       Idea:
       * CDS -> use start/stop
       * splice sites or exons -> acceptor/donor
       * translation -> begin/end
       * any range of bases -> first/last
       * interbase position (if two numbers are used) -> before/after
"""


import sys
import math
from itertools import izip_longest
from operator import itemgetter, attrgetter

import Bio
import Bio.Seq
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import seq3
from Bio import Restriction

from Modules import Retriever
from Modules import GenRecord
from Modules import Crossmap
from Modules import Parser
from Modules import Db
from Modules import Mutator
from Modules import Output
from Modules import Config


##############################################################################
# General utility functions, usefull outside of the Mutalyzer context.
##############################################################################


def _format_range(first, last):
    """
    Simplify a range to one position when applicable.

    @arg first: First coordinate of a range.
    @type first: integer
    @arg last: Second coordinate of a range.
    @type last: integer

    @return: {first}_{last} in case of a real range, {first} otherwise.
    @rtype: string
    """
    if first == last:
        return str(first)

    return '%i_%i' % (first, last)
#_format_range


# Used in: checkDeletionDuplication, checkInsertion
def _roll(s, first, last):
    """
    Determine the variability of a variant by looking at cyclic
    permutations. Not all cyclic permutations are tested at each time, it
    is sufficient to check ``aW'' if ``Wa'' matches (with ``a'' a letter,
    ``W'' a word) when rolling to the left for example.

    @arg s: A reference sequence.
    @type s: string
    @arg first: First position of the pattern in the reference sequence.
    @type first: int
    @arg last: Last position of the pattern in the reference sequence.
    @type last: int

    @return: tuple:
        - left  ; Amount of positions that the pattern can be shifted to
                  the left.
        - right ; Amount of positions that the pattern can be shifted to
                  the right.
    @rtype: tuple(int, int)
    """
    pattern = s[first - 1:last] # Extract the pattern
    pattern_length = len(pattern)

    # Keep rolling to the left as long as a cyclic permutation matches.
    minimum = first - 2
    j = pattern_length - 1
    while minimum > -1 and s[minimum] == pattern[j % pattern_length]:
        j -= 1
        minimum -= 1

    # Keep rolling to the right as long as a cyclic permutation matches.
    maximum = last
    j = 0
    while maximum < len(s) and s[maximum] == pattern[j % pattern_length]:
        j += 1
        maximum += 1

    return first - minimum - 2, maximum - last
#_roll


def _palinsnoop(s):
    """
    Check a sequence for a reverse-complement-palindromic prefix (and
    suffix). If one is detected, return the length of this prefix. If the
    string equals its reverse complement, return -1.

    @arg s: A nucleotide sequence.
    @type s: string

    @return: The number of elements that are palindromic or -1 if the string
             is a 'palindrome'.
    @rtype: string
    """
    s_revcomp = Bio.Seq.reverse_complement(s)

    for i in range(int(math.ceil(len(s) / 2.0))):
        if s[i] != s_revcomp[i]:
            # The first i elements are 'palindromic'.
            return i

    # Perfect 'palindrome'.
    return -1
#_palinsnoop


def _splice(s, splice_sites):
    """
    Construct the transcript or the coding sequence from a record and
    a list of splice sites.

    @arg s: A DNA sequence.
    @type s: string
    @arg splice_sites: A list of even length of integers.
    @type splice_sites: list

    @return: The concatenation of slices from the sequence that is present
             in the GenBank record.
    @rtype: string

    @todo: Assert length of splice_sites is even.
    """
    transcript = ''

    sites_iter = iter(splice_sites)
    for acceptor, donor in izip_longest(sites_iter, sites_iter):
        transcript += s[acceptor - 1:donor]

    return transcript
#_splice


# Todo: refactor
def __nsplice(string, splice_sites, CDS, orientation) :
    """
    Just like _splice(), but it only keeps the parts between CDS[0] and
    CDS[1] (in the right orientation).

    I guess we could easily do this as a separate step after _splice()?

    @todo: keep this function?
    @todo: documentation
    """

    transcript = ""
    if orientation == 1 :
        for i in range(0, len(splice_sites), 2) :
            if CDS[0] >= splice_sites[i] and CDS[0] <= splice_sites[i + 1] :
                transcript += string[CDS[0] - 1:splice_sites[i + 1]]
            else :
                if splice_sites[i] > CDS[0] :
                    transcript += \
                        string[splice_sites[i] - 1:splice_sites[i + 1]]
        #for
    #if
    else :
        for i in range(0, len(splice_sites), 2) :
            if CDS[1] >= splice_sites[i] and CDS[1] <= splice_sites[i + 1] :
                transcript += string[splice_sites[i] - 1:CDS[1]]
            else :
                if splice_sites[i] < CDS[1] :
                    transcript += \
                        string[splice_sites[i] - 1:splice_sites[i + 1]]
        #for
    #else

    return transcript
#__nsplice


def _cds_length(splice_sites):
    """
    Calculate the length of a CDS.

    @arg splice_sites: The coordinates of the CDS including internal splice
                       sites.
    @type splice_sites: list

    @return: Length of the CDS.
    @rtype: int

    @todo: Assert length of splice_sites is even.
    """
    l = 0

    sites_iter = iter(splice_sites)
    for acceptor, donor in izip_longest(sites_iter, sites_iter):
        l += donor - acceptor + 1

    return l
#_cds_length


def _is_dna(s):
    """
    Check whether a string is a DNA string.

    @arg s: Any string or Bio.Seq.Seq instance.
    @type s: string

    @return: True if the string is a DNA string, False otherwise.
    @rtype: boolean
    """
    for i in str(s):
        if not i in IUPAC.unambiguous_dna.letters:
            return False

    return True
#_is_dna


def _longest_common_prefix(s1, s2):
    """
    Calculate the longest common prefix of two strings.

    @arg s1: The first string.
    @type s1: string
    @arg s2: The second string.
    @type s2: string

    @return: The longest common prefix of s1 and s2.
    @rtype: string
    """
    pos = 0

    while pos < min(len(s1), len(s2)) and s1[pos] == s2[pos]:
        pos += 1

    return s1[:pos]
#_longest_common_prefix


def _longest_common_suffix(s1, s2):
    """
    Calculate the longest common suffix of two strings.

    @arg s1: The first string.
    @type s1: string
    @arg s2: The second string.
    @type s2: string

    @return: The longest common suffix of s1 and s2.
    @rtype: string
    """
    return _longest_common_prefix(s1[::-1], s2[::-1])[::-1]
#_longest_common_suffix


def _trim_common(s1, s2) :
    """
    Given two strings, trim their longest common prefix and suffix.

    @arg s1: A string.
    @type s1: string
    @arg s2: Another string.
    @type s2: string

    @return: A tuple of:
        - string: Trimmed version of s1.
        - string: Trimmed version of s2.
        - int:    Length of longest common prefix.
        - int:    Length of longest common suffix.

    @todo: More intelligently handle longest_common_prefix().
    """
    lcp = len(_longest_common_prefix(s1, s2))
    lcs = len(_longest_common_suffix(s1[lcp:], s2[lcp:]))
    return s1[lcp:len(s1) - lcs], s2[lcp:len(s2) - lcs], lcp, lcs
#_trim_common


def _over_splice_site(first, last, splice_sites):
    """
    Check wheter a genomic range {first}_{last} hits a splice site.

    @arg first: The first coordinate of the range in g. notation.
    @type first: int
    @arg last: The last coordinate of the range in g. notation.
    @type last: int
    @arg sites: A list of splice sites in g. notation.
    @type sites: list(int)

    @return: True if one or more splice sites are hit, False otherwise.
    @rtype: boolean
    """
    sites_iter = iter(splice_sites)
    for acceptor, donor in izip_longest(sites_iter, sites_iter):
        if first < acceptor and last >= acceptor:
            return True
        if donor and first <= donor and last > donor:
            return True

    return False
#_over_splice_site


##############################################################################
# End of general utility functions.
##############################################################################


##############################################################################
# Protein specific functionality.
##############################################################################


# Used in: _print_protein_html
def _insert_tag(s, pos1, pos2, tag1, tag2):
    """
    Insert two tags (tag1 and tag2) in string s at positions pos1 and pos2
    respectively if the positions are within the length of s. If not,
    either insert one tag or do nothing. If pos1 equals pos2, don't do
    anything either.

    @arg s: A sequence.
    @type s:
    @arg pos1: Position of tag1.
    @type pos1: int
    @arg pos2: Position of tag2.
    @type pos2: int
    @arg tag1: Content of tag1.
    @type tag1: string
    @arg tag2: Content of tag2.
    @type tag2: string

    @return: The original sequence, or a sequence with eiter tag1, tag2 or
             both tags inserted.
    @rtype: string
    """
    output = s
    block = len(s)

    # Only do something if pos1 != pos2.
    if pos1 != pos2:
        if 0 <= pos1 < block:
            # Insert tag1.
            output = output[:pos1] + tag1 + output[pos1:]
        if 0 <= pos2 < block:
            # Insert tag2.
            output = output[:-(block - pos2)] + tag2 \
                     + output[-(block - pos2):]

    return output
#_insert_tag


def _print_protein_html(s, first, last, O, where):
    """
    Make a fancy representation of a protein and put it in the Output
    object under the name 'where'. The representation contains HTML tags
    and is suitable for viewing in a monospaced font.

    @arg s: A protein sequence.
    @type s: string
    @arg first: First position to highlight.
    @type first: int
    @arg last: Last position to highlight.
    @type last: int
    @arg O: The Output object.
    @type O: Modules.Output.Output
    @arg where: Location in the {O} object to store the representation.
    @type where: string
    """
    if not s: return

    block = 10        # Each block consists of 10 amino acids.
    line = 6 * block  # Each line consists of 6 blocks.

    tag1 = '<b style="color:#FF0000">'  # Use this tag for highlighting.
    tag2 = '</b>'                       # And this one to end highlighting.

    # The maximum length for positions is the 10_log of the length of the
    # protein.
    m = int(math.floor(math.log(len(s), 10)) + 1)
    o = 1

    # Add the first position.
    output = '%s ' % str(o).rjust(m)

    for i in range(0, len(s), block):
        # Add the blocks.
        output += ' ' + _insert_tag(s[i:i + block], first - i, last - i,
                                    tag1, tag2)
        if not (i + block) % line and i + block < len(s):
            # One line done.
            o += line
            O.addOutput(where, output)
            # Add the position (while escaping any potential highlighting).
            output = '<tt style="color:000000;font-weight:normal">%s</tt> ' \
                     % str(o).rjust(m)

    # Add last line.
    O.addOutput(where, output)
#_print_protein_html


def in_frame_description(s1, s2) :
    """
    Give a description of an inframe difference of two proteins. Also give
    the position at which the proteins start to differ and the positions at
    which they are the same again.

    @arg s1: The original protein.
    @type s1: string
    @arg s2: The mutated protein.
    @type s2: string

    @return: A tuple of:
        - string ; Protein description of the change.
        - int    ; First position of the change.
        - int    ; Last position of the change in the first protein.
        - int    ; Last position of the change in the second protein.
    @rtype: tuple(string, int, int, int)

    @todo: More intelligently handle longest_common_prefix().
    """
    if s1 == s2:
        # Nothing happened.
        return ('p.(=)', 0, 0, 0)

    lcp = len(_longest_common_prefix(s1, s2))
    lcs = len(_longest_common_suffix(s1[lcp:], s2[lcp:]))
    s1_end = len(s1) - lcs
    s2_end = len(s2) - lcs

    # Insertion / Duplication / Extention.
    if not s1_end - lcp:
        if len(s1) == lcp:
            return ('p.(*%i%sext*%i)' % \
                    (len(s1) + 1, seq3(s2[len(s1)]), abs(len(s1) - len(s2))),
                    len(s1), len(s1), len(s2))
        inLen = s2_end - lcp

        if lcp - inLen >= 0 and s1[lcp - inLen:lcp] == s2[lcp:s2_end]:
            if inLen == 1:
                return ('p.(%s%idup)' % \
                        (seq3(s1[lcp - inLen]), lcp - inLen + 1),
                        lcp, lcp, lcp + 1)
            return ('p.(%s%i_%s%idup)' % \
                    (seq3(s1[lcp - inLen]),
                     lcp - inLen + 1, seq3(s1[lcp - 1]), lcp),
                    lcp, lcp, lcp + inLen)
        #if
        return ('p.(%s%i_%s%iins%s)' % \
                (seq3(s1[lcp - 1]), lcp, seq3(s1[lcp]),
                 lcp + 1, seq3(s2[lcp:s2_end])),
                lcp, lcp, s2_end)
    #if

    # Deletion / Inframe stop.
    if not s2_end - lcp:
        if len(s2) == lcp:
            return ('p.(%s%i*)' % (seq3(s1[len(s2)]), len(s2) + 1),
                    0, 0, 0)

        if lcp + 1 == s1_end:
            return ('p.(%s%idel)' % (seq3(s1[lcp]), lcp + 1),
                    lcp, lcp + 1, lcp)
        return ('p.(%s%i_%s%idel)' % \
                (seq3(s1[lcp]), lcp + 1, seq3(s1[s1_end - 1]), s1_end),
                lcp, s1_end, lcp)
    #if

    # Substitution.
    if s1_end == s2_end and s1_end == lcp + 1:
        return ('p.(%s%i%s)' % (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp])),
                lcp, lcp + 1, lcp + 1)

    # InDel.
    if lcp + 1 == s1_end:
        return ('p.(%s%idelins%s)' % \
                (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp:s2_end])),
                lcp, lcp + 1, s2_end)
    return ('p.(%s%i_%s%idelins%s)' % \
            (seq3(s1[lcp]), lcp + 1, seq3(s1[s1_end - 1]), s1_end,
             seq3(s2[lcp:s2_end])),
            lcp, s1_end, s2_end)
#in_frame_description


def out_of_frame_description(s1, s2) :
    """
    Give the description of an out of frame difference between two
    proteins. Give a description of an inframe difference of two proteins.
    Also give the position at which the proteins start to differ and the
    end positions (to be compatible with the in_frame_description function).

    @arg s1: The original protein.
    @type s1: string
    @arg s2: The mutated protein.
    @type s2: string

    @return: A tuple of:
        - string ; Protein description of the change.
        - int    ; First position of the change.
        - int    ; Last position of the first protein.
        - int    ; Last position of the second protein.
    @rtype: tuple(string, int, int, int)

    @todo: More intelligently handle longest_common_prefix().
    """
    lcp = len(_longest_common_prefix(s1, s2))

    if lcp == len(s2): # NonSense mutation.
        if lcp == len(s1): # Is this correct?
            return ('p.(=)', 0, 0, 0)
        return ('p.(%s%i*)' % (seq3(s1[lcp]), lcp + 1), lcp, len(s1), lcp)
    if lcp == len(s1) :
        return ('p.(*%i%sext*%i)' % \
                (len(s1) + 1, seq3(s2[len(s1)]), abs(len(s1) - len(s2))),
                len(s1), len(s1), len(s2))
    return ('p.(%s%i%sfs*%i)' % \
            (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp]), len(s2) - lcp + 1),
            lcp, len(s1), len(s2))
#out_of_frame_description


def _protein_description(cds_stop, s1, s2) :
    """
    Wrapper function for the in_frame_description() and
    out_of_frame_description() functions. It uses the value cds_stop to
    decide which one to call.

    @arg cds_stop: Position of the stop codon in c. notation (CDS length).
    @type cds_stop: int
    @arg s1: The original protein.
    @type s1: string
    @arg s2: The mutated protein.
    @type s2: string

    @return: A tuple of:
        - string ; Protein description of the change.
        - int    ; First position of the change.
        - int    ; Last position of the change in the first protein.
        - int    ; Last position of the change in the second protein.
    @rtype: tuple(string, int, int, int)
    """
    if cds_stop % 3:
        description = out_of_frame_description(str(s1), str(s2))
    else:
        description = in_frame_description(str(s1), str(s2))

    if not s2 or str(s1[0]) != str(s2[0]):
        # Mutation in start codon.
        return 'p.?', description[1], description[2], description[3]

    return description
#_protein_description


##############################################################################
# End of protein specific functionality.
##############################################################################


# Used in: _raw_variant
def _is_intronic(loc):
    """
    Check whether a location is intronic.

    @arg loc: A location from the Parser module.
    @type loc: pyparsing.ParseResults

    @return: True if the location is intronic, False otherwise.
    @rtype: boolean
    """
    if not loc:
        return False
    if not loc.PtLoc:
        return False
    if not loc.PtLoc.Offset:
        return False
    return True
#_is_intronic


# Used in: _coding_to_genomic
def _check_intronic_position(main, offset, transcript):
    """
    Check whether a c. position is really in an intron: The main coordinate
    must be a splice site and the offset coordinate must have the correct
    sign.

    @arg main: Main coordinate of the position.
    @type main: integer
    @arg offset: Offset coordinate of the position.
    @type offset: integer
    @arg transcript: Transcript under scrutiny.
    @type transcript: object

    @return: True if the combination (main, offset) is valid for this
             transcript, False otherwise.
    @rtype: boolean

    @todo: Use exceptions.
    """
    main_g = transcript.CM.x2g(main, 0)
    sites = transcript.CM.RNA

    if offset:
        oriented_offset = offset * transcript.CM.orientation
        try:
            i = sites.index(main_g)
            if not i % 2:
                # Splice acceptor, so sign must be -.
                if oriented_offset > 0:
                    return False
            else:
                # Splice donor, so sign must be +.
                if oriented_offset < 0:
                    return False
        except ValueError:
            # The main coordinate is not a splice site.
            return False

    return True
#_check_intronic_position


# Used in: _coding_to_genomic
def _get_offset(location) :
    """
    Convert the offset coordinate in a location (from the Parser) to an
    integer.

    @arg location: A location.
    @type location: pyparsing.ParseResults

    @return: Integer representation of the offset coordinate.
    @rtype: int
    """
    if location.Offset :
        if location.Offset == '?' : # This is highly debatable.
            return 0
        offset = int(location.Offset)
        if location.OffSgn == '-' :
            return -offset
        return offset

    return 0
#_get_offset


def __checkOptArg(ref, p1, p2, arg, O) :
    """
    Do several checks for the optional argument of a variant.


    @arg ref: The reference sequence
    @type ref: string
    @arg p1: Start position of the variant
    @type p1: integer
    @arg p2: End position of the variant
    @type p2: integer
    @arg arg: The optional argument
    @type arg:
    @arg O: The Output object
    @type O: object

    @return: True if the optional argument is correct, False otherwise.
    @rtype: boolean

    @todo: refactor
    @todo: Use exceptions.
    """
    if arg : # The argument is optional, if it is not present, it is correct.
        if arg.isdigit() :         # If it is a digit (3_9del7 for example),
            length = int(arg)      #   the digit must be equal to the length
            interval = p2 - p1 + 1 #   of the given range.
            if length != interval :
                O.addMessage(__file__, 3, "EARGLEN",
                    "The length (%i) differed from that of the range (%i)." % (
                    length, interval))
                return False
            #if
        #if
        else :
            if not _is_dna(arg) : # If it is not a digit, it muse be DNA.
                O.addMessage(__file__, 4, "ENODNA",
                    "Invalid letters in argument.")
                return False
            #if
            # And the DNA must match the reference sequence.
            ref_slice = str(ref[p1 - 1:p2])
            if ref_slice != str(arg) : # FIXME more informative.
                O.addMessage(__file__, 3, "EREF",
                    "%s not found at position %s, found %s instead." % (
                    arg, _format_range(p1, p2), ref_slice))
                return False
            #if
        #else
    #if
    return True
#__checkOptArg


def _add_batch_output(O):
    """
    Format the results to a batch output.

    Filter the mutalyzer output and reformat it for use in the batch system
    as output object 'batchDone'.

    @arg O: The Output object
    @type O: Modules.Output.Output

    @todo: More documentation.
    """
    goi, toi = O.getOutput("geneSymbol")[-1] # Two strings [can be empty]
    tList   = []                             # Temporary List
    tDescr  = []                             # Temporary Descr

    reference = O.getOutput("reference")[-1]
    recordType = O.getOutput("recordType")[0]
    descriptions = O.getOutput("NewDescriptions")
        #iName, jName, mType, cDescr, pDescr, gAcc, cAcc, pAcc,
        #fullDescr, fullpDescr

    if len(descriptions) == 0:
        #No descriptions generated [unlikely]
        return
    if O.Summary()[0]:
        #There were errors during the run, return.
        return
    for descr in descriptions:
        if goi in descr[0] and toi in descr[1]: # Gene and Transcript
            if tDescr:
                # Already inserted a value in the tDescr
                tDescr, tList = [], descriptions
                break
            tDescr = descr

    tList = descriptions

    var = O.getOutput("variant")[-1]

    # Generate output
    outputline = ""
    if tDescr: #Filtering worked, only one Description left
        (gName, trName, mType, cDescr,
            pDescr, gAcc, cAcc, pAcc, fullD, fullpD) = tDescr

        gene = "%s_v%.3i" % (gName, int(trName))

        outputline += "%s\t%s\t%s\t" % (reference, gene, var)

        #Add genomic Description
        outputline += "%s\t" % O.getOutput("gDescription")[0]

        #Add coding Description & protein Description
        outputline += "%s\t%s\t" % (cDescr, pDescr)

        gc = cDescr and "%s:%s" % (gene, cDescr)
        gp = pDescr and "%s:%s" % (gene, pDescr)

        #Add mutation with GeneSymbols
        outputline += "%s\t%s\t" % (gc, gp)

        #Add References, should get genomic ref from parsed data
        if recordType == "LRG":
            gAcc = reference
        if recordType == "GB":
            geno = ["NC", "NG", "AC", "NT", "NW", "NZ", "NS"]
            for g in geno:
                if reference.startswith(g):
                    gAcc = reference
                    break
        outputline += "%s\t%s\t%s\t" % (gAcc or "", cAcc or "", pAcc or "")

    else:
        outputline += "\t"*11

    #Add list of affected transcripts "|" seperator
    if tList:
        outputline += "%s\t" % "|".join(e[-2] for e in tList)
        outputline += "%s\t" % "|".join(e[-1] for e in tList)
    else:
        outputline += "\t"*2

    #Link naar additional info:
    #outputline+="http://localhost/mutalyzer2/redirect?mutationName=%s" %\
    #        "todovariant"


    O.addOutput("batchDone", outputline)
#_add_batch_output


def apply_substitution(position, original, substitute, mutator, record, O):
    """
    Do a semantic check for a substitution, do the actual substitution and
    give it a name.

    @arg position: Genomic location of the substitution.
    @type position: int
    @arg original: Nucleotide in the reference sequence.
    @type original: string
    @arg substitute: Nucleotide in the mutated sequence.
    @type substitute: string
    @arg mutator: A Mutator object.
    @type mutator: Modules.Mutator.Mutator
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg O: The Output object.
    @type O: Modules.Output.Output

    @todo: Exception instead of O.addMessage().
    """
    if not _is_dna(substitute):
        # This is not DNA.
        #O.addMessage(__file__, 4, "ENODNA", "Invalid letter in input")
        # todo: exception
        return

    if original == substitute:
        # This is not a real change.
        O.addMessage(__file__, 3, 'ENOVAR',
                     'No mutation given (%c>%c) at position %i.' % \
                     (original, substitute, position))

    mutator.subM(position, substitute)

    record.name(position, position, 'subst', mutator.orig[position - 1],
                substitute, None)
#apply_substitution


def apply_deletion_duplication(first, last, type, mutator, record, O):
    """
    Do a semantic check for a deletion or duplication, do the actual
    deletion/duplication and give it a name.

    @arg first: Genomic start position of the del/dup.
    @type first: int
    @arg last: Genomic end position of the del/dup.
    @type last: int
    @arg type: The variant type (del or dup).
    @type type: string
    @arg mutator: A Mutator object.
    @type mutator: Modules.Mutator.Mutator
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg O: The Output object.
    @type O: Modules.Output.Output

    @todo: Exception instead of O.addMessage().
    """
    roll = _roll(mutator.orig, first, last)
    shift = roll[1]

    # In the case of RNA, check if we roll over a splice site. If so, make
    # the roll shorter, just up to the splice site.
    if record.record.molType == 'n':
        sites = iter(record.record.geneList[0].transcriptList[0] \
                     .mRNA.positionList)
        for acceptor, donor in izip_longest(sites, sites):
            # Note that acceptor and donor splice sites both point to the
            # first, respectively last, position of the exon, so they are
            # both at different sides of the boundary.
            if last < acceptor and last + roll[1] >= acceptor:
                shift = acceptor - 1 - last
                break
            if last <= donor and last + roll[1] > donor:
                shift = donor - last
                break

    if shift:
        new_first = first + shift
        new_stop = last + shift
        O.addMessage(__file__, 2, 'WROLL',
            'Sequence "%s" at position %s was given, however, ' \
            'the HGVS notation prescribes that it should be "%s" at ' \
            'position %s.' % (
            mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
            _format_range(first, last),
            mutator.visualiseLargeString(str(mutator.orig[new_first - 1:new_stop])),
            _format_range(new_first, new_stop)))

    if shift != roll[1]:
        # The original roll was decreased because it crossed a splice site.
        incorrect_first = first + roll[1]
        incorrect_stop = last + roll[1]
        O.addMessage(__file__, 1, 'IROLLBACK',
            'Sequence "%s" at position %s was not corrected to "%s" at ' \
            'position %s, since they reside in different exons.' % (
            mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
            _format_range(first, last),
            mutator.visualiseLargeString(str(mutator.orig[incorrect_first - 1:incorrect_stop])),
            _format_range(incorrect_first, incorrect_stop)))

    if type == 'del':
        mutator.delM(first, last)
    else :
        mutator.dupM(first, last)

    record.name(first, last, type, '', '', (roll[0], shift))
#apply_deletion_duplication


def apply_inversion(first, last, mutator, record, O) :
    """
    Do a semantic check for an inversion, do the actual inversion, and give
    it a name.

    @arg first: Genomic start position of the inversion.
    @type first: int
    @arg last: Genomic end position of the inversion.
    @type last: int
    @arg mutator: A Mutator object.
    @type mutator: Modules.Mutator.Mutator
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg O: The Output object.
    @type O: Modules.Output.Output

    @todo: Exception instead of O.addMessage().
    """
    snoop = _palinsnoop(mutator.orig[first - 1:last])

    if snoop:
        # We have a reverse-complement-palindromic prefix.
        if snoop == -1 :
            # Actually, not just a prefix, but the entire selected sequence is
            # a 'palindrome'.
            O.addMessage(__file__, 2, 'WNOCHANGE',
                'Sequence "%s" at position %i_%i is a palindrome ' \
                '(its own reverse complement).' % (
                mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
                first, last))
            return
        else:
            O.addMessage(__file__, 2, 'WNOTMINIMAL',
                'Sequence "%s" at position %i_%i is a partial ' \
                'palindrome (the first %i nucleotide(s) are the reverse ' \
                'complement of the last one(s)), the HGVS notation ' \
                'prescribes that it should be "%s" at position %i_%i.' % (
                mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
                first, last, snoop,
                mutator.visualiseLargeString(
                    str(mutator.orig[first + snoop - 1: last - snoop])),
                first + snoop, last - snoop))
            first += snoop
            last -= snoop

    mutator.invM(first, last)

    if first == last:
        O.addMessage(__file__, 2, 'WWRONGTYPE', 'Inversion at position ' \
            '%i is actually a substitution.' % first_g)
        record.name(first, first, 'subst', mutator.orig[first - 1],
            Bio.Seq.reverse_complement(mutator.orig[first - 1]), None)
    else :
        record.name(first, last, 'inv', '', '', None)
#apply_inversion


def apply_insertion(before, after, s, mutator, record, O):
    """
    Do a semantic check for an insertion, do the actual insertion, and give
    it a name.

    @arg before: Genomic position before the insertion.
    @type before: int
    @arg after: Genomic position after the insertion.
    @type after: int
    @arg s: Nucleotides to be inserted.
    @type s: string
    @arg mutator: A Mutator object.
    @type mutator: Modules.Mutator.Mutator
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg O: The Output object.
    @type O: Modules.Output.Output

    @todo: Exception instead of O.addMessage().
    """
    if before + 1 != after:
        O.addMessage(__file__, 3, 'EINSRANGE',
            '%i and %i are not consecutive positions.' % (before, after))
        return

    if not s or not _is_dna(s):
        O.addMessage(__file__, 3, 'EUNKVAR', 'Although the syntax of this ' \
            'variant is correct, the effect can not be analysed.')
        return

    insertion_length = len(s)

    mutator.insM(before, s)
    new_before = mutator.shiftpos(before)
    new_stop = mutator.shiftpos(before) + insertion_length

    roll = _roll(mutator.mutated, new_before + 1, new_stop)
    shift = roll[1]

    # In the case of RNA, check if we roll over a splice site. If so, make
    # the roll shorter, just up to the splice site.
    if record.record.molType == 'n' :
        sites = iter(record.record.geneList[0].transcriptList[0] \
                     .mRNA.positionList)
        for acceptor, donor in izip_longest(sites, sites):
            # Note that acceptor and donor splice sites both point to the
            # first, respectively last, position of the exon, so they are
            # both at different sides of the boundary.
            if new_stop < acceptor and new_stop + roll[1] >= acceptor:
                shift = acceptor - 1 - new_stop
                break
            if new_stop <= donor and new_stop + roll[1] > donor:
                shift = donor - new_stop
                break

    if roll[0] + shift >= insertion_length:
        # Todo: Could there also be a IROLLBACK message in this case?
        O.addMessage(__file__, 2, 'WINSDUP',
            'Insertion of %s at position %i_%i was given, ' \
            'however, the HGVS notation prescribes that it should be a ' \
            'duplication of %s at position %i_%i.' % (
            s, before, before + 1,
            mutator.mutated[new_before + shift:new_stop + shift], before + shift,
            before + shift + insertion_length - 1))
        after += shift - 1
        before = after - insertion_length + 1
        record.name(before, after, 'dup', '', '',
                    (roll[0] + shift - insertion_length, 0))
    else:
        if shift:
            O.addMessage(__file__, 2, 'WROLL', 'Insertion of %s at position ' \
                '%i_%i was given, however, the HGVS notation prescribes ' \
                'that it should be an insertion of %s at position %i_%i.' % (
                s, before, before + 1,
                mutator.mutated[new_before + shift:new_stop + shift],
                new_before + shift, new_before + shift + 1))
        if shift != roll[1]:
            O.addMessage(__file__, 1, 'IROLLBACK',
                'Insertion of %s at position %i_%i was not corrected to an ' \
                'insertion of %s at position %i_%i, since they reside in ' \
                'different exons.' % (
                s, before, before + 1,
                mutator.mutated[new_before + roll[1]:new_stop + roll[1]],
                new_before + roll[1], new_before + roll[1] + 1))
        record.name(before, before + 1, 'ins',
                    mutator.mutated[new_before + shift:new_stop + shift], '',
                    (roll[0], shift))
#apply_insertion


def apply_delins(first, last, del, ins, mutator, record, output):
    """
    Do a semantic check for an delins, do the actual delins, and give
    it a name.

    @arg first: Genomic start position of the delins.
    @type first: int
    @arg last: Genomic end position of the delins.
    @type last: int
    @arg del: Sequence to delete (may be None, in which case it will be
              constructed from the reference sequence).
    @type del: string
    @arg ins: Sequence to insert.
    @type ins: string
    @arg mutator: A Mutator object.
    @type mutator: Modules.Mutator.Mutator
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg output: The Output object.
    @type output: Modules.Output.Output

    @todo: Exception instead of O.addMessage().
    """
    if not del:
        del = mutator.orig[first - 1:last]

    if str(del) == str(ins):
        output.addMessage(__file__, 2, 'WNOCHANGE',
                          'Sequence "%s" at position %i_%i is identical to ' \
                          'the variant.' % (
                mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
                              first, last))
        return

    del_trimmed, ins_trimmed, lcp, lcs = _trim_common(del, ins)

    if not len(del_trimmed):
        output.addMessage(__file__, 2, 'WWRONGTYPE', 'The given DelIns ' \
                          'is actually an insertion.')
        apply_insertion(first + lcp - 1, first + lcp, ins_trimmed, mutator,
                        record, output)
        return

    if len(del_trimmed) == 1 and len(ins_trimmed) == 1:
            output.addMessage(__file__, 2, 'WWRONGTYPE', 'The given DelIns ' \
                              'is actually a substitution.')
            apply_substitution(first + lcp, del_trimmed, ins_trimmed, mutator,
                               record, output)
            return

    if not len(ins_trimmed):
        output.addMessage(__file__, 2, 'WWRONGTYPE', 'The given DelIns ' \
                          'is actually a deletion.')
        apply_deletion_duplication(first + lcp, last - lcs, 'del',
                                   mutator, record, output)
        return

    if str(Bio.Seq.reverse_complement(del_trimmed)) == ins_trimmed:
        output.addMessage(__file__, 2, 'WWRONGTYPE', 'The given DelIns ' \
                          'is actually an inversion.')
        apply_inversion(first + lcp, last - lcs, mutator,
                        record, output)
        return

    if len(ins) != len(ins_trimmed):
        output.addMessage(__file__, 2, 'WNOTMINIMAL',
                'Sequence "%s" at position %i_%i has the same prefix or ' \
                'suffix as the inserted sequence "%s". The HGVS notation ' \
                'prescribes that it should be "%s" at position %i_%i.' % (
                mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
                first, last, ins, ins_trimmed, first + lcp, last - lcs))

    mutator.delinsM(first + lcp, last - lcs, ins_trimmed)

    record.name(first + lcp, last - lcs, 'delins', ins_trimmed, '', None)
#apply_delins


def _intronic_to_genomic(location, transcript):
    """
    Get genomic location from IVS location.

    @arg location: A location.
    @type location: pyparsing.ParseResults
    @arg transcript: todo
    @type transcript: todo

    @return: Genomic location represented by given IVS location.
    @rtype: int
    """
    ivs_number = int(location.IVSNumber)

    if ivs_number < 1 or ivs_number > transcript.CM.numberOfIntrons():
        # Todo: Exception?
        return None

    if location.OffSgn == '+':
        return transcript.CM.getSpliceSite(ivs_number * 2 - 1) + \
               transcript.CM.orientation * int(location.Offset)
    else:
        return transcript.CM.getSpliceSite(ivs_number * 2) - \
               transcript.CM.orientation * int(location.Offset)
#_intronic_to_genomic


def _exonic_to_genomic(location, transcript) :
    """
    Get genomic range from EX location.

    @arg location: A location.
    @type location: pyparsing.ParseResults
    @arg transcript: todo
    @type transcript: todo

    @return: A tuple of:
        - first: Genomic start location represented by given EX location.
        - last:  Genomic end location represented by given EX location.
    @rtype: tuple(int, int)

    @todo: We probably want to treat this as a-?_b+?, so take the centers of
           flanking exons.
    @todo: Exceptions instead of returning None?
    """
    first_exon = int(location.EXNumberStart)
    if first_exon < 1 or first_exon > transcript.CM.numberOfExons():
        return None
    first = transcript.CM.getSpliceSite(first_exon * 2 - 2)

    if location.EXNumberStop:
        last_exon = int(location.EXNumberStop)
        if last_exon < 1 or last_exon > transcript.CM.numberOfExons():
            return None
        last = transcript.CM.getSpliceSite(last_exon * 2 - 1)
    else:
        last = transcript.CM.getSpliceSite(first_exon * 2 - 1)

    return first, last
#_exonic_to_genomic


def _genomic_to_genomic(first_location, last_location):
    """
    Get genomic range from parsed genomic location.

    @arg first_location: The start location (g.) of the variant.
    @type first_location: pyparsing.ParseResults
    @arg last_location: The start location (g.) of the variant.
    @type last_location: pyparsing.ParseResults

    @return: A tuple of:
        - first: Genomic start location represented by given location.
        - last:  Genomic end location represented by given location.
    @rtype: tuple(int, int)

    @todo: Exceptions.
    """
    if not first_location.Main.isdigit():
        # For ? in a position.
        return None, None

    if not last_location.Main.isdigit():
        # For ? in a position.
        return None, None

    first = int(first_location.Main)
    last = int(last_position.Main)

    return first, last


def _coding_to_genomic(first_location, last_location, transcript):
    """
    Get genomic range from parsed c. location.

    @arg first_location: The start location (c.) of the variant.
    @type first_location: pyparsing.ParseResults
    @arg last_location: The start location (c.) of the variant.
    @type last_location: pyparsing.ParseResults
    @arg transcript: todo
    @type transcript: todo

    @return: A tuple of:
        - first: Genomic start location represented by given location.
        - last:  Genomic end location represented by given location.
    @rtype: tuple(int, int)

    @todo: Exceptions.
    """
    if not first_location.Main.isdigit():
        # For ? in a position.
        return None, None

    if not last_location.Main.isdigit():
        # For ? in a position.
        return None, None

    first_main = transcript.CM.main2int(first_location.MainSgn + \
                                        first_location.Main)
    first_offset = _get_offset(first_location)
    first = transcript.CM.x2g(first_main, first_offset)

    last_main = transcript.CM.main2int(last_location.MainSgn + \
                                       last_location.Main)
    last_offset = _get_offset(last_location)
    last = transcript.CM.x2g(last_main, last_offset)

    # Todo: Exceptions.
    if not _check_intronic_position(first_main, first_offset, transcript):
        return None, None
    if not _check_intronic_position(last_main, last_offset, transcript):
        return None, None

    if transcript.CM.orientation == -1:
        first, last = last, first

    return first, last
#_coding_to_genomic


def _raw_variant(mutator, variant, record, transcript, output):
    """
    Process a raw variant.

    @arg mutator: A Mutator object.
    @type mutator: Modules.Mutator.Mutator
    @arg variant: A parsed raw (simple, noncompound) variant.
    @type variant: pyparsing.ParseResults
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg transcript: A transcript object.
    @type transcript: Modules.GenRecord.Locus
    @arg output: The Output object.
    @type output: Modules.Output.Output

    @todo: Documentation.
    @todo: Exceptions.
    """
    if transcript and transcript.CM.orientation == -1:
        s1 = Bio.Seq.reverse_complement(variant.Arg1)
        s2 = Bio.Seq.reverse_complement(variant.Arg2)
    else:
        s1 = variant.Arg1
        s2 = variant.Arg2

    if variant.EXLoc:
        first, last = _exonic_to_genomic(variant.EXLoc, transcript)
        if not first:
            output.addMessage(__file__, 3, 'EPOS', 'Invalid EX position given.')
            return
        if last < first:
            # Todo: huh?
            first, last = last, first
    else:
        if variant.StartLoc:
            if variant.StartLoc.IVSLoc:
                if record.record.molType != 'g':
                    output.addMessage(__file__, 3, 'ENOINTRON', 'Intronic ' \
                        'position given for a non-genomic reference sequence.')
                    return
                first = _intronic_to_genomic(variant.StartLoc.IVSLoc, transcript)
                if not first:
                    output.addMessage(__file__, 3, 'EPOS',
                        'Invalid IVS position given.')
                    return
                last = first
                if variant.EndLoc and variant.EndLoc.IVSLoc:
                    # Todo: fixme
                    last = _intronic_to_genomic(variant.EndLoc.IVSLoc, transcript)
                    if last < first:
                        first, last = last, first
            else:
                if record.record.molType != 'g' and \
                       (_is_intronic(variant.StartLoc) or
                        _is_intronic(variant.EndLoc)):
                    output.addMessage(__file__, 3, 'ENOINTRON', 'Intronic ' \
                        'position given for a non-genomic reference sequence.')
                    return
                first_location = RawVar.StartLoc.PtLoc
                if RawVar.EndLoc:
                    last_location = RawVar.EndLoc.PtLoc
                else:
                    last_location = first_location
                if transcript:
                    first, last = _coding_to_genomic(first_location, last_location, transcript)
                else:
                    first, last = _genomic_to_genomic(first_location, last_location)
                if not first:
                    output.addMessage(__file__, 3, 'ESPLICE', 'Invalid intronic ' \
                        'position given.')
                    return
        else:
            output.addMessage(__file__, 4, 'EUNKNOWN', 'An unknown error occurred.')
            return

    if last < first:
        output.addMessage(__file__, 3, 'ERANGE', 'End position is smaller than ' \
                          'the begin position.')
        return

    if first < 1:
        output.addMessage(__file__, 4, 'ERANGE', 'Position %i is out of range.' %
                          first)
        return

    if last > len(mutator.orig):
        output.addMessage(__file__, 4, 'ERANGE', 'Position %s is out of range.' %
                          last)
        return

    if transcript and _over_splice_site(first, last, transcript.CM.RNA):
        output.addMessage(__file__, 2, 'WOVERSPLICE',
                          'Variant hits one or more splice sites.')

    if variant.MutationType in ['del', 'dup', 'subst', 'delins']:
        __checkOptArg(mutator.orig, first, last, s1, output)

    # Substitution.
    if variant.MutationType == 'subst':
        apply_substitution(first, s1, s2, mutator, record, output)

    # Deletion or duplication.
    if variant.MutationType in ['del', 'dup']:
        apply_deletion_duplication(first, last, variant.MutationType, mutator,
                                   record, output)

    # Inversion.
    if variant.MutationType == 'inv':
        apply_inversion(first, last, mutator, record, output)

    # Insertion.
    if variant.MutationType == 'ins':
        apply_insertion(first, last, s1, mutator, record, output)

    # DelIns.
    if variant.MutationType == 'delins':
        apply_delins(first, last, s1, s2, mutator, record, output)
#_raw_variant


def __ppp(mutator, description, record, output):
    """
    @arg mutator: A Mutator object.
    @type mutator: Modules.Mutator.Mutator
    @arg description: Parsed HGVS variant description.
    @type description: pyparsing.ParseResults
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg output: The Output object.
    @type output: Modules.Output.Output

    @todo: Documentation.
    @todo: Exceptions.
    """
    if not description.RawVar and not description.SingleAlleleVarSet:
        # Nothing to do. Exception?
        return

    if description.RefType == 'r':
        output.addMessage(__file__, 4, "ERNA", "Descriptions on RNA level " \
                          "are not supported.")
        return

    if description.RefType in ['c', 'n']:

        gene, transcript = None, None
        gene_symbol, transcript_id = O.getOutput('geneSymbol')[-1]

        if description.LrgAcc:
            # LRG case, pick the top gene.
            gene = record.record.geneList[0]
            if transcript_id:
                transcript = gene.findLocus(transcript_id)
                if not transcript:
                    output.addMessage(__file__, 4, "ENOTRANSCRIPT",
                        "Multiple transcripts found for gene %s. Please " \
                        "choose from: %s" %(gene.name,
                            ", ".join(gene.listLoci())))
            else:
                # No transcript id given.
                if len(gene.transcriptList) == 1:
                    # No transcript given, only 1 found, pick that.
                    transcript = gene.transcriptList[0]
                else:
                    output.addMessage(__file__, 4, "ENOTRANSCRIPT",
                        "No transcript given for gene %s. Please " \
                        "choose from: %s" %(gene.name,
                            ", ".join(gene.listLoci())))

        else:
            # Not an LRG, find our gene manually.
            genes = record.record.listGenes()
            transcript_id = transcript_id and "%.3i" % int(transcript_id)

            if gene_symbol in genes:
                # We found our gene.
                gene = record.record.findGene(gene_symbol)
            elif (len(genes) == 1) and not(gene_symbol):
                # No gene given and there is only one gene in the record.
                # Todo: message?
                gene = record.record.geneList[0]
            else:
                output.addMessage(__file__, 4, "EINVALIDGENE",
                    "Gene %s not found. Please choose from: %s" % (
                    gene_symbol, ", ".join(genes)))

            if gene:
                # Find transcript.
                transcripts = gene.listLoci()
                if transcript_id in transcripts:
                    # Found our transcript.
                    transcript = gene.findLocus(transcript_id)
                elif (len(transcripts) == 1) and not(transcript_id):
                    # No transcript given and there is only one transcript for
                    # this gene.
                    transcript = gene.transcriptList[0]
                else:
                    output.addMessage(__file__, 4, "ENOTRANSCRIPT",
                        "Multiple transcripts found for gene %s. Please " \
                        "choose from: %s" %(gene.name,
                        ", ".join(gene.listLoci())))

        # Add selected gene symbol to output
        output.addOutput('geneSymbol', (gene and gene.name or '',
                                        transcript and transcript.name or ''))

        # Return if no transcript is selected
        if not transcript:
            # Skip all BatchJobs with the same preColon data.
            output.addOutput('BatchFlags',
                             ('S2', output.getOutput('preColon')[-1]))
            # Explicit return in case of an error.
            return

    else:
        # Not description.RefType in ['c', 'n'].
        transcript = None

    if transcript and not transcript.transcribe:
        return

    if description.SingleAlleleVarSet:
        for var in description.SingleAlleleVarSet:
            _raw_variant(mutator, var.RawVar, record, transcript, output)
    else:
        _raw_variant(mutator, description.RawVar, record, transcript, output)

    if not transcript:
        # Genomic given or error with transcript.
        return

    if not record.record.geneList:
        # EST
        return

    # Add exon table to output.
    for i in range(0, transcript.CM.numberOfExons() * 2, 2):
        acceptor = transcript.CM.getSpliceSite(i)
        donor = transcript.CM.getSpliceSite(i + 1)
        output.addOutput('exonInfo', [acceptor, donor,
                                      transcript.CM.g2c(acceptor),
                                      transcript.CM.g2c(donor)])

    # Add CDS info to output.
    cds_stop = transcript.CM.info()[2]
    output.addOutput('cdsStart_g', transcript.CM.x2g(1, 0))
    output.addOutput('cdsStart_c', 1)
    output.addOutput('cdsStop_g', transcript.CM.x2g(cds_stop, 0))
    output.addOutput('cdsStop_c', cds_stop)

    # Add transcript info to output.
    if transcript.transcribe:
        output.addOutput('myTranscriptDescription', transcript.description)
        output.addOutput('origMRNA',
            str(_splice(mutator.orig, transcript.mRNA.positionList)))
        output.addOutput('mutatedMRNA',
            str(_splice(mutator.mutated,
                        mutator.newSplice(transcript.mRNA.positionList))))

    # Add protein prediction to output.
    if transcript.translate:
        cds_original = Seq(str(_splice(mutator.orig, transcript.CDS.positionList)),
                           IUPAC.unambiguous_dna)
        cds_variant = Seq(str(__nsplice(mutator.mutated,
                                        mutator.newSplice(transcript.mRNA.positionList),
                                        mutator.newSplice(transcript.CDS.location),
                                        transcript.CM.orientation)),
                          IUPAC.unambiguous_dna)

        #output.addOutput('origCDS', cds_original)

        if transcript.CM.orientation == -1:
            cds_original = Bio.Seq.reverse_complement(cds_original)
            cds_variant = Bio.Seq.reverse_complement(cds_variant)

        if '*' in cds_original.translate(table=transcript.txTable)[:-1]:
            output.addMessage(__file__, 3, 'ESTOP',
                              'In frame stop codon found.')
            return

        protein_original = cds_original.translate(table=transcript.txTable,
                                                  to_stop=True)
        protein_variant = cds_variant.translate(table=transcript.txTable,
                                                to_stop=True)

        # Note: addOutput('origCDS', ...) was first before the possible
        #       reverse complement operation above.
        output.addOutput('origCDS', cds_original)
        output.addOutput("newCDS", cds_variant[:(len(str(protein_variant)) + 1) * 3])

        output.addOutput('oldprotein', protein_original + '*')

        if not protein_variant or protein_variant[0] != 'M':
            # Todo: Protein differences are not color-coded,
            # use something like below in _protein_description().
            _print_protein_html(protein_original + '*', 0, 0, output,
                                'oldProteinFancy')
            if str(cds_variant[0:3]) in \
                   Bio.Data.CodonTable.unambiguous_dna_by_id \
                   [transcript.txTable].start_codons:
                output.addOutput('newprotein', '?')
                _print_protein_html('?', 0, 0, output, 'newProteinFancy')
                output.addOutput('altStart', str(cds_variant[0:3]))
                if str(protein_original[1:]) != str(protein_variant[1:]):
                    output.addOutput('altProtein',
                                     'M' + protein_variant[1:] + '*')
                    _print_protein_html('M' + protein_variant[1:] + '*', 0, 0,
                                        output, 'altProteinFancy')
            else :
                output.addOutput('newprotein', '?')
                _print_protein_html('?', 0, 0, O, 'newProteinFancy')

        else:
            cds_length = _cds_length(
                mutator.newSplice(transcript.CDS.positionList))
            descr, first, last_original, last_variant = \
                   _protein_description(cds_length, protein_original,
                                        protein_variant)

            output.addOutput('myProteinDescription', descr)

            _print_protein_html(protein_original + '*', first, last_original,
                                output, 'oldProteinFancy')
            if str(protein_original) != str(protein_variant):
                output.addOutput('newprotein', protein_variant + '*')
                _print_protein_html(protein_variant + '*', first, last_variant,
                                    output, 'newProteinFancy')
#__ppp


def process(description, config, output):
    """
    @return: A GenRecord object.
    @rtype: Modules.GenRecord.GenRecord

    @todo: documentation
    """
    output.addOutput('inputvariant', description)

    parser = Parser.Nomenclatureparser(output)
    parsed_description = parser.parse(description)

    # Todo: remove?
    del parser

    if not parsed_description:
        # Parsing went wrong.
        return None

    if parsed_description.Version:
        record_id = parsed_description.RefSeqAcc + '.' + parsed_description.Version
    else:
        record_id = parsed_description.RefSeqAcc

    gene_symbol = transcript_id = ''

    database = Db.Cache(config.Db)
    if parsed_description.LrgAcc:
        filetype = 'LRG'
        record_id = parsed_description.LrgAcc
        transcript_id = parsed_description.LRGTranscriptID
        retriever = Retriever.LRGRetriever(config.Retriever, output, database)
    else:
        if parsed_description.Gene:
            gene_symbol = parsed_description.Gene.GeneSymbol or ''
            transcript_id = parsed_description.Gene.TransVar or ''
            if parsed_description.Gene.ProtIso:
                output.addMessage(__file__, 4, 'EPROT', 'Indexing by ' \
                                  'protein isoform is not supported.')
        retriever = Retriever.GenBankRetriever(C.Retriever, output, database)
        filetype = 'GB'

    # Add recordType to output for output formatting.
    output.addOutput('recordType', filetype)

    output.addOutput('reference', record_id)

    # Note: geneSymbol[0] is used as a filter for batch runs.
    output.addOutput('geneSymbol', (gene_symbol, transcript_id))

    # Note: preColon is used to filter out Batch entries that will result in
    # identical errors.
    output.addOutput('preColon', description.split(':')[0])
    output.addOutput('variant', description.split(':')[-1])

    retrieved_record = retriever.loadrecord(record_id)

    if not retrieved_record:
        return None

    # Todo: remove?
    del retriever
    del database

    record = GenRecord.GenRecord(output, config.GenRecord)
    record.record = retrieved_record
    record.checkRecord()

    mutator = Mutator.Mutator(record.record.seq, config.Mutator, output)

    # Note: The GenRecord instance is carrying the sequence in .record.seq.
    #       So is the Mutator instance in .mutator.orig.

    __ppp(mutator, parsed_description, record, output)

    # Protein.
    for gene in record.record.geneList:
        for transcript in gene.transcriptList:
            if not ';' in transcript.description \
                   and transcript.CDS and transcript.translate:
                cds_original = Seq(str(_splice(mutator.orig, transcript.CDS.positionList)),
                                   IUPAC.unambiguous_dna)
                cds_variant = Seq(str(__nsplice(mutator.mutated,
                                                mutator.newSplice(transcript.mRNA.positionList),
                                                mutator.newSplice(transcript.CDS.location),
                                                transcript.CM.orientation)),
                                  IUPAC.unambiguous_dna)
                if transcript.CM.orientation == -1:
                    cds_original = Bio.Seq.reverse_complement(cds_original)
                    cds_variant = Bio.Seq.reverse_complement(cds_variant)

                #if '*' in cds_original.translate()[:-1]:
                #    output.addMessage(__file__, 3, "ESTOP",
                #                      "In frame stop codon found.")
                #    return
                ##if

                if not len(cds_original) % 3:
                    try:
                        # FIXME this is a bit of a rancid fix.
                        protein_original = cds_original.translate(table=transcript.txTable,
                                                                  cds=True,
                                                                  to_stop=True)
                    except Bio.Data.CodonTable.TranslationError:
                        output.addMessage(__file__, 4, "ETRANS", "Original " \
                                          "CDS could not be translated.")
                        return record
                    protein_variant = cds_variant.translate(table=transcript.txTable,
                                                            to_stop=True)
                    cds_length = _cds_length(mutator.newSplice(transcript.CDS.positionList))
                    transcript.proteinDescription = _protein_description(cds_length,
                                                                         protein_original,
                                                                         protein_variant)[0]
                else:
                    output.addMessage(__file__, 2, "ECDS", "CDS length is " \
                        "not a multiple of three in gene %s, transcript " \
                        "variant %s." % (gene.name, transcript.name))
                    transcript.proteinDescription = '?'

    reference = output.getOutput('reference')[-1]
    if ';' in record.record.description:
        generated_description = '[' + record.record.description + ']'
    else:
        generated_description = record.record.description

    output.addOutput('genomicDescription', '%s:%c.%s' % \
                     (reference, record.record.molType, generated_description))
    output.addOutput('gDescription', '%c.%s' % \
                     (record.record.molType, generated_description))
    output.addOutput('molType', record.record.molType)

    if record.record.chromOffset:
        if ';' in record.record.chromDescription:
            chromosomal_description = '[' + record.record.chromDescription + ']'
        else:
            chromosomal_description = record.record.chromDescription
        output.addOutput('genomicChromDescription', '%s:%c.%s' % \
                         (record.record.recordId,
                          record.record.molType, chromosomal_description))

    # Now we add variant descriptions for all transcripts, including protein
    # level descriptions. In the same loop, we also create the legend.

    for gene in record.record.geneList:
        for transcript in sorted(gene.transcriptList, key=attrgetter('name')):

            # Note: I don't think genomic_id is ever used, because it is
            # always ''.
            coding_description = ''
            protein_description = ''
            full_description = ''
            full_protein_description = ''
            genomic_id = coding_id = protein_id = ''

            if ';' in transcript.description:
                generated_description = '[' + transcript.description + ']'
            else:
                generated_description = transcript.description

            if record.record._sourcetype == 'LRG':
                if transcript.name:
                    full_description = '%st%s:%c.%s' % \
                                       (reference, transcript.name,
                                        transcript.molType,
                                        generated_description)
                    output.addOutput('descriptions', full_description)
                else:
                    output.addOutput('descriptions', gene.name)
            else:
                full_description = '%s(%s_v%s):%c.%s' % \
                                   (reference, gene.name, transcript.name,
                                    transcript.molType,
                                    generated_description)
                output.addOutput('descriptions', full_description)

            if transcript.molType == 'c':
                coding_description = 'c.%s' % generated_description
                protein_description = transcript.proteinDescription
                if record.record._sourcetype == 'LRG':
                    full_protein_description = '%sp%s:%s' % \
                                               (reference, transcript.name,
                                                protein_description)
                else:
                    full_protein_description = '%s(%s_i%s):%s' % \
                                               (reference, gene.name,
                                                transcript.name,
                                                protein_description)

                coding_id, protein_id = \
                           transcript.transcriptID, transcript.proteinID
                output.addOutput('protDescriptions',
                                 full_protein_description)

            # The 'NewDescriptions' field is currently not used.
            output.addOutput('NewDescriptions',
                             (gene.name, transcript.name,
                              transcript.molType, coding_description,
                              protein_description, genomic_id, coding_id,
                              protein_id, full_description,
                              full_protein_description))

            # Now add to the legend, but exclude nameless transcripts.
            if not transcript.name:
                continue

            output.addOutput('legends',
                             ['%s_v%s' % (gene.name, transcript.name),
                              transcript.transcriptID, transcript.locusTag,
                              transcript.transcriptProduct,
                              transcript.linkMethod])

            if transcript.translate:
                output.addOutput('legends',
                                 ['%s_i%s' % (gene.name, transcript.name),
                                  transcript.proteinID, transcript.locusTag,
                                  transcript.proteinProduct,
                                  transcript.linkMethod])

    # Add GeneSymbol and Transcript Var to the Output object for batch.
    if parsed_description.Gene:
        output.addOutput('geneOfInterest',
                         dict(parsed_description.Gene.items()))
    else:
        output.addOutput('geneOfInterest', dict())

    _add_batch_output(output)

    output.addOutput('original', str(mutator.orig))
    output.addOutput('mutated', str(mutator.mutated))

    # Todo: remove?
    del MUU

    return record
#process


def main(cmd):
    """
    Command line interface to the name checker.

    @todo: documentation
    """
    C = Config.Config()
    O = Output.Output(__file__, C.Output)

    O.addMessage(__file__, -1, "INFO", "Received variant " + cmd)

    RD = process(cmd, C, O)

    O.addMessage(__file__, -1, "INFO", "Finished processing variant " + cmd)

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
        visualisation = O.getOutput("visualisation")
        if visualisation :
            for i in range(len(visualisation)) :
                if i and not i % 3 :
                    print
                print visualisation[i]
            #for
            print
        #if

        reference = O.getOutput("reference")[-1]
        for i in O.getOutput("descriptions") :
            print i
        print
        for i in O.getOutput("protDescriptions") :
            print i
        print

        if RD.record and RD.record._sourcetype == "LRG": #LRG record
            from collections import defaultdict
            toutput = defaultdict(list)
            poutput = defaultdict(list)
            for i in RD.record.geneList:
                for j in i.transcriptList:
                    d = j.description
                    d = ';' in d and '['+d+']' or d
                    if j.name:
                        toutput[i.name].append(
                            "%st%s:%c.%s" % (reference, j.name, j.molType, d))
                    else:
                        pass
                    if j.molType == 'c':
                        poutput[i.name].append(
                                "%sp%s:%s" % (reference, j.name,
                                    j.proteinDescription))
                        poutput[i.name].sort()
                toutput[i.name].sort()

            #Transcript Notation
            print "Following transcripts were affected:"
            for key, values in toutput.items():
                print key
                for value in values:
                    print "\t"+value

            #Protein Notation
            print "\nFollowing proteins were affected:"
            for key, values in poutput.items():
                print key
                for value in values:
                    print "\t"+value
            #for
        #if
        else :
            for i in RD.record.geneList :
                for j in i.transcriptList :
                    if ';' in j.description :
                        print "%s(%s_v%s):%c.[%s]" % (reference, i.name, j.name,
                                                      j.molType, j.description)
                    else :
                        print "%s(%s_v%s):%c.%s" % (reference, i.name, j.name,
                                                    j.molType, j.description)
                        if (j.molType == 'c') :
                            print "%s(%s_i%s):%s" % (reference, i.name, j.name,
                                                     j.proteinDescription)
                    #else
                #for
            #for
        #else

        #Genomic Notation
        rdrd = RD.record.description
        gdescr = ';' in rdrd and '['+rdrd+']' or rdrd
        print "\nGenomic notation:\n\t%s:g.%s" % (reference, gdescr)
        print O.getOutput("genomicChromDescription")

        op = O.getOutput("oldprotein")
        if op :
            print "\nOld protein:"
            #__bprint(op[0], O)
            for i in O.getOutput("oldProteinFancy") :
                print i
            print
        #if
        np = O.getOutput("newprotein")
        if np :
            print "\nNew protein:"
            #__bprint(np[0], O)
            for i in O.getOutput("newProteinFancy") :
                print i
            print
        #if
        ap = O.getOutput("altProtein")
        if ap :
            print "\nAlternative protein using start codon %s:" % \
                O.getOutput("altstart")[0]
            #__bprint(ap[0], O)
            for i in O.getOutput("altProteinFancy") :
                print i
            print
        #if

        for i in O.getOutput("exonInfo") :
            print i
        print
        print O.getOutput("cdsStart")
        print O.getOutput("cdsStop")
        print

        for i in O.getOutput("legends") :
            print i

        print
        print "Restriction sites:"
        for i in O.getOutput("restrictionSites") :
            print i

        print "+++ %s" % O.getOutput("myTranscriptDescription")

    #if
    ### OUTPUT BLOCK ###
    del O
#main


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Please provide a variant'
        sys.exit(1)
    main(sys.argv[1])
