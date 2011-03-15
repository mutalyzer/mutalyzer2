#!/usr/bin/env python

"""
The nomenclature checker.

@todo: Use exceptions for failure handling.
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


def _format_range(pos1, pos2):
    """
    Simplify a range to one position when applicable.

    @arg pos1: First coordinate of a range.
    @type pos1: integer
    @arg pos2: Second coordinate of a range.
    @type pos2: integer

    @return: pos1_pos2 in case of a real range, pos1 otherwise.
    @rtype: string
    """
    if pos1 == pos2:
        return str(pos1)

    return '%i_%i' % (pos1, pos2)
#_format_range


# Used in: checkDeletionDuplication, checkInsertion
def _roll(s, start, stop):
    """
    Determine the variability of a variant by looking at cyclic
    permutations. Not all cyclic permutations are tested at each time, it
    is sufficient to check ``aW'' if ``Wa'' matches (with ``a'' a letter,
    ``W'' a word) when rolling to the left for example.

    @arg s: A reference sequence.
    @type s: string
    @arg start: Start position of the pattern in the reference sequence.
    @type start: int
    @arg stop: End position of the pattern in the reference sequence.
    @type stop: int

    @return: tuple:
        - left  ; Amount of positions that the pattern can be shifted to
                  the left.
        - right ; Amount of positions that the pattern can be shifted to
                  the right.
    @rtype: tuple(int, int)
    """
    pattern = s[start - 1:stop] # Extract the pattern
    pattern_length = len(pattern)

    # Keep rolling to the left as long as a cyclic permutation matches.
    minimum = start - 2
    j = pattern_length - 1
    while minimum > -1 and s[minimum] == pattern[j % pattern_length]:
        j -= 1
        minimum -= 1

    # Keep rolling to the right as long as a cyclic permutation matches.
    maximum = stop
    j = 0
    while maximum < len(s) and s[maximum] == pattern[j % pattern_length]:
        j += 1
        maximum += 1

    return start - minimum - 2, maximum - stop
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


def _over_splice_site(pos1, pos2, splice_sites):
    """
    Check wheter a genomic range (pos1_pos2) hits a splice site.

    @arg pos1: The first coordinate of the range in g. notation.
    @type pos1: int
    @arg pos2: The first coordinate of the range in g. notation.
    @type pos2: int
    @arg sites: A list of splice sites in g. notation.
    @type sites: list(int)

    @return: True if one or more splice sites are hit, False otherwise.
    @rtype: boolean
    """
    sites_iter = iter(splice_sites)
    for acceptor, donor in izip_longest(sites_iter, sites_iter):
        if pos1 < acceptor and pos2 >= acceptor:
            return True
        if donor and pos1 <= donor and pos2 > donor:
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


def _print_protein_html(s, pos1, pos2, O, where):
    """
    Make a fancy representation of a protein and put it in the Output
    object under the name 'where'. The representation contains HTML tags
    and is suitable for viewing in a monospaced font.

    @arg s: A protein sequence.
    @type s: string
    @arg pos1: First position to highlight.
    @type pos1: int
    @arg pos2: Last position to highlight.
    @type pos2: int
    @arg O: The Output object.
    @type O: object
    @arg where: Location in the Output object to store the representation.
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
        output += ' ' + _insert_tag(s[i:i + block], pos1 - i, pos2 - i,
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
        - int    ; Start position of the change.
        - int    ; End position of the change in the first protein.
        - int    ; End position of the change in the second protein.
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
        - int    ; Start position of the change.
        - int    ; End position of the first protein.
        - int    ; End position of the second protein.
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
        - int    ; Start position of the change.
        - int    ; End position of the change in the first protein.
        - int    ; End position of the change in the second protein.
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


# Used in: __rv
def _is_intronic_position(loc):
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
#_intronic_position


# Used in: __normal2g
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


def _get_offset(loc) :
    """
    Convert the offset coordinate in a location (from the Parser) to an
    integer.

    @arg loc: A location.
    @type loc: pyparsing.ParseResults

    @return: Integer representation of the offset coordinate.
    @rtype: int
    """
    if loc.Offset :
        if loc.Offset == '?' : # This is highly debatable.
            return 0
        offset = int(loc.Offset)
        if loc.OffSgn == '-' :
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


def apply_deletion_duplication(start, end, type, mutator, record, O):
    """
    Do a semantic check for a deletion or duplication, do the actual
    deletion/duplication and give it a name.

    @arg start: Genomic start position of the del/dup.
    @type start: int
    @arg end: Genomic end position of the del/dup.
    @type end: int
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
    roll = _roll(mutator.orig, start, end)
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
            if end < acceptor and end + roll[1] >= acceptor:
                shift = acceptor - 1 - end
                break
            if end <= donor and end + roll[1] > donor:
                shift = donor - end
                break

    if shift:
        new_start = start + shift
        new_stop = end + shift
        O.addMessage(__file__, 2, 'WROLL',
            'Sequence "%s" at position %s was given, however, ' \
            'the HGVS notation prescribes that it should be "%s" at ' \
            'position %s.' % (
            mutator.visualiseLargeString(str(mutator.orig[start - 1:end])),
            _format_range(start, end),
            mutator.visualiseLargeString(str(mutator.orig[new_start - 1:new_stop])),
            _format_range(new_start, new_stop)))

    if shift != roll[1]:
        # The original roll was decreased because it crossed a splice site.
        incorrect_start = start + roll[1]
        incorrect_stop = end + roll[1]
        O.addMessage(__file__, 1, 'IROLLBACK',
            'Sequence "%s" at position %s was not corrected to "%s" at ' \
            'position %s, since they reside in different exons.' % (
            mutator.visualiseLargeString(str(mutator.orig[start - 1:end])),
            _format_range(start, end),
            mutator.visualiseLargeString(str(mutator.orig[incorrect_start - 1:incorrect_stop])),
            _format_range(incorrect_start, incorrect_stop)))

    if type == 'del':
        mutator.delM(start, end)
    else :
        mutator.dupM(start, end)

    record.name(start, end, type, '', '', (roll[0], shift))
#apply_deletion_duplication


def apply_inversion(start, end, mutator, record, O) :
    """
    Do a semantic check for an inversion, do the actual inversion, and give
    it a name.

    @arg start: Genomic start position of the inversion.
    @type start: int
    @arg end: Genomic end position of the inversion.
    @type end: int
    @arg mutator: A Mutator object.
    @type mutator: Modules.Mutator.Mutator
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg O: The Output object.
    @type O: Modules.Output.Output

    @todo: Exception instead of O.addMessage().
    """
    snoop = _palinsnoop(mutator.orig[start - 1:end])

    if snoop:
        # We have a reverse-complement-palindromic prefix.
        if snoop == -1 :
            # Actually, not just a prefix, but the entire selected sequence is
            # a 'palindrome'.
            O.addMessage(__file__, 2, 'WNOCHANGE',
                'Sequence "%s" at position %i_%i is a palindrome ' \
                '(its own reverse complement).' % (
                mutator.visualiseLargeString(str(mutator.orig[start - 1:end])),
                start, end))
            return
        else:
            O.addMessage(__file__, 2, 'WNOTMINIMAL',
                'Sequence "%s" at position %i_%i is a partial ' \
                'palindrome (the first %i nucleotide(s) are the reverse ' \
                'complement of the last one(s)), the HGVS notation ' \
                'prescribes that it should be "%s" at position %i_%i.' % (
                mutator.visualiseLargeString(str(mutator.orig[start - 1:end])),
                start, end, snoop,
                mutator.visualiseLargeString(
                    str(mutator.orig[start + snoop - 1: end - snoop])),
                start + snoop, end - snoop))
            start += snoop
            end -= snoop

    mutator.invM(start, end)

    if start == end:
        O.addMessage(__file__, 2, 'WWRONGTYPE', 'Inversion at position ' \
            '%i is actually a substitution.' % start_g)
        record.name(start, start, 'subst', mutator.orig[start - 1],
            Bio.Seq.reverse_complement(mutator.orig[start - 1]), None)
    else :
        record.name(start, end, 'inv', '', '', None)
#apply_inversion


def apply_insertion(start, end, s, mutator, record, O):
    """
    Do a semantic check for an insertion, do the actual insertion, and give
    it a name.

    @arg start: Genomic start position of the insertion.
    @type start: int
    @arg end: Genomic end position of the insertion.
    @type end: int
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
    if start + 1 != end:
        O.addMessage(__file__, 3, 'EINSRANGE',
            '%i and %i are not consecutive positions.' % (start, end))
        return

    if not s or not _is_dna(s):
        O.addMessage(__file__, 3, 'EUNKVAR', 'Although the syntax of this ' \
            'variant is correct, the effect can not be analysed.')
        return

    insertion_length = len(s)

    mutator.insM(start, s)
    new_start = mutator.shiftpos(start)
    new_stop = mutator.shiftpos(start) + insertion_length

    roll = _roll(mutator.mutated, new_start + 1, new_stop)
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
            s, start, start + 1,
            mutator.mutated[new_start + shift:new_stop + shift], start + shift,
            start + shift + insertion_length - 1))
        end += shift - 1
        start = end - insertion_length + 1
        record.name(start, end, 'dup', '', '',
                    (roll[0] + shift - insertion_length, 0))
    else:
        if shift:
            O.addMessage(__file__, 2, 'WROLL', 'Insertion of %s at position ' \
                '%i_%i was given, however, the HGVS notation prescribes ' \
                'that it should be an insertion of %s at position %i_%i.' % (
                s, start, start + 1,
                mutator.mutated[new_start + shift:new_stop + shift],
                new_start + shift, new_start + shift + 1))
        if shift != roll[1]:
            O.addMessage(__file__, 1, 'IROLLBACK',
                'Insertion of %s at position %i_%i was not corrected to an ' \
                'insertion of %s at position %i_%i, since they reside in ' \
                'different exons.' % (
                s, start, start + 1,
                mutator.mutated[new_start + roll[1]:new_stop + roll[1]],
                new_start + roll[1], new_start + roll[1] + 1))
        record.name(start, start + 1, 'ins',
                    mutator.mutated[new_start + shift:new_stop + shift], '',
                    (roll[0], shift))
#apply_insertion


def __ivs2g(location, transcript) :
    """
    @todo: documentation
    """

    ivsNumber = int(location.IVSNumber)

    if ivsNumber < 1 or ivsNumber > transcript.CM.numberOfIntrons() :
        return None

    if location.OffSgn == '+' :
        return transcript.CM.getSpliceSite(ivsNumber * 2 - 1) + \
            transcript.CM.orientation * int(location.Offset)
    return transcript.CM.getSpliceSite(ivsNumber * 2) - \
        transcript.CM.orientation * int(location.Offset)
#__ivs2g

def __ex2g(location, transcript) :
    """
    @todo: documentation
    """

    numberOfExons = transcript.CM.numberOfExons()

    exNumberStart = int(location.EXNumberStart)
    if exNumberStart < 1 or exNumberStart > transcript.CM.numberOfExons() :
        return None
    start_g = transcript.CM.getSpliceSite(exNumberStart * 2 - 2)

    if location.EXNumberStop :
        exNumberStop = int(location.EXNumberStop)
        if exNumberStop < 1 or exNumberStop > transcript.CM.numberOfExons() :
            return None
        stop_g = transcript.CM.getSpliceSite(exNumberStop * 2 - 1)
    else :
        stop_g = transcript.CM.getSpliceSite(exNumberStart * 2 - 1)

    return start_g, stop_g
#__ex2g

def __normal2g(RawVar, transcript) :
    """
    @todo: documentation
    """

    if not RawVar.StartLoc.PtLoc.Main.isdigit() :
        return None, None # For ? in a position.

    start_g = int(RawVar.StartLoc.PtLoc.Main)
    end_g = start_g
    if RawVar.EndLoc :
        if not RawVar.EndLoc.PtLoc.Main.isdigit() : # For ? in a position.
            return None, None
        #end_g = transcript.CM.main2int(
        #    RawVar.EndLoc.PtLoc.MainSgn + RawVar.EndLoc.PtLoc.Main)
        end_g = int(RawVar.EndLoc.PtLoc.Main)
    #if


    # If it is not, convert it to g. notation.
    if transcript :
        start_main = transcript.CM.main2int(RawVar.StartLoc.PtLoc.MainSgn + \
                                            RawVar.StartLoc.PtLoc.Main)
        #if not RawVar.StartLoc.PtLoc.Offset.isdigit() :
        #    return

        start_offset = _get_offset(RawVar.StartLoc.PtLoc)

        if not _check_intronic_position(start_main, start_offset, transcript) :
            return None, None

        start_g = transcript.CM.x2g(start_main, start_offset)
        end_g = start_g
        if RawVar.EndLoc :
            end_main = transcript.CM.main2int(RawVar.EndLoc.PtLoc.MainSgn + \
                                           RawVar.EndLoc.PtLoc.Main)
            #if not RawVar.EndLoc.PtLoc.Offset.isdigit() :
            #    return
            end_offset = _get_offset(RawVar.EndLoc.PtLoc)
            if not _check_intronic_position(end_main, end_offset, transcript) :
                return None, None
            end_g = transcript.CM.x2g(end_main, end_offset)
        #if
        if transcript.CM.orientation == -1 :
            start_g, end_g = end_g, start_g
    #if

    return start_g, end_g
#__normal2g

def __rv(MUU, RawVar, GenRecordInstance, parts, O, transcript) :
    """
    @todo: documentation
    """

    # FIXME check this
    # First assume that the variant is given in g. notation.
    #print RawVar.StartLoc.PtLoc.MainSgn + RawVar.StartLoc.PtLoc.Main
    #print __PtLoc2offset(RawVar.StartLoc.PtLoc)

    Arg1 = RawVar.Arg1
    Arg2 = RawVar.Arg2

    if RawVar.EXLoc :
        start_g, end_g = __ex2g(RawVar.EXLoc, transcript)
        if not start_g :
            O.addMessage(__file__, 3, "EPOS", "Invalid EX position given.")
            return
        #if
        if end_g < start_g : # FIXME
            start_g, end_g = end_g, start_g
    #if
    else :
        if RawVar.StartLoc :
            if RawVar.StartLoc.IVSLoc :
                if GenRecordInstance.record.molType != 'g' :
                    O.addMessage(__file__, 3, "ENOINTRON", "Intronic " \
                        "position given for a non-genomic reference sequence.")
                    return
                start_g = __ivs2g(RawVar.StartLoc.IVSLoc, transcript)
                if not start_g :
                    O.addMessage(__file__, 3, "EPOS",
                        "Invalid IVS position given.")
                    return
                #if
                end_g = start_g
                if RawVar.EndLoc and RawVar.EndLoc.IVSLoc : # FIXME
                    end_g = __ivs2g(RawVar.EndLoc.IVSLoc, transcript)
                    if end_g < start_g :
                        start_g, end_g = end_g, start_g
                #if
            #if
            else :
                if GenRecordInstance.record.molType != 'g' and \
                   (_is_intronic_position(RawVar.StartLoc) or
                    _is_intronic_position(RawVar.EndLoc)) :
                    O.addMessage(__file__, 3, "ENOINTRON", "Intronic " \
                        "position given for a non-genomic reference sequence.")
                    return
                start_g, end_g = __normal2g(RawVar, transcript)
                if not start_g :
                    O.addMessage(__file__, 3, "ESPLICE", "Invalid intronic " \
                        "position given.")
                    return
            #else
        #if
        else :
            O.addMessage(__file__, 4, "EUNKNOWN", "An unknown error occurred.")
            return
        #else
    #else
    if end_g < start_g :
        O.addMessage(__file__, 3, "ERANGE", "End position is smaller than " \
                     "the begin position.")
        return
    #if

    if start_g < 1 :
        O.addMessage(__file__, 4, "ERANGE", "Position %i is out of range." %
                     start_g)
        return
    #if
    if end_g > len(MUU.orig) :
        O.addMessage(__file__, 4, "ERANGE", "Position %s is out of range." %
                     end_g)
        return
    #if

    if transcript and transcript.CM.orientation == -1 :
        Arg1 = Bio.Seq.reverse_complement(RawVar.Arg1)
        Arg2 = Bio.Seq.reverse_complement(RawVar.Arg2)

    if transcript and _over_splice_site(start_g, end_g, transcript.CM.RNA) :
        O.addMessage(__file__, 2, "WOVERSPLICE",
            "Variant hits one or more splice sites.")

    if RawVar.MutationType in ["del", "dup", "subst", "delins"] :
        __checkOptArg(MUU.orig, start_g, end_g, Arg1, O)

    if RawVar.MutationType == "subst" :
        checkSubstitution(start_g, Arg1, Arg2, MUU, GenRecordInstance, O)
    if RawVar.MutationType in ["del", "dup"] :
        checkDeletionDuplication(start_g, end_g, RawVar.MutationType, MUU,
                                 GenRecordInstance, O)
    if RawVar.MutationType == "inv" :
        checkInversion(start_g, end_g, MUU, GenRecordInstance, O)
    if RawVar.MutationType == "ins" :
        checkInsertion(start_g, end_g, Arg1, MUU, GenRecordInstance, O)


    # DelIns.
    if RawVar.MutationType == "delins" :
        if not Arg1 :
            Arg1 = MUU.orig[start_g - 1:end_g]

        if str(Arg1) == str(Arg2) :
            O.addMessage(__file__, 2, "WNOCHANGE",
                "Sequence \"%s\" at position %i_%i is identical to the " \
                "variant." % (
                MUU.visualiseLargeString(str(MUU.orig[start_g - 1:end_g])),
                start_g, end_g))
            return
        #if

        del_part, ins_part, lcp, lcs = _trim_common(Arg1, Arg2)
        if not len(del_part) :
            O.addMessage(__file__, 2, "WWRONGTYPE", "The given DelIns " \
                         "is actually an insertion.")
            checkInsertion(start_g + lcp - 1, start_g + lcp, ins_part, MUU,
                           GenRecordInstance, O)
            return
        #if
        if len(del_part) == 1 and len(ins_part) == 1 :
            O.addMessage(__file__, 2, "WWRONGTYPE", "The given DelIns " \
                         "is actually a substitution.")
            checkSubstitution(start_g + lcp, del_part, ins_part, MUU,
                              GenRecordInstance, O)
            return
        #if
        if not len(ins_part) :
            O.addMessage(__file__, 2, "WWRONGTYPE", "The given DelIns " \
                         "is actually a deletion.")
            checkDeletionDuplication(start_g + lcp, end_g - lcs, "del",
                                     MUU, GenRecordInstance, O)
            return
        #if
        if str(Bio.Seq.reverse_complement(del_part)) == ins_part :
            O.addMessage(__file__, 2, "WWRONGTYPE", "The given DelIns " \
                         "is actually an inversion.")
            checkInversion(start_g + lcp, end_g - lcs, MUU,
                           GenRecordInstance, O)
            return
        #if
        if len(Arg2) != len(ins_part) :
            O.addMessage(__file__, 2, "WNOTMINIMAL",
                "Sequence \"%s\" at position %i_%i has the same prefix or " \
                "suffix as the inserted sequence \"%s\". The HGVS notation " \
                "prescribes that it should be \"%s\" at position %i_%i." % (
                MUU.visualiseLargeString(str(MUU.orig[start_g - 1:end_g])),
                start_g, end_g, Arg2, ins_part, start_g + lcp, end_g - lcs))

        MUU.delinsM(start_g + lcp, end_g - lcs, ins_part)

        GenRecordInstance.name(start_g + lcp, end_g - lcs, "delins", ins_part,
            "", None)
    #if
#__rv

def __ppp(MUU, parts, GenRecordInstance, O) :
    """
    @todo: documentation
    """
    if parts.RawVar or parts.SingleAlleleVarSet :
        if parts.RefType == 'r' :
            O.addMessage(__file__, 4, "ERNA", "Descriptions on RNA level " \
                "are not supported.")
        if parts.RefType in ['c', 'n'] :
            GS, W = None, None
            goi, toi = O.getOutput("geneSymbol")[-1]
            if parts.LrgAcc:                   # LRG
                GS = GenRecordInstance.record.geneList[0] #LRG pick top gene
                if toi:
                    W = GS.findLocus(toi)
                    if not W:
                        O.addMessage(__file__, 4, "ENOTRANSCRIPT",
                            "Multiple transcripts found for gene %s. Please " \
                            "choose from: %s" %(GS.name,
                                ", ".join(GS.listLoci())))
                else:                       # No transcript id given
                    if len(GS.transcriptList) == 1:
                        #No transcript given, only 1 found
                        W = GS.transcriptList[0]
                    else:
                        O.addMessage(__file__, 4, "ENOTRANSCRIPT",
                            "No transcript given for gene %s. Please " \
                            "choose from: %s" %(GS.name,
                                ", ".join(GS.listLoci())))

            #if
            else:
                # gene of interest
                genes = GenRecordInstance.record.listGenes()
                toi = toi and "%.3i" % int(toi)

                if goi in genes: #we found our gene
                    GS = GenRecordInstance.record.findGene(goi)
                elif (len(genes) == 1) and not(goi):
                    #There is only one gene in the Record, message?
                    GS = GenRecordInstance.record.geneList[0]
                else:
                    O.addMessage(__file__, 4, "EINVALIDGENE",
                        "Gene %s not found. Please choose from: %s" % (
                        goi, ", ".join(genes)))

                if GS:
                    #Find Transcript
                    transcripts = GS.listLoci()
                    if toi in transcripts:
                        W = GS.findLocus(toi)
                    elif (len(transcripts) == 1) and not(toi):
                        W = GS.transcriptList[0]
                    else:
                        O.addMessage(__file__, 4, "ENOTRANSCRIPT",
                            "Multiple transcripts found for gene %s. Please " \
                            "choose from: %s" %(GS.name,
                            ", ".join(GS.listLoci())))
            #else

            # Add seletcted geneSymbol to output
            O.addOutput("geneSymbol", (GS and GS.name or "", W and W.name or ""))

            # Return if no transcript is selected
            if not W:
                #Skip all BatchJobs with the same preColon data
                O.addOutput("BatchFlags", ("S2",
                    O.getOutput("preColon")[-1]))
                return None #Explicit return in case of an error
        #if
        else :
            W = None
        #if W and not W.location :
        #    W = None
        if W and not W.transcribe :
            return

        if parts.SingleAlleleVarSet:
            for i in parts.SingleAlleleVarSet :
                __rv(MUU, i.RawVar, GenRecordInstance, parts, O, W)
        else :
            __rv(MUU, parts.RawVar, GenRecordInstance, parts, O, W)


        if not W : # Genomic given or error with transcript
            return
        if not GenRecordInstance.record.geneList : # EST
            return

        for i in range(0, W.CM.numberOfExons() * 2, 2) :
            exonStart = W.CM.getSpliceSite(i)
            exonStop = W.CM.getSpliceSite(i + 1)
            O.addOutput("exonInfo", [exonStart, exonStop,
                W.CM.g2c(exonStart), W.CM.g2c(exonStop)])

        O.addOutput("cdsStart_g", W.CM.x2g(1, 0))
        O.addOutput("cdsStart_c", 1)
        cdsStop = W.CM.info()[2]
        O.addOutput("cdsStop_g", W.CM.x2g(cdsStop, 0))
        O.addOutput("cdsStop_c", cdsStop)

        if W.transcribe :
            O.addOutput("myTranscriptDescription", W.description)

            O.addOutput("origMRNA",
                str(_splice(MUU.orig, W.mRNA.positionList)))
            O.addOutput("mutatedMRNA",
                str(_splice(MUU.mutated, MUU.newSplice(W.mRNA.positionList))))
        #if


        if W.translate :
            cds = Seq(str(_splice(MUU.orig, W.CDS.positionList)),
                      IUPAC.unambiguous_dna)
            cdsm = Seq(str(__nsplice(MUU.mutated,
                                     MUU.newSplice(W.mRNA.positionList),
                                     MUU.newSplice(W.CDS.location),
                                     W.CM.orientation)),
                       IUPAC.unambiguous_dna)
            O.addOutput("origCDS", cds)

            if W.CM.orientation == -1 :
                cds = Bio.Seq.reverse_complement(cds)
                cdsm = Bio.Seq.reverse_complement(cdsm)
            #if

            if '*' in cds.translate(table = W.txTable)[:-1] :
                O.addMessage(__file__, 3, "ESTOP", "In frame stop codon found.")
                return
            #if
            orig = cds.translate(table = W.txTable, to_stop = True)
            O.addOutput("oldprotein", orig + '*')
            trans = cdsm.translate(table = W.txTable, to_stop = True)
            O.addOutput("newCDS", cdsm[:(len(str(trans)) + 1) * 3])

            if not trans or trans[0] != 'M' :
                # Todo: Protein differences are not color-coded,
                # use something like below in _protein_description().
                _print_protein_html(orig + '*', 0, 0, O, "oldProteinFancy")
                if str(cdsm[0:3]) in \
                    Bio.Data.CodonTable.unambiguous_dna_by_id[
                        W.txTable].start_codons :
                    O.addOutput("newprotein", '?')
                    _print_protein_html('?', 0, 0, O, "newProteinFancy")
                    O.addOutput("altStart", str(cdsm[0:3]))
                    if str(orig[1:]) != str(trans[1:]) :
                        O.addOutput("altProtein", 'M' + trans[1:] + '*')
                        _print_protein_html('M' + trans[1:] + '*', 0, 0, O, "altProteinFancy")
                #if
                else :
                    O.addOutput("newprotein", '?')
                    _print_protein_html('?', 0, 0, O, "newProteinFancy")
                #else
            else :
                cdsLen = _cds_length(MUU.newSplice(W.CDS.positionList))
                descr = _protein_description(cdsLen, orig, trans)
                O.addOutput("myProteinDescription", descr[0])

                _print_protein_html(orig + '*', descr[1], descr[2], O,
                    "oldProteinFancy")
                if str(orig) != str(trans) :
                    O.addOutput("newprotein", trans + '*')
                    _print_protein_html(trans + '*', descr[1], descr[3], O,
                        "newProteinFancy")
            #else
        #if
    #if
#__ppp

def process(cmd, C, O) :
    """
    @todo: documentation
    """
    parser = Parser.Nomenclatureparser(O)
    O.addOutput("inputvariant", cmd)
    ParseObj = parser.parse(cmd)
    del parser
    if not ParseObj :
        #Parsing went wrong
        return None     #Excplicit return of None in case of an error

    if ParseObj.Version :
        RetrieveRecord = ParseObj.RefSeqAcc + '.' + ParseObj.Version
    else :
        RetrieveRecord = ParseObj.RefSeqAcc

    D = Db.Cache(C.Db)
    if ParseObj.LrgAcc :
        filetype = "LRG"
        RetrieveRecord = ParseObj.LrgAcc
        geneSymbol = ("", ParseObj.LRGTranscriptID)
        retriever = Retriever.LRGRetriever(C.Retriever, O, D)
    else :
        if ParseObj.Gene:
            geneSymbol = (ParseObj.Gene.GeneSymbol or "",
                    ParseObj.Gene.TransVar or "")
            if ParseObj.Gene.ProtIso :
                O.addMessage(__file__, 4, "EPROT", "Indexing by protein " \
                    "isoform is not supported.")
        else:
            geneSymbol = ("", "")
        retriever = Retriever.GenBankRetriever(C.Retriever, O, D)
        filetype = "GB"

    # Store the recordType for output formatting
    O.addOutput("recordType", filetype)

    # Note concerning objects in outputObject, example:
    # O.getOutput('reference')[-1] countains the last added value
    # O.getOutput('reference')[0] countains the first added value
    # These can refer to the same element
    O.addOutput("reference", RetrieveRecord)

    # The geneSymbol[0] is used as a filter for batch runs
    O.addOutput("geneSymbol", geneSymbol) #tuple(Gene, TransV)

    # preColon is used to filter out Batch entries
    # that will result in identical errors
    O.addOutput("preColon", cmd.split(":")[0])
    O.addOutput("variant", cmd.split(":")[-1])

    record = retriever.loadrecord(RetrieveRecord)
    #if record and record.version and not '.' in RetrieveRecord : #FIXME
    #    O.addOutput("reference", RetrieveRecord + '.' + record.version)
    #else :

    if not record :
        return
    del retriever
    del D

    GenRecordInstance = GenRecord.GenRecord(O, C.GenRecord)
    GenRecordInstance.record = record
    GenRecordInstance.checkRecord()
    #NOTE:  GenRecordInstance is carrying the sequence in   .record.seq
    #       so is the Mutator.Mutator instance MUU          .orig

    MUU = Mutator.Mutator(GenRecordInstance.record.seq, C.Mutator, O)
    __ppp(MUU, ParseObj, GenRecordInstance, O)

    # PROTEIN
    for i in GenRecordInstance.record.geneList :
        #if i.location :
        for j in i.transcriptList :
            if not ';' in j.description and j.CDS and j.translate :
                cds = Seq(str(_splice(MUU.orig, j.CDS.positionList)),
                          IUPAC.unambiguous_dna)
                cdsm = Seq(str(__nsplice(MUU.mutated,
                                         MUU.newSplice(j.mRNA.positionList),
                                         MUU.newSplice(j.CDS.location),
                                         j.CM.orientation)),
                           IUPAC.unambiguous_dna)
                if j.CM.orientation == -1 :
                    cds = Bio.Seq.reverse_complement(cds)
                    cdsm = Bio.Seq.reverse_complement(cdsm)
                #if

                #if '*' in cds.translate()[:-1] :
                #    O.addMessage(__file__, 3, "ESTOP",
                #                 "In frame stop codon found.")
                #    return
                ##if

                if not len(cds) % 3 :
                    try : # FIXME this is a bit of a rancid fix.
                        orig = cds.translate(table = j.txTable, cds = True,
                            to_stop = True)
                    except Bio.Data.CodonTable.TranslationError :
                        O.addMessage(__file__, 4, "ETRANS", "Original " \
                            "CDS could not be translated.")
                        return GenRecordInstance
                    trans = cdsm.translate(table = j.txTable,
                        to_stop = True)

                    cdsLen = _cds_length(MUU.newSplice(j.CDS.positionList))
                    j.proteinDescription = _protein_description(cdsLen, orig,
                        trans)[0]
                #if
                else :
                    O.addMessage(__file__, 2, "ECDS", "CDS length is " \
                        "not a multiple of three in gene %s, transcript " \
                        "variant %s." % (i.name, j.name))
                    j.proteinDescription = '?'
    # /PROTEIN

    reference = O.getOutput("reference")[-1]
    if ';' in GenRecordInstance.record.description :
        descr = '['+GenRecordInstance.record.description+']'
    else:
        descr = GenRecordInstance.record.description

    O.addOutput("genomicDescription", "%s:%c.%s" % (reference,
        GenRecordInstance.record.molType, descr))
    O.addOutput("gDescription", "%c.%s" % (
        GenRecordInstance.record.molType, descr))
    O.addOutput("molType", GenRecordInstance.record.molType)

    if GenRecordInstance.record.chromOffset :
        if ';' in GenRecordInstance.record.chromDescription :
            chromDescr = '['+GenRecordInstance.record.chromDescription+']'
        else:
            chromDescr = GenRecordInstance.record.chromDescription

        O.addOutput("genomicChromDescription", "%s:%c.%s" % (
            GenRecordInstance.record.recordId,
            GenRecordInstance.record.molType, chromDescr))
    #if

    if GenRecordInstance.record._sourcetype == "LRG": #LRG record
        for i in GenRecordInstance.record.geneList:
            for j in sorted(i.transcriptList, key = attrgetter("name")) :
                (iName, jName, mType, cDescr, pDescr,
                        gAcc, cAcc, pAcc, fullDescr, fullpDescr) =\
                    (i.name, j.name, j.molType, "", "", "", "", "", "", "")

                if ';' in j.description:
                    descr = '['+j.description+']'
                else:
                    descr = j.description

                if j.name:
                    fullDescr =\
                        "%st%s:%c.%s" % (reference, j.name, j.molType, descr)
                    O.addOutput("descriptions", fullDescr)
                #if
                else:
                    O.addOutput("descriptions", (i.name))

                if j.molType == 'c':
                    cDescr = "c.%s" % descr
                    pDescr = j.proteinDescription
                    fullpDescr = "%sp%s:%s" % (reference, j.name, pDescr)
                    O.addOutput("protDescriptions", fullpDescr)
                    cAcc, pAcc = j.transcriptID, j.proteinID
                #if

                O.addOutput("NewDescriptions", (
                    iName, jName, mType, cDescr, pDescr, gAcc,
                    cAcc, pAcc, fullDescr, fullpDescr))
            #for
        #for
    #if
    else :
        for i in GenRecordInstance.record.geneList :
            for j in sorted(i.transcriptList, key = attrgetter("name")) :
                (iName, jName, mType, cDescr, pDescr,
                        gAcc, cAcc, pAcc, fullDescr, fullpDescr) =\
                    (i.name, j.name, j.molType, "", "", "", "", "", "", "")

                if ';' in j.description :
                    descr = '['+j.description+']'
                else:
                    descr = j.description

                fullDescr = "%s(%s_v%s):%c.%s" % (reference,\
                        iName, jName, mType, descr)
                O.addOutput("descriptions", fullDescr)

                if (j.molType == 'c') :
                    cDescr = "c.%s" % descr
                    pDescr = j.proteinDescription
                    fullpDescr = "%s(%s_i%s):%s" % (
                        reference, iName, jName, pDescr)
                    O.addOutput("protDescriptions", fullpDescr)
                    cAcc, pAcc = j.transcriptID, j.proteinID
                #if

                O.addOutput("NewDescriptions", (
                    iName, jName, mType, cDescr, pDescr, gAcc,
                    cAcc, pAcc, fullDescr, fullpDescr))
            #for
        #for
    #else


    # LEGEND
    for i in GenRecordInstance.record.geneList :
        for j in sorted(i.transcriptList, key = attrgetter("name")) :

            if not j.name: continue #Exclude nameless transcripts

            O.addOutput("legends", ["%s_v%s" % (i.name, j.name),
                        j.transcriptID, j.locusTag,
                        j.transcriptProduct, j.linkMethod])
            if j.translate :
                O.addOutput("legends", ["%s_i%s" % (i.name, j.name),
                    j.proteinID, j.locusTag,
                    j.proteinProduct, j.linkMethod])
        #for

    #Add GeneSymbol and Transcript Var to the Output object for batch
    if ParseObj.Gene:
        O.addOutput("geneOfInterest", dict(ParseObj.Gene.items()))
    else:
        O.addOutput("geneOfInterest", dict())

    _add_batch_output(O)

    O.addOutput("original", str(MUU.orig))
    O.addOutput("mutated", str(MUU.mutated))
    del MUU

    return GenRecordInstance
    #if
#process

def main(cmd) :
    """
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

if __name__ == "__main__" :
    if len(sys.argv) > 1:
        main(sys.argv[1])
#if
