#!/usr/bin/python

"""
Prototype of a module that can generate a HGVS description of the variant(s)
leading from one sequence to an other.

@requires: Bio.Seq
"""

import collections
from Bio import Seq
from Bio.SeqUtils import seq3
from Bio.Data import CodonTable

from mutalyzer.util import longest_common_prefix, longest_common_suffix
from mutalyzer.util import palinsnoop, roll
from mutalyzer import models


# Maximum size of the LCS matrix
MAX_MATRIX_SIZE = 8000000


class LCS(object):
    """
    Class that calculates a Longest Common Substring matrix once and provides
    a function to find the LCS of a submatrix.
    """

    __delim = ' '

    def __init__(self, s1, s2, lcp, s1_end, s2_end, DNA=False):
        """
        Initialise the class.

        @arg s1: A string.
        @type s1: str
        @arg s2: A string.
        @type s2: str
        @arg lcp: The length of the longest common prefix of {s1} and {s2}.
        @type lcp: int
        @arg s1_end: End of the substring in {s1}.
        @type s1_end:
        @arg s2_end: End of the substring in {s2}.
        @type s2_end: int
        @arg DNA:
        @type DNA: bool
        """
        self.__lcp = lcp
        self.__s1 = s1[self.__lcp:s1_end]
        self.__s2 = s2[self.__lcp:s2_end]
        self.__s2_len = s2_end - lcp
        self.__matrix = self.LCSMatrix(self.__s1, self.__s2)

        self.__s2_rc = None
        self.__matrix_rc = None
        if DNA:
            self.__s2_rc = Seq.reverse_complement(s2[self.__lcp:s2_end])
            self.__matrix_rc = self.LCSMatrix(self.__s1, self.__s2_rc)
        #if
    #__init__

    def __str__(self):
        """
        Return a graphical representation of the LCS matrix, mainly for
        debugging.

        @returns: A graphical representation of the LCS matrix.
        @rtype: str
        """
        return self.visMatrix((0, len(self.__s1)), (0, len(self.__s2)))
    #__str__

    def visMatrix(self, r1, r2, rc=False):
        """
        Return a graphical representation of the LCS matrix, mainly for
        debugging.

        @returns: A graphical representation of the LCS matrix.
        @rtype: str
        """
        nr1 = r1[0] - self.__lcp, r1[1] - self.__lcp
        nr2 = r2[0] - self.__lcp, r2[1] - self.__lcp

        M = self.__matrix
        s2 = self.__s2
        if rc:
            M = self.__matrix_rc
            s2 = self.__s2_rc

        out = self.__delim.join(self.__delim + '-' + s2[nr2[0]:nr2[1]]) + '\n'
        for i in range(nr1[0], nr1[1] + 1):
            out += (('-' + self.__s1)[i] + self.__delim +
                self.__delim.join(map(lambda x: str(M[i][x]),
                range(nr2[0], nr2[1] + 1))) + '\n')

        return out
    #visMatrix

    def LCSMatrix(self, s1, s2):
        """
        Calculate the Longest Common Substring matrix.

        @arg s1: A string.
        @type s1: str
        @arg s2: A string.
        @type s2: str

        @returns: A matrix with the LCS of {s1}[i], {s2}[j] at position i, j.
        @rval: list[list[int]]
        """
        y_max = len(s1) + 1
        x_max = len(s2) + 1
        M = [[0] * x_max for i in range(y_max)]

        for x in range(1, y_max):
            for y in range(1, x_max):
                if s1[x - 1] == s2[y - 1]:
                    M[x][y] = M[x - 1][y - 1] + 1

        return M
    #LCSMatrix

    def findMax(self, r1, r2, rc=False):
        """
        Find the LCS in a submatrix of {M}, given by the ranges {r1} and {r2}.

        @arg r1: A range for the first dimension of M.
        @type r1: tuple(int, int)
        @arg r2: A range for the second dimension of M.
        @type r2: tuple(int, int)
        @arg rc: Use the reverse complement matrix.
        @type rc: bool

        @returns:
        @rtype: tuple(int, int, int)
        """
        longest, x_longest, y_longest = 0, 0, 0
        nr1 = r1[0] - self.__lcp, r1[1] - self.__lcp
        nr2 = r2[0] - self.__lcp, r2[1] - self.__lcp

        M = self.__matrix
        if rc:
            M = self.__matrix_rc
            nr2 = self.__s2_len - nr2[1], self.__s2_len - nr2[0]
        #if

        for i in range(nr1[0], nr1[1] + 1):
            x_relative = i - nr1[0]

            for j in range(nr2[0], nr2[1] + 1):
                y_relative = j - nr2[0]
                realVal = min(M[i][j], x_relative, y_relative)

                if realVal > longest:
                    longest, x_longest, y_longest = (realVal, x_relative,
                        y_relative)
            #for
        #for

        return x_longest, y_longest, longest
    #findMax
#LCS

def makeFSTables(table_id):
    """
    For every pair of amino acids, calculate the set of possible amino acids in
    a different reading frame. Do this for both alternative reading frames (+1
    and +2).

    @arg table_id: Coding table ID.
    @type table_id: int
    @returns: Two dictionaries for the two alternative reading frames.
    @rtype: tuple(dict, dict)
    """
    # Make the forward translation table.
    table = dict(CodonTable.unambiguous_dna_by_id[table_id].forward_table)
    for i in CodonTable.unambiguous_dna_by_id[table_id].stop_codons:
        table[i] = '*'

    # Make the reverse translation table.
    reverse_table = collections.defaultdict(list)
    for i in table:
        reverse_table[table[i]].append(i)

    # Make the frame shift tables.
    FS1 = collections.defaultdict(set)
    FS2 = collections.defaultdict(set)
    for AA_i in reverse_table:
        for AA_j in reverse_table:
            for codon_i in reverse_table[AA_i]:
                for codon_j in reverse_table[AA_j]:
                    FS1[AA_i + AA_j].add(table[(codon_i + codon_j)[1:4]]) # +1.
                    FS2[AA_i + AA_j].add(table[(codon_i + codon_j)[2:5]]) # +2.
                #for
    return FS1, FS2
#makeFSTables

def __makeOverlaps(peptide):
    """
    Make a list of overlapping 2-mers of {peptide} in order of appearance.

    @arg peptide: A peptide sequence.
    @type peptide: str
    @returns: All 2-mers of {peptide} in order of appearance.
    @rtype: list(str)
    """
    return map(lambda x: peptide[x:x+2], range(len(peptide) - 1))
#__makeOverlaps

def __options(pList, peptidePrefix, FS, output):
    """
    Enumerate all peptides that could result from a frame shift.

    @arg pList: List of overlapping 2-mers of a peptide.
    @type pList: list(str)
    @arg peptidePrefix: Prefix of a peptide in the alternative reading frame.
    @type peptidePrefix: str
    @arg FS: Frame shift table.
    @type FS: dict
    @arg output: List of peptides, should be empty initially.
    @type output: list(str)
    """
    if not pList:
        output.append(peptidePrefix)
        return
    #if
    for i in FS[pList[0]]:
        __options(pList[1:], peptidePrefix + i, FS, output)
#__options

def enumFS(peptide, FS):
    """
    Enumerate all peptides that could result from a frame shift.

    @arg peptide: Original peptide sequence.
    @type peptide: str
    @arg FS: Frame shift table.
    @type FS: dict
    """
    output = []

    __options(__makeOverlaps(peptide), "", FS, output)
    return output
#enumFS

def fitFS(peptide, altPeptide, FS):
    """
    Check whether peptide {altPeptide} is a possible frame shift of peptide
    {peptide}.

    @arg peptide: Original peptide sequence.
    @type peptide: str
    @arg altPeptide: Observed peptide sequence.
    @type altPeptide: str
    @arg FS: Frame shift table.
    @type FS: dict
    """
    if len(peptide) < len(altPeptide):
        return False

    pList = __makeOverlaps(peptide)

    for i in range(len(altPeptide)):
        if not altPeptide[i] in FS[pList[i]]:
            return False
    return True
#fitFS

class RawVar(models.RawVar):
    """
    Container for a raw variant.

    To use this class correctly, do not supply more than the minimum amount of
    data. The {description()} function may not work properly if too much
    information is given.

    Example: if {end} is initialised for a substitution, a range will be
      retuned, resulting in a description like: 100_100A>T
    """

    def __init__(self, DNA=True, start=0, start_offset=0, end=0, end_offset=0,
        type="none", deleted="", inserted="", shift=0, startAA="", endAA="",
        term=0):
        """
        Initialise the class with the appropriate values.

        @arg start: Start position.
        @type start: int
        @arg start_offset:
        @type start_offset: int
        @arg end: End position.
        @type end: int
        @arg end_offset:
        @type end_offset: int
        @arg type: Variant type.
        @type type: str
        @arg deleted: Deleted part of the reference sequence.
        @type deleted: str
        @arg inserted: Inserted part.
        @type inserted: str
        @arg shift: Amount of freedom.
        @type shift: int
        """
        # TODO: Will this container be used for all variants, or only genomic?
        #       start_offset and end_offset may be never used.
        self.DNA = DNA
        self.start = start
        self.start_offset = start_offset
        self.end = end
        self.end_offset = end_offset
        self.type = type
        self.deleted = deleted
        self.inserted = inserted
        self.shift = shift
        self.startAA = startAA
        self.endAA = endAA
        self.term = term
        self.hgvs = self.description()
        self.hgvsLength = self.descriptionLength()
    #__init__

    def __DNADescription(self):
        """
        Give the HGVS description of the raw variant stored in this class.

        Note that this function relies on the absence of values to make the
        correct description. Also see the comment in the class definition.

        @returns: The HGVS description of the raw variant stored in this class.
        @rtype: str
        """
        if not self.start:
            return "="

        descr = "%i" % self.start

        if self.end:
            descr += "_%i" % self.end

        if self.type != "subst":
            descr += "%s" % self.type

            if self.inserted:
                return descr + "%s" % self.inserted
            return descr
        #if

        return descr + "%s>%s" % (self.deleted, self.inserted)
    #__DNADescription

    def __proteinDescription(self):
        """
        Give the HGVS description of the raw variant stored in this class.

        Note that this function relies on the absence of values to make the
        correct description. Also see the comment in the class definition.

        @returns: The HGVS description of the raw variant stored in this class.
        @rtype: str
        """
        if self.type == "unknown":
            return "?"
        if not self.start:
            return "="

        descr = ""
        if not self.deleted:
            if self.type == "ext":
                descr += '*'
            else:
                descr += "%s" % seq3(self.startAA)
        #if
        else:
            descr += "%s" % seq3(self.deleted)
        descr += "%i" % self.start
        if self.end:
            descr += "_%s%i" % (seq3(self.endAA), self.end)
        if self.type not in ["subst", "stop", "ext", "fs"]:
            descr += self.type
        if self.inserted:
            descr += "%s" % seq3(self.inserted)

        if self.type == "stop":
            return descr + '*'
        if self.term:
            return descr + "%s*%i" % (self.type, self.term)
        return descr
    #__proteinDescription

    def __DNADescriptionLength(self):
        """
        Give the standardised length of the HGVS description of the raw variant
        stored in this class.

        Note that this function relies on the absence of values to make the
        correct description. Also see the comment in the class definition.

        @returns: The standardised length of the HGVS description of the raw
            variant stored in this class.
        @rtype: int
        """
        if not self.start : # `=' or `?'
            return 1

        descrLen = 1 # Start position.

        if self.end : # '_' and end position.
            descrLen += 2

        if self.type != "subst":
            descrLen += len(self.type)

            if self.inserted:
                return descrLen + len(self.inserted)
            return descrLen
        #if

        return 4 # Start position, '>' and end position.
    #__DNAdescriptionLength

    def __proteinDescriptionLength(self):
        """
        Give the standardised length of the HGVS description of the raw variant
        stored in this class.

        Note that this function relies on the absence of values to make the
        correct description. Also see the comment in the class definition.

        @returns: The standardised length of the HGVS description of the raw
            variant stored in this class.
        @rtype: int
        """
        if not self.start: # =
            return 1

        descrLen = 1      # Start position.
        if not self.deleted and self.type == "ext":
            descrLen += 1 # *
        else:
            descrLen += 3 # One amino acid.
        if self.end:
            descrLen += 5 # `_' + one amino acid + end position.
        if self.type not in ["subst", "stop", "ext", "fs"]:
            descrLen += len(self.type)
        if self.inserted:
            descrLen += 3 * len(self.inserted)
        if self.type == "stop":
            return descrLen + 1 # *
        if self.term:
            return descrLen + len(self.type) + 2 # `*' + length until stop.
        return descrLen
    #__proteinDescriptionLength

    def description(self):
        """
        """
        if self.DNA:
            return self.__DNADescription()
        return self.__proteinDescription()
    #description

    def descriptionLength(self):
        """
        Give the standardised length of the HGVS description of the raw variant
        stored in this class.

        @returns: The standardised length of the HGVS description of the raw
            variant stored in this class.
        @rtype: int
        """
        if self.DNA:
            return self.__DNADescriptionLength()
        return self.__proteinDescriptionLength()
    #descriptionLength
#RawVar

def alleleDescription(allele):
    """
    Convert a list of raw variants to an HGVS allele description.

    @arg allele: A list of raw variants representing an allele description.
    @type allele: list(RawVar)

    @returns: The HGVS description of {allele}.
    @rval: str
    """
    if len(allele) > 1:
        return "[%s]" % ';'.join(map(lambda x : x.hgvs, allele))
    return allele[0].hgvs
#alleleDescription

def alleleDescriptionLength(allele):
    """
    Calculate the standardised length of an HGVS allele description.

    @arg allele: A list of raw variants representing an allele description.
    @type allele: list(RawVar)

    @returns: The standardised length of the HGVS description of {allele}.
    @rval: int
    """
    # NOTE: Do we need to count the ; and [] ?
    return sum(map(lambda x : x.hgvsLength, allele))
#alleleDescriptionLength

def printpos(s, start, end, fill = 0):
    """
    For debugging purposes.
    """
    # TODO: See if this can partially replace or be merged with the
    #       visualisation in the __mutate() function of mutator.py
    fs = 10 # Flank size.

    return "%s %s%s %s" % (s[start - fs:start], s[start:end], '-' * fill,
        s[end:end + fs])
#printpos

def DNA_description(M, s1, s2, s1_start, s1_end, s2_start, s2_end):
    """
    Give an allele description of the change from {s1} to {s2} in the range
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

    @returns: A list of RawVar objects, representing the allele.
    @rval: list(RawVar)
    """
    # TODO: Instead of copying this function and adjusting it to make it work
    #       for proteins, consider disabling parts like the inversion.
    # TODO: Think about frameshift descriptions.

    # Nothing happened.
    if s1 == s2:
        return [RawVar()]

    # Insertion / Duplication.
    if s1_start == s1_end:
        ins_length = s2_end - s2_start
        shift5, shift3 = roll(s2, s2_start + 1, s2_end)
        shift = shift5 + shift3

        s1_start += shift3
        s1_end += shift3
        s2_start += shift3
        s2_end += shift3

        if s2_start - ins_length >= 0 and \
            s1[s1_start - ins_length:s1_start] == s2[s2_start:s2_end]:

            if ins_length == 1:
                return [RawVar(start=s1_start, type="dup", shift=shift)]
            return [RawVar(start=s1_start - ins_length + 1, end=s1_end,
                type="dup", shift=shift)]
        #if
        return [RawVar(start=s1_start, end=s1_start + 1,
            inserted=s2[s2_start:s2_end], type="ins", shift=shift)]
    #if

    # Deletion.
    if s2_start == s2_end:
        shift5, shift3 = roll(s1, s1_start + 1, s1_end)
        shift = shift5 + shift3

        s1_start += shift3 + 1
        s1_end += shift3

        if s1_start == s1_end:
            return [RawVar(start=s1_start, type="del", shift=shift)]
        return [RawVar(start=s1_start, end=s1_end, type="del", shift=shift)]
    #if

    # Substitution.
    if s1_start + 1 == s1_end and s2_start + 1 == s2_end:
        return [RawVar(start=s1_start + 1, deleted=s1[s1_start],
            inserted=s2[s2_start], type="subst")]

    # Simple InDel.
    if s1_start + 1 == s1_end:
        return [RawVar(start=s1_start + 1, inserted=s2[s2_start:s2_end],
            type="delins")]

    # TODO: Refactor the code after this point.

    # At this stage, we either have an inversion, an indel or a Compound
    # variant.
    s1_end_f, s2_end_f, lcs_f_len = M.findMax((s1_start, s1_end),
        (s2_start, s2_end))
    s1_end_r, s2_end_r, lcs_r_len = M.findMax((s1_start, s1_end),
        (s2_start, s2_end), rc=True)

    # Palindrome snooping.
    trim = palinsnoop(s1[s1_start + s1_end_r - lcs_r_len:s1_start + s1_end_r])
    if trim == -1 :   # Full palindrome.
        lcs_r_len = 0 # s1_end_r and s2_end_r should not be used after this.

    # Inversion or Compound variant.
    default = [RawVar(start=s1_start + 1, end=s1_end,
        inserted=s2[s2_start:s2_end], type="delins")]

    if not (lcs_f_len or lcs_r_len) : # Optimisation, not really needed.
        return default

    # Inversion.
    if lcs_f_len <= lcs_r_len:
        if trim > 0 : # Partial palindrome.
            s1_end_r -= trim
            s2_end_r -= trim
            lcs_r_len -= 2 * trim
        #if

        # Simple Inversion.
        if s2_end - s2_start == lcs_r_len and s1_end - s1_start == lcs_r_len:
            return [RawVar(start=s1_start + 1, end=s1_end, type="inv")]

        r1_len = s1_end_r - lcs_r_len
        r2_len = s1_end - s1_start - s1_end_r
        m1_len = s2_end_r - lcs_r_len
        m2_len = s2_end - s2_start - s2_end_r

        # The flanks of the inversion (but not both) can be empty, so we
        # generate descriptions conditionally.
        leftRv = []
        rightRv = []
        if r1_len or m2_len:
            lcs = len(longest_common_suffix(s1[s1_start:s1_start + r1_len],
                s2[s2_start:s2_start + m2_len]))
            leftRv = DNA_description(M, s1, s2,
                s1_start, s1_start + r1_len - lcs,
                s2_start, s2_start + m2_len - lcs)
        #if
        if r2_len or m1_len:
            lcp = len(longest_common_prefix(s1[s1_end - r2_len:s1_end],
                s2[s2_end - m1_len:s2_end]))
            rightRv = DNA_description(M, s1, s2,
                s1_end - r2_len + lcp, s1_end, s2_end - m1_len + lcp, s2_end)
        #if

        partial = leftRv + [RawVar(start=s1_start + r1_len + 1,
            end=s1_end - r2_len, type="inv")] + rightRv
    #if

    # Compound variant.
    else:
        r1_len = s1_end_f - lcs_f_len
        r2_len = s1_end - s1_start - s1_end_f
        m1_len = s2_end_f - lcs_f_len
        m2_len = s2_end - s2_start - s2_end_f

        partial = DNA_description(M, s1, s2, s1_start, s1_start + r1_len,
            s2_start, s2_start + m1_len) + DNA_description(M, s1, s2,
            s1_end - r2_len, s1_end, s2_end - m2_len, s2_end)
    #else

    if alleleDescriptionLength(partial) <= alleleDescriptionLength(default):
        return partial
    return default
#DNA_description

def protein_description(M, s1, s2, s1_start, s1_end, s2_start, s2_end):
    """
    Give an allele description of the change from {s1} to {s2} in the range
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

    @returns: A list of RawVar objects, representing the allele.
    @rval: list(RawVar)
    """
    if s1 == '?' or s2 == '?':
        return [RawVar(DNA=False, type="unknown")]

    # One of the sequences is missing.
    if not (s1 and s2):
        return [RawVar(DNA=False)]

    # Nothing happened.
    if s1 == s2:
        return [RawVar(DNA=False)]

    # Substitution.
    if s1_start + 1 == s1_end and  s2_start + 1 == s2_end:
        return [RawVar(DNA=False, start=s1_start + 1, deleted=s1[s1_start],
            inserted=s2[s2_start], type="subst")]

    # Insertion / Duplication / Extention.
    if s1_start == s1_end:
        if len(s1) == s1_start + 1:
            return [RawVar(DNA=False, start=s1_start + 1,
                inserted=s2[s2_start], term=abs(len(s1) - len(s2)),
                type="ext")]
        ins_length = s2_end - s2_start
        shift5, shift3 = roll(s2, s2_start + 1, s2_end)
        shift = shift5 + shift3

        s1_start += shift3
        s1_end += shift3
        s2_start += shift3
        s2_end += shift3

        if s2_start - ins_length >= 0 and \
            s1[s1_start - ins_length:s1_start] == s2[s2_start:s2_end]:

            if ins_length == 1:
                return [RawVar(DNA=False, start=s1_start,
                    startAA=s1[s1_start - 1], type="dup", shift=shift)]
            return [RawVar(DNA=False, start=s1_start - ins_length + 1,
                startAA=s1[s1_start - ins_length], end=s1_end,
                endAA=s1[s1_end - 1], type="dup", shift=shift)]
        #if
        return [RawVar(DNA=False, start=s1_start, startAA=s1[s1_start - 1],
            end=s1_start + 1, endAA=s1[s1_end], inserted=s2[s2_start:s2_end],
            type="ins", shift=shift)]
    #if

    # Deletion / Inframe stop.
    if s2_start == s2_end:
        if len(s2) == s2_end + 1:
            return [RawVar(DNA=False, start=s1_start + 1, startAA=s1[s1_start],
                type="stop")]
        shift5, shift3 = roll(s1, s1_start + 1, s1_end)
        shift = shift5 + shift3

        s1_start += shift3 + 1
        s1_end += shift3

        if s1_start == s1_end:
            return [RawVar(DNA=False, start=s1_start, startAA=s1[s1_start - 1],
                 type="del", shift=shift)]
        return [RawVar(DNA=False, start=s1_start, startAA=s1[s1_start - 1],
            end=s1_end, endAA=s1[s1_end - 1], type="del", shift=shift)]
    #if

    # Simple InDel.
    if s1_start + 1 == s1_end:
        return [RawVar(DNA=False, start=s1_start + 1, startAA=s1[s1_start],
            inserted=s2[s2_start:s2_end], type="delins")]

    # Frameshift.
    if len(s1) == s1_end + 1 and len(s2) == s2_end + 1:
        # TODO: the FS tables should be made for all coding tables in advance.
        FS1, FS2 = makeFSTables(1) # Standard coding table.

        if (fitFS(s1[s1_start + 1:], s2[s2_start + 1:], FS1) or
            fitFS(s1[s1_start + 1:], s2[s2_start + 1:], FS2) or
            fitFS(s2[s2_start + 1:], s1[s1_start + 2:], FS1) or
            fitFS(s2[s2_start + 1:], s1[s1_start + 2:], FS2)): 
            return [RawVar(DNA=False, start=s1_start + 1, deleted=s1[s1_start],
                inserted=s2[s2_start], term=len(s2) - s2_start, type="fs")]
    #if

    # At this point, we have an indel.
    s1_end_f, s2_end_f, lcs_f_len = M.findMax((s1_start, s1_end),
        (s2_start, s2_end))

    # InDel.
    default = [RawVar(DNA=False, start=s1_start + 1, startAA=s1[s1_start],
        end=s1_end, endAA=s1[s1_end], inserted=s2[s2_start:s2_end],
        type="delins")]

    if not lcs_f_len : # Optimisation, not really needed.
        return default

    r1_len = s1_end_f - lcs_f_len
    r2_len = s1_end - s1_start - s1_end_f
    m1_len = s2_end_f - lcs_f_len
    m2_len = s2_end - s2_start - s2_end_f

    partial = protein_description(M, s1, s2, s1_start, s1_start + r1_len,
        s2_start, s2_start + m1_len) + protein_description(M, s1, s2,
        s1_end - r2_len, s1_end, s2_end - m2_len, s2_end)

    if alleleDescriptionLength(partial) <= alleleDescriptionLength(default):
        return partial
    return default
#protein_description

def describe(original, mutated, DNA=True):
    """
    Convenience function for DNA_description().

    @arg original:
    @type original: str
    @arg mutated:
    @type mutated: str

    @returns: A list of RawVar objects, representing the allele.
    @rval: list(RawVar)
    """
    s1 = str(original)
    s2 = str(mutated)
    lcp = len(longest_common_prefix(s1, s2))
    lcs = len(longest_common_suffix(s1[lcp:], s2[lcp:]))
    s1_end = len(s1) - lcs
    s2_end = len(s2) - lcs

    if (s1_end - lcp) * (s2_end - lcp) > MAX_MATRIX_SIZE:
        return

    if not DNA:
        M = LCS(s1, s2, lcp, s1_end, s2_end)
        return protein_description(M, s1, s2, lcp, s1_end, lcp, s2_end)
    #if
    M = LCS(s1, s2, lcp, s1_end, s2_end, DNA=True)
    return DNA_description(M, s1, s2, lcp, s1_end, lcp, s2_end)
#describe
