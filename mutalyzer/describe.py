"""
Prototype of a module that can generate a HGVS description of the variant(s)
leading from one sequence to an other.

@requires: Bio.Seq
"""


from __future__ import unicode_literals

import collections
from Bio.SeqUtils import seq3
from Bio.Data import CodonTable

from mutalyzer.util import longest_common_prefix, longest_common_suffix
from mutalyzer.util import palinsnoop, roll, reverse_complement
from mutalyzer import models

from extractor import extractor


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
    @type peptide: unicode
    @returns: All 2-mers of {peptide} in order of appearance.
    @rtype: list(unicode)
    """
    return map(lambda x: peptide[x:x+2], range(len(peptide) - 1))
#__makeOverlaps

def __options(pList, peptidePrefix, FS, output):
    """
    Enumerate all peptides that could result from a frame shift.

    @arg pList: List of overlapping 2-mers of a peptide.
    @type pList: list(unicode)
    @arg peptidePrefix: Prefix of a peptide in the alternative reading frame.
    @type peptidePrefix: unicode
    @arg FS: Frame shift table.
    @type FS: dict
    @arg output: List of peptides, should be empty initially.
    @type output: list(unicode)
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
    @type peptide: unicode
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
    @type peptide: unicode
    @arg altPeptide: Observed peptide sequence.
    @type altPeptide: unicode
    @arg FS: Frame shift table.
    @type FS: dict
    """
    # Todo: This is a temporary fix to prevent crashing on frameshift
    #     detection (I think bug #124).
    return False

    if len(peptide) < len(altPeptide):
        return False

    pList = __makeOverlaps(peptide)

    for i in range(len(altPeptide)):
        if not altPeptide[i] in FS[pList[i]]:
            return False
    return True
#fitFS

def findFS(peptide, altPeptide, FS):
    """
    Find the longest part of {altPeptide} that fits in {peptide} in a certain
    frame given by {FS}.

    @arg peptide: Original peptide sequence.
    @type peptide: unicode
    @arg altPeptide: Observed peptide sequence.
    @type altPeptide: unicode
    @arg FS: Frame shift table.
    @type FS: dict

    @returns: The length and the offset in {peptide} of the largest frameshift.
    @rtype: tuple(int, int)
    """
    pList = __makeOverlaps(peptide)
    maxFS = 0
    fsStart = 0

    for i in range(len(pList))[::-1]:
        for j in range(min(i + 1, len(altPeptide))):
            if not altPeptide[::-1][j] in FS[pList[i - j]]:
                break
        if j >= maxFS:
            maxFS = j
            fsStart = i - j + 2
        #if
    #for

    return maxFS - 1, fsStart
#findFS

class Seq(object):
    """
    Container for an inserted sequence.
    """
    #def __init__(self, sequence="", reference="", start=0, end=0,
    def __init__(self, sequence="", start=0, end=0, reverse=False):
        """
        """
        self.sequence = sequence
        #self.reference = reference
        self.start = start
        self.end = end
        self.reverse = reverse

        self.type = "trans"
        if self.sequence:
            self.type = "ins"
    #__init__

    def __str__(self):
        if self.sequence:
            return self.sequence

        inverted = "inv" if self.reverse else ""
        return "{}_{}{}".format(self.start, self.end, inverted)
    #__str__

    def __nonzero__(self):
         return bool(self.sequence)

    #def dump(self):
    #    """
    #    Debug function.
    #    """
    #    if self.sequence:
    #        return self.sequence
    #    return self.reference[self.start:self.end]
    ##dump
#Seq

class SeqList(list):
    def __str__(self):
        representation = ';'.join(map(str, self))

        if len(self) > 1:
            return "[{}]".format(representation)
        return representation
    #__str__

    def __nonzero__(self):
        return bool(self[0])
#SeqList

class RawVar(models.RawVar):
    """
    Container for a raw variant.
    """

    def __init__(self, DNA=True, start=0, start_offset=0, end=0, end_offset=0,
            sample_start=0, sample_start_offset=0, sample_end=0,
            sample_end_offset=0, type="none", deleted=SeqList([Seq()]),
            inserted=SeqList([Seq()]), shift=0, startAA="", endAA="", term=0):
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
        @arg sample_start: Start position.
        @type sample_start: int
        @arg sample_start_offset:
        @type sample_start_offset: int
        @arg sample_end: End position.
        @type sample_end: int
        @arg sample_end_offset:
        @type sample_end_offset: int
        @arg type: Variant type.
        @type type: unicode
        @arg deleted: Deleted part of the reference sequence.
        @type deleted: unicode
        @arg inserted: Inserted part.
        @type inserted: object
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
        self.sample_start = sample_start
        self.sample_start_offset = sample_start_offset
        self.sample_end = sample_end
        self.sample_end_offset = sample_end_offset
        self.type = type
        self.deleted = deleted
        self.inserted = inserted
        self.shift = shift
        self.startAA = startAA
        self.endAA = endAA
        self.term = term
        self.update()
        #self.hgvs = self.description()
        #self.hgvsLength = self.descriptionLength()
    #__init__

    def __DNADescription(self):
        """
        Give the HGVS description of the raw variant stored in this class.

        @returns: The HGVS description of the raw variant stored in this class.
        @rtype: unicode
        """
        if not self.start:
            return "="

        descr = "{}".format(self.start)

        if self.start != self.end:
            descr += "_{}".format(self.end)

        if self.type != "subst":
            descr += "{}".format(self.type)

            if self.type in ("ins", "delins"):
                return descr + "{}".format(self.inserted)
            return descr
        #if

        return descr + "{}>{}".format(self.deleted, self.inserted)
    #__DNADescription

    def __proteinDescription(self):
        """
        Give the HGVS description of the raw variant stored in this class.

        Note that this function relies on the absence of values to make the
        correct description. Also see the comment in the class definition.

        @returns: The HGVS description of the raw variant stored in this class.
        @rtype: unicode
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
                descr += "{}".format(seq3(self.startAA))
        #if
        else:
            descr += "{}".format(seq3(self.deleted))
        descr += "{}".format(self.start)
        if self.end:
            descr += "_{}{}".format(seq3(self.endAA), self.end)
        if self.type not in ["subst", "stop", "ext", "fs"]: # fs is not a type
            descr += self.type
        if self.inserted:
            descr += "{}".format(seq3(self.inserted))

        if self.type == "stop":
            return descr + '*'
        if self.term:
            return descr + "fs*{}".format(self.term)
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
        if not self.start: # `=' or `?'
            return 1

        descrLen = 1 # Start position.

        if self.end: # '_' and end position.
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

    def update(self):
        """
        """
        self.hgvs = self.description()
        self.hgvsLength = self.descriptionLength()
    #update

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
#DescribeRawVar

def allele_description(allele):
    """
    Convert a list of raw variants to an HGVS allele description.

    @arg allele: A list of raw variants representing an allele description.
    @type allele: list(DescribeRawVar)

    @returns: The HGVS description of {allele}.
    @rval: unicode
    """
    if len(allele) > 1:
        return "[{}]".format(';'.join(map(lambda x: x.hgvs, allele)))
    return allele[0].hgvs
#allele_description

def alleleDescriptionLength(allele):
    """
    Calculate the standardised length of an HGVS allele description.

    @arg allele: A list of raw variants representing an allele description.
    @type allele: list(DescribeRawVar)

    @returns: The standardised length of the HGVS description of {allele}.
    @rval: int
    """
    # NOTE: Do we need to count the ; and [] ?
    return sum(map(lambda x: x.hgvsLength, allele))
#alleleDescriptionLength

def printpos(s, start, end, fill=0):
    """
    For debugging purposes.
    """
    # TODO: See if this can partially replace or be merged with the
    #       visualisation in the __mutate() function of mutator.py
    fs = 10 # Flank size.

    return "{} {}{} {}".format(s[start - fs:start], s[start:end], '-' * fill,
        s[end:end + fs])
#printpos

def var2RawVar(s1, s2, var, seq_list=[], DNA=True):
    """
    """
    # Unknown.
    if s1 == '?' or s2 == '?':
        return [RawVar(DNA=DNA, type="unknown")]

    # Insertion / Duplication.
    if var.reference_start == var.reference_end:
        ins_length = var.sample_end - var.sample_start
        shift5, shift3 = roll(s2, var.sample_start + 1, var.sample_end)
        shift = shift5 + shift3

        var.reference_start += shift3
        var.reference_end += shift3
        var.sample_start += shift3
        var.sample_end += shift3

        if (var.sample_start - ins_length >= 0 and
            s1[var.reference_start - ins_length:var.reference_start] ==
            s2[var.sample_start:var.sample_end]):

            return RawVar(DNA=DNA, start=var.reference_start - ins_length + 1,
                end=var.reference_end, type="dup", shift=shift,
                sample_start=var.sample_start + 1, sample_end=var.sample_end,
                inserted=SeqList([Seq(sequence=s2[
                var.sample_start:var.sample_end])]))
        #if
        return RawVar(DNA=DNA, start=var.reference_start,
            end=var.reference_start + 1,
            inserted=seq_list or
            SeqList([Seq(sequence=s2[var.sample_start:var.sample_end])]),
            type="ins", shift=shift, sample_start=var.sample_start + 1,
            sample_end=var.sample_end)
    #if

    # Deletion.
    if var.sample_start == var.sample_end:
        shift5, shift3 = roll(s1, var.reference_start + 1, var.reference_end)
        shift = shift5 + shift3

        var.reference_start += shift3
        var.reference_end += shift3

        return RawVar(DNA=DNA, start=var.reference_start + 1,
            end=var.reference_end, type="del", shift=shift,
            sample_start=var.sample_start, sample_end=var.sample_end + 1,
            deleted=SeqList([Seq(sequence=s1[
                var.reference_start:var.reference_end])]))
    #if

    # Substitution.
    if (var.reference_start + 1 == var.reference_end and
        var.sample_start + 1 == var.sample_end):

        return RawVar(DNA=DNA, start=var.reference_start + 1,
            end=var.reference_end, sample_start=var.sample_start + 1,
            sample_end=var.sample_end, type="subst",
            deleted=SeqList([Seq(sequence=s1[var.reference_start])]),
            inserted=SeqList([Seq(sequence=s2[var.sample_start])]))
    #if

    # Inversion.
    if var.type & extractor.REVERSE_COMPLEMENT:
        trim = palinsnoop(s1[var.reference_start:var.reference_end])

        if trim > 0: # Partial palindrome.
            var.reference_end -= trim
            var.sample_end -= trim
        #if

        return RawVar(DNA=DNA, start=var.reference_start + 1,
            end=var.reference_end, type="inv",
            sample_start=var.sample_start + 1, sample_end=var.sample_end,
            deleted=SeqList([Seq(sequence=s1[
                var.reference_start:var.reference_end])]),
            inserted=SeqList([Seq(sequence=s2[
                var.sample_start:var.reference_end])]))
    #if

    # InDel.
    return RawVar(DNA=DNA, start=var.reference_start + 1,
        end=var.reference_end, deleted=SeqList([Seq(sequence=s1[
                var.reference_start:var.reference_end])]), inserted=seq_list or
        SeqList([Seq(sequence=s2[var.sample_start:var.sample_end])]),
        type="delins", sample_start=var.sample_start + 1,
        sample_end=var.sample_end)
#var2RawVar

def describe(s1, s2, DNA=True):
    """
    Give an allele description of the change from {s1} to {s2}.

    arg s1: Sequence 1.
    type s1: unicode
    arg s2: Sequence 2.
    type s2: unicode

    @returns: A list of DescribeRawVar objects, representing the allele.
    @rval: list(DescribeRawVar)
    """
    description = []

    if not DNA:
        FS1, FS2 = makeFSTables(1)
        longestFSf = max(findFS(s1, s2, FS1), findFS(s1, s2, FS2))
        longestFSr = max(findFS(s2, s1, FS1), findFS(s2, s1, FS2))

        if longestFSf > longestFSr:
            print s1[:longestFSf[1]], s1[longestFSf[1]:]
            print s2[:len(s2) - longestFSf[0]], s2[len(s2) - longestFSf[0]:]
            s1_part = s1[:longestFSf[1]]
            s2_part = s2[:len(s2) - longestFSf[0]]
            term = longestFSf[0]
        #if
        else:
            print s1[:len(s1) - longestFSr[0]], s1[len(s1) - longestFSr[0]:]
            print s2[:longestFSr[1]], s2[longestFSr[1]:]
            s1_part = s1[:len(s1) - longestFSr[0]]
            s2_part = s2[:longestFSr[1]]
            term = len(s2) - longestFSr[1]
        #else

        s1_part = s1
        s2_part = s2
        for variant in extractor.extract(unicode(s1_part), len(s1_part),
            unicode(s2_part), len(s2_part), 1):
            description.append(var2RawVar(s1, s2, variant, DNA=DNA))

        if description:
            description[-1].term = term + 2
            description[-1].update()
        #if
    #if
    else: # DNA description extraction, the only thing that `works'.
        in_transposition = 0

        for variant in extractor.extract(unicode(s1), len(s1), unicode(s2), len(s2),
                0):
            #print variant.type, variant.reference_start, variant.reference_end, variant.sample_start, variant.sample_end
            #print variant.type & extractor.TRANSPOSITION_OPEN, variant.type & extractor.TRANSPOSITION_CLOSE

            if variant.type & extractor.TRANSPOSITION_OPEN:
                if not in_transposition:
                    seq_list = SeqList()
                in_transposition += 1
            #if

            if in_transposition:
                if variant.type & extractor.IDENTITY:
                    seq_list.append(Seq(#reference=s1,
                        start=variant.sample_start + 1, end=variant.sample_end,
                        reverse=False))
                elif variant.type & extractor.REVERSE_COMPLEMENT:
                    seq_list.append(Seq(#reference=s1,
                        start=variant.sample_start + 1, end=variant.sample_end,
                        reverse=True))
                else:
                    seq_list.append(Seq(
                        sequence=s2[variant.sample_start:variant.sample_end]))
            #if
            elif not (variant.type & extractor.IDENTITY):
               description.append(var2RawVar(s1, s2, variant, DNA=DNA))

            if variant.type & extractor.TRANSPOSITION_CLOSE:
                in_transposition -= 1

                if not in_transposition:
                    description.append(var2RawVar(s1, s2, variant, seq_list,
                        DNA=DNA))
                #for i in seq_list:
                #    print i.dump()
            #if
        #for
    #else

    # Nothing happened.
    if not description:
        return [RawVar(DNA=DNA)]

    return description
#describe
