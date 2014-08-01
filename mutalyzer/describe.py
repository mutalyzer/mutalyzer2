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


def printpos(s, start, end, fill=0):
    """
    For debugging purposes.
    """
    # TODO: See if this can partially replace or be merged with the
    #       visualisation in the __mutate() function of mutator.py
    fs = 10 # Flank size.

    return "{} {}{} {}".format(s[start - fs:start], s[start:end], '-' * fill,
        s[end:end + fs])


def make_fs_tables(table_id):
    """
    For every pair of amino acids, calculate the set of possible amino acids in
    a different reading frame. Do this for both alternative reading frames (+1
    and +2).

    :arg table_id: Coding table ID.
    :type table_id: int

    :returns: Two dictionaries for the two alternative reading frames.
    :rtype: tuple(dict, dict)
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
    fs1 = collections.defaultdict(set)
    fs2 = collections.defaultdict(set)
    for aa_i in reverse_table:
        for aa_j in reverse_table:
            for codon_i in reverse_table[aa_i]:
                for codon_j in reverse_table[aa_j]:
                    fs1[aa_i + aa_j].add(table[(codon_i + codon_j)[1:4]]) # +1.
                    fs2[aa_i + aa_j].add(table[(codon_i + codon_j)[2:5]]) # +2.
                #for
    return fs1, fs2
#make_fs_tables

def _peptide_overlaps(peptide):
    """
    Make a list of overlapping 2-mers of {peptide} in order of appearance.

    :arg peptide: A peptide sequence.
    :type peptide: unicode

    :returns: All 2-mers of {peptide} in order of appearance.
    :rtype: list(unicode)
    """
    return map(lambda x: peptide[x:x+2], range(len(peptide) - 1))
#_peptide_overlaps

def _options(peptide_overlaps, peptide_prefix, fs, output):
    """
    Enumerate all peptides that could result from a frame shift.

    :arg peptide_overlaps: List of overlapping 2-mers of a peptide.
    :type peptide_overlaps: list(unicode)
    :arg peptide_prefix: Prefix of a peptide in the alternative reading frame.
    :type peptide_prefix: unicode
    :arg fs: Frame shift table.
    :type fs: dict
    :arg output: List of peptides, should be empty initially.
    :type output: list(unicode)
    """
    if not peptide_overlaps:
        output.append(peptide_prefix)
        return
    #if
    for i in fs[peptide_overlaps[0]]:
        _options(peptide_overlaps[1:], peptide_prefix + i, fs, output)
#_options

def enum_fs(peptide, fs):
    """
    Enumerate all peptides that could result from a frame shift.

    :arg peptide: Original peptide sequence.
    :type peptide: unicode
    :arg fs: Frame shift table.
    :type fs: dict

    :returns: List of peptides that could result from a frame shift.
    :rtype: list(unicode)
    """
    output = []

    _options(_peptide_overlaps(peptide), "", fs, output)
    return output
#enum_fs

def fit_fs(peptide, alternative_peptide, fs):
    """
    Check whether peptide {alternative_peptide} is a possible frame shift of
    peptide {peptide}.

    :arg peptide: Original peptide sequence.
    :type peptide: unicode
    :arg alternative_peptide: Observed peptide sequence.
    :type alternative_peptide: unicode
    :arg fs: Frame shift table.
    :type fs: dict

    :returns: Whether {alternative_peptide} is a frameshift of {peptide}.
    :rtype: bool
    """
    # Todo: This is a temporary fix to prevent crashing on frameshift
    #     detection (I think bug #124).
    return False

    if len(peptide) < len(alternative_peptide):
        return False

    peptide_overlaps = _peptide_overlaps(peptide)

    for i in range(len(alternative_peptide)):
        if not alternative_peptide[i] in fs[peptide_overlaps[i]]:
            return False
    return True
#fit_fs

def find_fs(peptide, alternative_peptide, fs):
    """
    Find the longest part of {alternative_peptide} that fits in {peptide} in a
    certain frame given by {fs}.

    :arg peptide: Original peptide sequence.
    :type peptide: unicode
    :arg alternative_peptide: Observed peptide sequence.
    :type alternative_peptide: unicode
    :arg fs: Frame shift table.
    :type fs: dict

    :returns: The length and the offset in {peptide} of the largest frameshift.
    :rtype: tuple(int, int)
    """
    peptide_overlaps = _peptide_overlaps(peptide)
    max_fs = 0
    fs_start = 0

    for i in range(len(peptide_overlaps))[::-1]:
        for j in range(min(i + 1, len(alternative_peptide))):
            if not alternative_peptide[::-1][j] in fs[peptide_overlaps[i - j]]:
                break
        if j >= max_fs:
            max_fs = j
            fs_start = i - j + 2
        #if
    #for

    return max_fs - 1, fs_start
#find_fs

class Seq(object):
    """
    Container for an inserted sequence.
    """
    def __init__(self, sequence="", start=0, end=0, reverse=False):
        """
        """
        self.sequence = sequence
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
#Seq

class SeqList(list):
    def __str__(self):
        representation = ';'.join(map(str, self))

        if len(self) > 1:
            return "[{}]".format(representation)
        return representation
    #__str__
#SeqList

class HGVSVar(object):
    def update(self):
        self.hgvs = str(self)
        self.hgvs_length = len(self)
    #update
#HGVSVar

class DNAVar(models.DNAVar, HGVSVar):
    """
    Container for a DNA variant.
    """

    def __init__(self, start=0, start_offset=0, end=0, end_offset=0,
            sample_start=0, sample_start_offset=0, sample_end=0,
            sample_end_offset=0, type="none", deleted=SeqList([Seq()]),
            inserted=SeqList([Seq()]), shift=0):
        """
        Initialise the class with the appropriate values.

        :arg start: Start position.
        :type start: int
        :arg start_offset:
        :type start_offset: int
        :arg end: End position.
        :type end: int
        :arg end_offset:
        :type end_offset: int
        :arg sample_start: Start position.
        :type sample_start: int
        :arg sample_start_offset:
        :type sample_start_offset: int
        :arg sample_end: End position.
        :type sample_end: int
        :arg sample_end_offset:
        :type sample_end_offset: int
        :arg type: Variant type.
        :type type: unicode
        :arg deleted: Deleted part of the reference sequence.
        :type deleted: unicode
        :arg inserted: Inserted part.
        :type inserted: object
        :arg shift: Amount of freedom.
        :type shift: int
        """
        # TODO: Will this container be used for all variants, or only genomic?
        #       start_offset and end_offset may be never used.
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
        self.update()
    #__init__

    def __str__(self):
        """
        Give the HGVS description of the raw variant stored in this class.

        :returns: The HGVS description of the raw variant stored in this class.
        :rtype: unicode
        """
        if self.type == "unknown":
            return "?"
        if self.type == "none":
            return "="

        description = "{}".format(self.start)

        if self.start != self.end:
            description += "_{}".format(self.end)

        if self.type != "subst":
            description += "{}".format(self.type)

            if self.type in ("ins", "delins"):
                return description + "{}".format(self.inserted)
            return description
        #if

        return description + "{}>{}".format(self.deleted, self.inserted)
    #__str__

    def __len__(self):
        """
        Give the standardised length of the HGVS description of the raw variant
        stored in this class.

        Note that this function relies on the absence of values to make the
        correct description. Also see the comment in the class definition.

        :returns: The standardised length of the HGVS description of the raw
            variant stored in this class.
        :rtype: int
        """
        if self.type in ("none", "unknown"): # `=' or `?'
            return 1

        description_length = 1 # Start position.

        if self.start != self.end: # '_' and end position.
            description_length += 2

        if self.type != "subst":
            description_length += len(self.type)

            if self.type in ("ins", "delins"):
                return description_length + len(self.inserted)
            return description_length
        #if

        return 4 # Start position, '>' and end position.
    #__len__
#DNAVar

class ProteinVar(models.ProteinVar, HGVSVar):
    """
    Container for a raw variant.
    """

    def __init__(self, start=0, end=0, sample_start=0, sample_end=0,
            type="none", deleted=SeqList([Seq()]), inserted=SeqList([Seq()]),
            shift=0, term=0):
        """
        Initialise the class with the appropriate values.

        :arg start: Start position.
        :type start: int
        :arg end: End position.
        :type end: int
        :arg sample_start: Start position.
        :type sample_start: int
        :arg sample_end: End position.
        :type sample_end: int
        :arg type: Variant type.
        :type type: str
        :arg deleted: Deleted part of the reference sequence.
        :type deleted: str
        :arg inserted: Inserted part.
        :type inserted: object
        :arg shift: Amount of freedom.
        :type shift: int
        :arg term:
        :type term:
        """
        self.start = start
        self.end = end
        self.sample_start = sample_start
        self.sample_end = sample_end
        self.type = type
        self.deleted = deleted
        self.inserted = inserted
        self.shift = shift
        self.term = term
        self.update()
    #__init__

    def __str__(self):
        """
        Give the HGVS description of the raw variant stored in this class.

        Note that this function relies on the absence of values to make the
        correct description. Also see the comment in the class definition.

        :returns: The HGVS description of the raw variant stored in this class.
        :rtype: unicode
        """
        if self.type == "unknown":
            return "?"
        if self.type == "none":
            return "="

        description = ""
        if not self.deleted:
            if self.type == "ext":
                description += '*'
            else:
                description += "{}".format(seq3(self.start_aa))
        #if
        else:
            description += "{}".format(seq3(self.deleted))
        description += "{}".format(self.start)
        if self.end:
            description += "_{}{}".format(seq3(self.end_aa), self.end)
        if self.type not in ["subst", "stop", "ext", "fs"]: # fs is not a type
            description += self.type
        if self.inserted:
            description += "{}".format(seq3(self.inserted))

        if self.type == "stop":
            return description + '*'
        if self.term:
            return description + "fs*{}".format(self.term)
        return description
    #__str__

    def __len__(self):
        """
        Give the standardised length of the HGVS description of the raw variant
        stored in this class.

        Note that this function relies on the absence of values to make the
        correct description. Also see the comment in the class definition.

        :returns: The standardised length of the HGVS description of the raw
            variant stored in this class.
        :rtype: int
        """
        if not self.start: # =
            return 1

        description_length = 1      # Start position.
        if not self.deleted and self.type == "ext":
            description_length += 1 # *
        else:
            description_length += 3 # One amino acid.
        if self.end:
            description_length += 5 # `_' + one amino acid + end position.
        if self.type not in ["subst", "stop", "ext", "fs"]:
            description_length += len(self.type)
        if self.inserted:
            description_length += 3 * len(self.inserted)
        if self.type == "stop":
            return description_length + 1 # *
        if self.term:
            return description_length + len(self.type) + 2 # `*' + until stop.
        return description_length
    #__len__
#ProteinVar

class Allele(list):
    def __str__(self):
        """
        Convert a list of raw variants to an HGVS allele description.

        :returns: The HGVS description of {allele}.
        :rtype: str
        """
        if len(self) > 1:
            return "[{}]".format(';'.join(map(lambda x: x.hgvs, self)))
        return self[0].hgvs
    #__str__

    def length(self):
        """
        Calculate the standardised length of an HGVS allele description.

        :returns: The standardised length of the HGVS description of {allele}.
        :rtype: int
        """
        # NOTE: Do we need to count the ; and [] ?
        return sum(map(lambda x: x.hgvs_length, self))
    #length
#Allele


def var_to_rawvar(s1, s2, var, seq_list=[], container=DNAVar):
    """
    """
    # Unknown.
    if s1 == '?' or s2 == '?':
        return [container(type="unknown")]

    # Insertion / Duplication.
    if var.reference_start == var.reference_end:
        ins_length = var.sample_end - var.sample_start
        shift5, shift3 = roll(s2, var.sample_start + 1, var.sample_end)
        shift = shift5 + shift3

        var.reference_start += shift3
        var.reference_end += shift3
        var.sample_start += shift3
        var.sample_end += shift3

        if not seq_list and (var.sample_start - ins_length >= 0 and
            s1[var.reference_start - ins_length:var.reference_start] ==
            s2[var.sample_start:var.sample_end]):

            return container(start=var.reference_start - ins_length + 1,
                end=var.reference_end, type="dup", shift=shift,
                sample_start=var.sample_start + 1, sample_end=var.sample_end,
                inserted=SeqList([Seq(sequence=s2[
                var.sample_start:var.sample_end])]))
        #if
        return container(start=var.reference_start,
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

        return container(start=var.reference_start + 1,
            end=var.reference_end, type="del", shift=shift,
            sample_start=var.sample_start, sample_end=var.sample_end + 1,
            deleted=SeqList([Seq(sequence=s1[
                var.reference_start:var.reference_end])]))
    #if

    # Substitution.
    if (var.reference_start + 1 == var.reference_end and
        var.sample_start + 1 == var.sample_end):

        return container(start=var.reference_start + 1,
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

        return container(start=var.reference_start + 1,
            end=var.reference_end, type="inv",
            sample_start=var.sample_start + 1, sample_end=var.sample_end,
            deleted=SeqList([Seq(sequence=s1[
                var.reference_start:var.reference_end])]),
            inserted=SeqList([Seq(sequence=s2[
                var.sample_start:var.reference_end])]))
    #if

    # InDel.
    return container(start=var.reference_start + 1,
        end=var.reference_end, deleted=SeqList([Seq(sequence=s1[
                var.reference_start:var.reference_end])]), inserted=seq_list or
        SeqList([Seq(sequence=s2[var.sample_start:var.sample_end])]),
        type="delins", sample_start=var.sample_start + 1,
        sample_end=var.sample_end)
#var_to_rawvar

def describe(s1, s2, dna=True):
    """
    Give an allele description of the change from {s1} to {s2}.

    :arg s1: Sequence 1.
    :type s1: unicode
    :arg s2: Sequence 2.
    :type s2: unicode

    :returns: A list of RawVar objects, representing the allele.
    :rtype: list(RawVar)
    """
    description = Allele()

    if not dna:
        fs1, fs2 = make_fs_tables(1)
        longest_fs_f = max(find_fs(s1, s2, fs1), find_fs(s1, s2, fs2))
        longest_fs_r = max(find_fs(s2, s1, fs1), find_fs(s2, s1, fs2))

        if longest_fs_f > longest_fs_r:
            print s1[:longest_fs_f[1]], s1[longest_fs_f[1]:]
            print s2[:len(s2) - longest_fs_f[0]], \
                s2[len(s2) - longest_fs_f[0]:]
            s1_part = s1[:longest_fs_f[1]]
            s2_part = s2[:len(s2) - longest_fs_f[0]]
            term = longest_fs_f[0]
        #if
        else:
            print s1[:len(s1) - longest_fs_r[0]], \
                s1[len(s1) - longest_fs_r[0]:]
            print s2[:longest_fs_r[1]], s2[longest_fs_r[1]:]
            s1_part = s1[:len(s1) - longest_fs_r[0]]
            s2_part = s2[:longest_fs_r[1]]
            term = len(s2) - longest_fs_r[1]
        #else

        s1_part = s1
        s2_part = s2
        for variant in extractor.extract(unicode(s1_part), len(s1_part),
            unicode(s2_part), len(s2_part), 1):
            description.append(var_to_rawvar(s1, s2, variant, container=ProteinVar))

        if description:
            description[-1].term = term + 2
            description[-1].update()
        #if
    #if
    else: # DNA description extraction, the only thing that `works'.
        in_transposition = 0

        for variant in extractor.extract(unicode(s1), len(s1), unicode(s2), len(s2),
                0):
            print variant.type, variant.reference_start, variant.reference_end, variant.sample_start, variant.sample_end, variant.transposition_start, variant.transposition_end
            print variant.type & extractor.TRANSPOSITION_OPEN, variant.type & extractor.TRANSPOSITION_CLOSE

            if variant.type & extractor.TRANSPOSITION_OPEN:
                if not in_transposition:
                    seq_list = SeqList()
                in_transposition += 1
            #if

            if in_transposition:
                if variant.type & extractor.IDENTITY:
                    seq_list.append(Seq(#reference=s1,
                        start=variant.transposition_start + 1,
                        end=variant.transposition_end, reverse=False))
                elif variant.type & extractor.REVERSE_COMPLEMENT:
                    seq_list.append(Seq(#reference=s1,
                        start=variant.transposition_start + 1,
                        end=variant.transposition_end, reverse=True))
                else:
                    seq_list.append(Seq(
                        sequence=s2[variant.sample_start:variant.sample_end]))
            #if
            elif not (variant.type & extractor.IDENTITY):
                description.append(var_to_rawvar(s1, s2, variant))

            if variant.type & extractor.TRANSPOSITION_CLOSE:
                in_transposition -= 1

                if not in_transposition:
                    description.append(var_to_rawvar(s1, s2, variant, seq_list))
                #for i in seq_list:
                #    print i.dump()
            #if
        #for
    #else

    # Nothing happened.
    if not description:
        return Allele(RawVar())

    return description
#describe
