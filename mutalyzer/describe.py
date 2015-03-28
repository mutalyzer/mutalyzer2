"""
Prototype of a module that can generate a HGVS description of the variant(s)
leading from one sequence to an other.

@requires: Bio.Seq
"""


from __future__ import unicode_literals

import collections

from Bio.Data import CodonTable

from mutalyzer.util import palinsnoop, roll
from mutalyzer.variant import ISeq, ISeqList, DNAVar, ProteinVar, Allele

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


def var_to_rawvar(s1, s2, var, seq_list=[], container=DNAVar,
        weight_position=1):
    """
    """
    # Unknown.
    if s1 == '?' or s2 == '?':
        return [container(type="unknown", weight_position=weight_position)]

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

            # NOTE: We may want to omit the inserted / deleted sequence and
            # use the ranges instead.
            return container(start=var.reference_start - ins_length + 1,
                end=var.reference_end, type="dup", shift=shift,
                sample_start=var.sample_start + 1, sample_end=var.sample_end,
                inserted=ISeqList([ISeq(sequence=s2[
                var.sample_start:var.sample_end],
                    weight_position=weight_position)]),
                weight_position=weight_position)
        #if
        return container(start=var.reference_start,
            end=var.reference_start + 1,
            inserted=seq_list or
            ISeqList([ISeq(sequence=s2[var.sample_start:var.sample_end],
                weight_position=weight_position)]),
            type="ins", shift=shift, sample_start=var.sample_start + 1,
            sample_end=var.sample_end, weight_position=weight_position)
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
            deleted=ISeqList([ISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
            weight_position=weight_position)
    #if

    # Substitution.
    if (var.reference_start + 1 == var.reference_end and
        var.sample_start + 1 == var.sample_end):

        return container(start=var.reference_start + 1,
            end=var.reference_end, sample_start=var.sample_start + 1,
            sample_end=var.sample_end, type="subst",
            deleted=ISeqList([ISeq(sequence=s1[var.reference_start],
                weight_position=weight_position)]),
            inserted=ISeqList([ISeq(sequence=s2[var.sample_start],
                weight_position=weight_position)]),
            weight_position=weight_position)
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
            deleted=ISeqList([ISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
            inserted=ISeqList([ISeq(sequence=s2[
                var.sample_start:var.reference_end],
                weight_position=weight_position)]),
            weight_position=weight_position)
    #if

    # InDel.
    return container(start=var.reference_start + 1,
        end=var.reference_end, deleted=ISeqList([ISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
        inserted=seq_list or
        ISeqList([ISeq(sequence=s2[var.sample_start:var.sample_end],
            weight_position=weight_position)]),
        type="delins", sample_start=var.sample_start + 1,
        sample_end=var.sample_end, weight_position=weight_position)
#var_to_rawvar

def describe_dna(s1, s2):
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
    in_transposition = 0

    extracted = extractor.extract(s1.encode('utf-8'), len(s1),
                                  s2.encode('utf-8'), len(s2), 0)
    for variant in extracted.variants:
       # print (variant.type, variant.reference_start,
       #     variant.reference_end, variant.sample_start,
       #     variant.sample_end, variant.transposition_start,
       #     variant.transposition_end)
       # print (variant.type & extractor.TRANSPOSITION_OPEN, variant.type &
       #     extractor.TRANSPOSITION_CLOSE)

        if variant.type & extractor.TRANSPOSITION_OPEN:
            if not in_transposition:
                seq_list = ISeqList()
            in_transposition += 1
        #if

        if in_transposition:
            if variant.type & extractor.IDENTITY:
                seq_list.append(ISeq(start=variant.transposition_start + 1,
                    end=variant.transposition_end, reverse=False,
                    weight_position=extracted.weight_position))
            elif variant.type & extractor.REVERSE_COMPLEMENT:
                seq_list.append(ISeq(start=variant.transposition_start + 1,
                    end=variant.transposition_end, reverse=True,
                    weight_position=extracted.weight_position))
            else:
                seq_list.append(ISeq(
                    sequence=s2[variant.sample_start:variant.sample_end],
                    weight_position=extracted.weight_position))
        #if
        elif not (variant.type & extractor.IDENTITY):
            description.append(var_to_rawvar(s1, s2, variant,
                weight_position=extracted.weight_position))

        if variant.type & extractor.TRANSPOSITION_CLOSE:
            in_transposition -= 1

            if not in_transposition:
                description.append(var_to_rawvar(s1, s2, variant, seq_list,
                    weight_position=extracted.weight_position))
        #if
    #for

    if not description:
        return Allele([DNAVar()])
    return description
#describe_dna

def describe_protein(s1, s2):
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

    fs1, fs2 = make_fs_tables(1)
    longest_fs_f = max(find_fs(s1, s2, fs1), find_fs(s1, s2, fs2))
    longest_fs_r = max(find_fs(s2, s1, fs1), find_fs(s2, s1, fs2))

    if longest_fs_f > longest_fs_r:
        print s1[:longest_fs_f[1]], s1[longest_fs_f[1]:]
        print s2[:len(s2) - longest_fs_f[0]], s2[len(s2) - longest_fs_f[0]:]
        s1_part = s1[:longest_fs_f[1]]
        s2_part = s2[:len(s2) - longest_fs_f[0]]
        term = longest_fs_f[0]
    #if
    else:
        print s1[:len(s1) - longest_fs_r[0]], s1[len(s1) - longest_fs_r[0]:]
        print s2[:longest_fs_r[1]], s2[longest_fs_r[1]:]
        s1_part = s1[:len(s1) - longest_fs_r[0]]
        s2_part = s2[:longest_fs_r[1]]
        term = len(s2) - longest_fs_r[1]
    #else

    s1_part = s1
    s2_part = s2
    for variant in extractor.extract(s1_part.encode('utf-8'), len(s1_part),
                                     s2_part.encode('utf-8'), len(s2_part), 1):
        description.append(var_to_rawvar(s1, s2, variant,
            container=ProteinVar))

    if description:
        description[-1].term = term + 2

    return description
#describe_protein
