"""
Generate a HGVS description of the variant(s) leading from one sequence to an
other.
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
    #       visualisation in the _visualise() function of mutator.py
    fs = 10 # Flank size.

    return '{} {}{} {}'.format(s[start - fs:start], s[start:end], '-' * fill,
        s[end:end + fs])


def var_to_rawvar(s1, s2, var, seq_list=[], container=DNAVar,
        weight_position=1):
    """
    Convert a variant from the extractor module to one of the RawVar
    subclasses.

    :arg unicode s1: Reference sequence.
    :arg unicode s2: Sample sequence.
    :arg str var: Variant from the extractor module.
    :arg str seq_list: Container for an inserted sequence.
    :arg str container: Destination container.
    :arg str weight_position: Weight of a position.
    """
    # Unknown.
    if s1 == '?' or s2 == '?':
        return [container(type='unknown', weight_position=weight_position)]

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
                end=var.reference_end, type='dup', shift=shift,
                sample_start=var.sample_start + 1, sample_end=var.sample_end,
                inserted=ISeqList([ISeq(sequence=s2[
                var.sample_start:var.sample_end],
                    weight_position=weight_position)]),
                weight_position=weight_position)

        return container(start=var.reference_start,
            end=var.reference_start + 1,
            inserted=seq_list or
            ISeqList([ISeq(sequence=s2[var.sample_start:var.sample_end],
                weight_position=weight_position)]),
            type='ins', shift=shift, sample_start=var.sample_start + 1,
            sample_end=var.sample_end, weight_position=weight_position)

    # Deletion.
    if var.sample_start == var.sample_end:
        shift5, shift3 = roll(s1, var.reference_start + 1, var.reference_end)
        shift = shift5 + shift3

        var.reference_start += shift3
        var.reference_end += shift3

        return container(start=var.reference_start + 1,
            end=var.reference_end, type='del', shift=shift,
            sample_start=var.sample_start, sample_end=var.sample_end + 1,
            deleted=ISeqList([ISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
            weight_position=weight_position)

    # Substitution.
    if (var.reference_start + 1 == var.reference_end and
            var.sample_start + 1 == var.sample_end):
        return container(start=var.reference_start + 1,
            end=var.reference_end, sample_start=var.sample_start + 1,
            sample_end=var.sample_end, type='subst',
            deleted=ISeqList([ISeq(sequence=s1[var.reference_start],
                weight_position=weight_position)]),
            inserted=ISeqList([ISeq(sequence=s2[var.sample_start],
                weight_position=weight_position)]),
            weight_position=weight_position)

    # Inversion.
    if var.type & extractor.REVERSE_COMPLEMENT:
        trim = palinsnoop(s1[var.reference_start:var.reference_end])

        if trim > 0: # Partial palindrome.
            var.reference_end -= trim
            var.sample_end -= trim

        return container(start=var.reference_start + 1,
            end=var.reference_end, type='inv',
            sample_start=var.sample_start + 1, sample_end=var.sample_end,
            deleted=ISeqList([ISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
            inserted=ISeqList([ISeq(sequence=s2[
                var.sample_start:var.reference_end],
                weight_position=weight_position)]),
            weight_position=weight_position)

    # InDel.
    return container(start=var.reference_start + 1,
        end=var.reference_end, deleted=ISeqList([ISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
        inserted=seq_list or
        ISeqList([ISeq(sequence=s2[var.sample_start:var.sample_end],
            weight_position=weight_position)]),
        type='delins', sample_start=var.sample_start + 1,
        sample_end=var.sample_end, weight_position=weight_position)


def describe_dna(s1, s2):
    """
    Give an allele description of the change from {s1} to {s2}.

    :arg unicode s1: Sequence 1.
    :arg unicode s2: Sequence 2.

    :returns list(RawVar): A list of RawVar objects, representing the allele.
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
        elif not (variant.type & extractor.IDENTITY):
            description.append(var_to_rawvar(s1, s2, variant,
                weight_position=extracted.weight_position))

        if variant.type & extractor.TRANSPOSITION_CLOSE:
            in_transposition -= 1

            if not in_transposition:
                description.append(var_to_rawvar(s1, s2, variant, seq_list,
                    weight_position=extracted.weight_position))

    if not description:
        return Allele([DNAVar()])
    return description


def describe_protein(s1, s2):
    """
    Give an allele description of the change from {s1} to {s2}.

    :arg unicode s1: Sequence 1.
    :arg unicode s2: Sequence 2.

    :returns list(RawVar): A list of RawVar objects, representing the allele.
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
    else:
        print s1[:len(s1) - longest_fs_r[0]], s1[len(s1) - longest_fs_r[0]:]
        print s2[:longest_fs_r[1]], s2[longest_fs_r[1]:]
        s1_part = s1[:len(s1) - longest_fs_r[0]]
        s2_part = s2[:longest_fs_r[1]]
        term = len(s2) - longest_fs_r[1]

    s1_part = s1
    s2_part = s2
    for variant in extractor.extract(s1_part.encode('utf-8'), len(s1_part),
            s2_part.encode('utf-8'), len(s2_part), 1):
        description.append(var_to_rawvar(s1, s2, variant,
            container=ProteinVar))

    if description:
        description[-1].term = term + 2

    return description
