"""
Tests for the mutalyzer.mutator module.
"""


from __future__ import unicode_literals

import pytest
import random
from Bio.Seq import Seq

from mutalyzer.mutator import Mutator


@pytest.fixture
def length():
    # Many tests depend on this being at least 30.
    return 30


@pytest.fixture
def sequence(length):
    return Seq(''.join(random.choice('ACGT') for _ in range(length)))


@pytest.fixture
def mutator(output, sequence):
    return Mutator(sequence, output)


def test_shift_no_change(length, mutator):
    """
    No change, no shifts.
    """
    # Numbering is 1-based
    for i in range(1, length + 1):
        assert mutator.shift(i) == i


def test_shift_del_example(mutator):
    """
    Example of g.2del.
    """
    mutator.deletion(2, 2)
    assert mutator.shift(1) == 1
    assert mutator.shift(2) == 2
    assert mutator.shift(3) == 2


@pytest.mark.parametrize('length', [10])
@pytest.mark.parametrize('d', range(1, 11))
def test_shift_del(length, mutator, d):
    """
    Starting from the deleted position (not included), shift -1.
    """
    mutator.deletion(d, d)
    for p in range(1, d + 1):
        assert mutator.shift(p) == p
    for p in range(d + 1, length + 1):
        assert mutator.shift(p) == p - 1


@pytest.mark.parametrize('length', [10])
@pytest.mark.parametrize('d', range(1, 10))
def test_shift_del2(length, mutator, d):
    """
    Starting from the deleted positions (not included), shift -2.
    """
    mutator.deletion(d, d + 1)
    for p in range(1, d + 2):
        assert mutator.shift(p) == p
    for p in range(d + 2, length + 1):
        assert mutator.shift(p) == p - 2


def test_shift_ins_example(mutator):
    """
    Example of g.2_3insA.
    """
    mutator.insertion(2, 'A')
    assert mutator.shift(1) == 1
    assert mutator.shift(2) == 2
    assert mutator.shift(3) == 4


@pytest.mark.parametrize('length', [10])
@pytest.mark.parametrize('i', range(11))
def test_shift_ins(length, mutator, i):
    """
    Starting from the interbase insertion position, shift +1.
    """
    mutator.insertion(i, 'T')
    for p in range(1, i + 1):
        assert mutator.shift(p) == p
    for p in range(i + 1, length + 1):
        assert mutator.shift(p) == p + 1


@pytest.mark.parametrize('length', [10])
@pytest.mark.parametrize('i', range(11))
def test_shift_ins2(length, mutator, i):
    """
    Starting from the interbase insertion position, shift +2.
    """
    mutator.insertion(i, 'TT')
    for p in range(1, i + 1):
        assert mutator.shift(p) == p
    for p in range(i + 1, length + 1):
        assert mutator.shift(p) == p + 2


def test_shift_sites_no_change(mutator):
    """
    No change, no shifts.

    Note: Splice sites come in pairs (acceptor and donor site) and the
    numbers are the first, respectively last, position in the exon.

    So in this example we have:     ---======----======-----===---
                                       |    |    |    |     | |
                                       4    9   14    19   25 27
    """
    sites = [4, 9, 14, 19, 25, 27]
    assert mutator.shift_sites(sites) == sites


def test_shift_sites_acc_del_before(mutator):
    """
    Deletion in intron directly before exon.

    @note: This hits a splice site, so we don't really support it.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.deletion(13, 13)   # g.13del
    assert mutator.shift_sites(sites) == [4, 9, 13, 16, 24, 26]


def test_shift_sites_acc_del_after(mutator):
    """
    Deletion at first exon position.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.deletion(14, 14)   # g.14del
    assert mutator.shift_sites(sites) == [4, 9, 14, 16, 24, 26]


def test_shift_sites_don_del_before(mutator):
    """
    Deletion at last exon position.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.deletion(17, 17)   # g.17del
    assert mutator.shift_sites(sites) == [4, 9, 14, 16, 24, 26]


def test_shift_sites_don_del_after(mutator):
    """
    Deletion in intron directly after exon.

    @note: This hits a splice site, so we don't really support it.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.deletion(18, 18)   # g.18del
    assert mutator.shift_sites(sites) == [4, 9, 14, 17, 24, 26]


def test_shift_sites_acc_del2_before(mutator):
    """
    Deletion of 2 in intron directly before exon.

    @note: This hits a splice site, so we don't really support it.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.deletion(12, 13)   # g.12_13del
    assert mutator.shift_sites(sites) == [4, 9, 12, 15, 23, 25]


def test_shift_sites_acc_del2_on(mutator):
    """
    Deletion of 2 in intron/exon.

    @note: This hits a splice site, so we don't really support it.
    """
    return

    sites = [4, 9, 14, 17, 25, 27]
    mutator.deletion(13, 14)   # g.13_14del
    assert mutator.shift_sites(sites) == [4, 9, 13, 15, 23, 25]


def test_shift_sites_acc_del2_after(mutator):
    """
    Deletion of 2 at first exon position.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.deletion(14, 15)   # g.14_15del
    assert mutator.shift_sites(sites) == [4, 9, 14, 15, 23, 25]


def test_shift_sites_don_del2_before(mutator):
    """
    Deletion of 2 at last exon positions.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.deletion(16, 17)   # g.16_17del
    assert mutator.shift_sites(sites) == [4, 9, 14, 15, 23, 25]


def test_shift_sites_don_del2_on(mutator):
    """
    Deletion of 2 in exon/intron.

    @note: This hits a splice site, so we don't really support it.
    """
    return

    sites = [4, 9, 14, 17, 25, 27]
    mutator.deletion(17, 18)   # g.17_18del
    assert mutator.shift_sites(sites) == [4, 9, 14, 16, 23, 25]


def test_shift_sites_don_del2_after(mutator):
    """
    Deletion of 2 in intron directly after exon.

    @note: This hits a splice site, so we don't really support it.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.deletion(18, 19)   # g.18_19del
    assert mutator.shift_sites(sites) == [4, 9, 14, 17, 23, 25]


def test_shift_sites_acc_ins_before(mutator):
    """
    Insertion 1 position before intron/exon boundary.

    @note: This hits a splice site, so we don't really support it.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(12, 'A')   # g.12_13insA
    assert mutator.shift_sites(sites) == [4, 9, 15, 18, 26, 28]


def test_shift_sites_acc_ins_on(mutator):
    """
    Insertion in intron/exon boundary.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(13, 'A')   # g.13_14insA
    assert mutator.shift_sites(sites) == [4, 9, 14, 18, 26, 28]


def test_shift_sites_first_acc_ins_on(mutator):
    """
    Insertion in first intron/exon boundary not be included.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(3, 'A')   # g.3_4insA
    assert mutator.shift_sites(sites) == [5, 10, 15, 18, 26, 28]


def test_shift_sites_acc_ins_after(mutator):
    """
    Insertion 1 position after intron/exon boundary.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(14, 'A')   # g.14_15insA
    assert mutator.shift_sites(sites) == [4, 9, 14, 18, 26, 28]


def test_shift_sites_don_ins_before(mutator):
    """
    Insertion 1 position before exon/intron boundary.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(16, 'A')   # g.16_17insA
    assert mutator.shift_sites(sites) == [4, 9, 14, 18, 26, 28]


def test_shift_sites_don_ins_on(mutator):
    """
    Insertion in exon/intron boundary.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(17, 'A')   # g.17_18insA
    assert mutator.shift_sites(sites) == [4, 9, 14, 18, 26, 28]


def test_shift_sites_last_don_ins_on(mutator):
    """
    Insertion in last exon/intron boundary should not be included.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(27, 'A')   # g.27_28insA
    assert mutator.shift_sites(sites) == [4, 9, 14, 17, 25, 27]


def test_shift_sites_don_ins_after(mutator):
    """
    Insertion 1 position after exon/intron boundary.

    @note: This hits a splice site, so we don't really support it.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(18, 'A')   # g.18_19insA
    assert mutator.shift_sites(sites) == [4, 9, 14, 17, 26, 28]


def test_shift_sites_acc_ins2_before(mutator):
    """
    Insertion of 2 1 position before intron/exon boundary.

    @note: This hits a splice site, so we don't really support it.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(12, 'AT')   # g.12_13insAT
    assert mutator.shift_sites(sites) == [4, 9, 16, 19, 27, 29]


def test_shift_sites_first_acc_ins2_on(mutator):
    """
    Insertion of 2 in last exon/intron boundary should not be included.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(3, 'AT')   # g.3_4insAT
    assert mutator.shift_sites(sites) == [6, 11, 16, 19, 27, 29]


def test_shift_sites_acc_ins2_after(mutator):
    """
    Insertion of 2 1 position after intron/exon boundary.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(14, 'AT')   # g.14_15insAT
    assert mutator.shift_sites(sites) == [4, 9, 14, 19, 27, 29]


def test_shift_sites_don_ins2_before(mutator):
    """
    Insertion of 2 1 position before exon/intron boundary.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(16, 'AT')   # g.16_17insAT
    assert mutator.shift_sites(sites) == [4, 9, 14, 19, 27, 29]


def test_shift_sites_last_don_ins2_on(mutator):
    """
    Insertion of 2 in last exon/intron boundary should not be included.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(27, 'AT')   # g.27_28insAT
    assert mutator.shift_sites(sites) == [4, 9, 14, 17, 25, 27]


def test_shift_sites_don_ins2_after(mutator):
    """
    Insertion of 2 1 position after exon/intron boundary.

    @note: This hits a splice site, so we don't really support it.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(18, 'AT')   # g.18_19insAT
    assert mutator.shift_sites(sites) == [4, 9, 14, 17, 27, 29]


def test_shift_sites_acc_ins3_before(mutator):
    """
    Insertion of 3 1 position before intron/exon boundary.

    @note: This hits a splice site, so we don't really support it.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(12, 'ATT')   # g.12_13insATT
    assert mutator.shift_sites(sites) == [4, 9, 17, 20, 28, 30]


def test_shift_sites_acc_ins3_on(mutator):
    """
    Insertion of 3 in intron/exon boundary.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(13, 'ATT')   # g.13_14insATT
    assert mutator.shift_sites(sites) == [4, 9, 14, 20, 28, 30]


def test_shift_sites_first_acc_ins3_on(mutator):
    """
    Insertion of 3 in first intron/exon boundary should not be included.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(3, 'ATT')   # g.3_4insATT
    assert mutator.shift_sites(sites) == [7, 12, 17, 20, 28, 30]


def test_shift_sites_acc_ins3_after(mutator):
    """
    Insertion of 3 1 position after intron/exon boundary.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(14, 'ATT')   # g.14_15insATT
    assert mutator.shift_sites(sites) == [4, 9, 14, 20, 28, 30]


def test_shift_sites_don_ins3_before(mutator):
    """
    Insertion of 3 1 position before exon/intron boundary.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(16, 'ATT')   # g.16_17insATT
    assert mutator.shift_sites(sites) == [4, 9, 14, 20, 28, 30]


def test_shift_sites_don_ins3_on(mutator):
    """
    Insertion of 3 in exon/intron boundary.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(17, 'ATT')   # g.17_18insATT
    assert mutator.shift_sites(sites) == [4, 9, 14, 20, 28, 30]


def test_shift_sites_last_don_ins3_on(mutator):
    """
    Insertion of 3 in last exon/intron boundary should not be included.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(27, 'ATT')   # g.27_28insATT
    assert mutator.shift_sites(sites) == [4, 9, 14, 17, 25, 27]


def test_shift_sites_don_ins3_after(mutator):
    """
    Insertion of 3 1 position after exon/intron boundary.

    @note: This hits a splice site, so we don't really support it.
    """
    sites = [4, 9, 14, 17, 25, 27]
    mutator.insertion(18, 'ATT')   # g.18_19insATT
    assert mutator.shift_sites(sites) == [4, 9, 14, 17, 28, 30]


def test_shift_sites_adj_del_before1(mutator):
    """
    Adjacent exons: deletion at second-last position of first exon.

    @note: In this example we have adjacent exons (like e.g. in RNA),
    which looks like this (the square brackets [ and ] are part of the
    exons):
                ---[====][======][========]---
                   |   /  \    /  \       |
                   4  9   10  17  18      27
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.deletion(16, 16)   # g.16del
    assert mutator.shift_sites(sites) == [4, 9, 10, 16, 17, 26]


def test_shift_sites_adj_del_before(mutator):
    """
    Adjacent exons: deletion at last position of first exon.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.deletion(17, 17)   # g.17del
    assert mutator.shift_sites(sites) == [4, 9, 10, 16, 17, 26]


def test_shift_sites_adj_del_after(mutator):
    """
    Adjacent exons: deletion at first position of second exon.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.deletion(18, 18)   # g.18del
    assert mutator.shift_sites(sites) == [4, 9, 10, 17, 18, 26]


def test_shift_sites_adj_del_after1(mutator):
    """
    Adjacent exons: deletion at second position of second exon.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.deletion(19, 19)   # g.19del
    assert mutator.shift_sites(sites) == [4, 9, 10, 17, 18, 26]


def test_shift_sites_adj_ins_before(mutator):
    """
    Adjacent exons: insertion 1 position before exon/exon boundary.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.insertion(16, 'A')   # g.16_17insA
    assert mutator.shift_sites(sites) == [4, 9, 10, 18, 19, 28]


def test_shift_sites_adj_ins_on(mutator):
    """
    Adjacent exons: insertion at exon/exon boundary.

    @note: This insertion could be seen as being
           1) at the end of the first exon, or
           2) at the start of the second exon.
           Both would probably be 'correct', but we would like consistent
           results. Therefore, we stick to the first option.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.insertion(17, 'A')   # g.17_18insA
    assert mutator.shift_sites(sites) == [4, 9, 10, 18, 19, 28]


def test_shift_sites_adj_ins_after(mutator):
    """
    Adjacent exons: insertion 1 position after exon/exon boundary.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.insertion(18, 'A')   # g.18_19insA
    assert mutator.shift_sites(sites) == [4, 9, 10, 17, 18, 28]


def test_shift_sites_adj_del2_before1(mutator):
    """
    Adjacent exons: deletion of 2 at second-last position of first exon.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.deletion(15, 16)   # g.15_16del
    assert mutator.shift_sites(sites) == [4, 9, 10, 15, 16, 25]


def test_shift_sites_adj_del2_before(mutator):
    """
    Adjacent exons: deletion of 2 at last position of first exon.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.deletion(16, 17)   # g.16_17del
    assert mutator.shift_sites(sites) == [4, 9, 10, 15, 16, 25]


def test_shift_sites_adj_del2_on(mutator):
    """
    Adjacent exons: deletion of 2 at exon/exon boundary.

    @todo: This is a special case of bug #????. Once fixed, the two
           exons will be joined to one new exon.
    """
    return

    sites = [4, 9, 10, 17, 18, 27]
    mutator.deletion(17, 18)   # g.17_18del
    assert mutator.shift_sites(sites) == [4, 9, 10, 16, 17, 25]


def test_shift_sites_adj_del2_after(mutator):
    """
    Adjacent exons: deletion of 2 at first position of second exon.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.deletion(18, 19)   # g.18_19del
    assert mutator.shift_sites(sites) == [4, 9, 10, 17, 18, 25]


def test_shift_sites_adj_del2_after1(mutator):
    """
    Adjacent exons: deletion of 2 at second position of second exon.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.deletion(19, 20)   # g.19_20del
    assert mutator.shift_sites(sites) == [4, 9, 10, 17, 18, 25]


def test_shift_sites_adj_ins2_before(mutator):
    """
    Adjacent exons: insertion of 2 1 position before exon/exon boundary.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.insertion(16, 'AT')   # g.16_17insAT
    assert mutator.shift_sites(sites) == [4, 9, 10, 19, 20, 29]


def test_shift_sites_adj_ins2_on(mutator):
    """
    Adjacent exons: insertion of 2 at exon/exon boundary.

    @note: This insertion could be seen as being
           1) at the end of the first exon, or
           2) at the start of the second exon.
           Both would probably be 'correct', but we would like consistent
           results. Therefore, we stick to the first option.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.insertion(17, 'AT')   # g.17_18insAT
    assert mutator.shift_sites(sites) == [4, 9, 10, 19, 20, 29]


def test_shift_sites_adj_ins2_after(mutator):
    """
    Adjacent exons: insertion of 2 1 position after exon/exon boundary.
    """
    sites = [4, 9, 10, 17, 18, 27]
    mutator.insertion(18, 'AT')   # g.18_19insAT
    assert mutator.shift_sites(sites) == [4, 9, 10, 17, 18, 29]


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_del(mutator):
    """
    Simple deletion 2del.
    """
    mutator.deletion(2, 2)
    assert unicode(mutator.mutated) == unicode(Seq('ACGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_largedel(mutator):
    """
    Simple large deletion 2_7del.
    """
    mutator.deletion(2, 7)
    assert unicode(mutator.mutated) == unicode(Seq('AG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_ins(mutator):
    """
    Simple insertion 2_3insA.
    """
    mutator.insertion(2, 'A')
    assert unicode(mutator.mutated) == unicode(Seq('ATACGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_largeins(mutator):
    """
    Simple large insertion 2_3insATCG.
    """
    mutator.insertion(2, 'ATCG')
    assert unicode(mutator.mutated) == unicode(Seq('ATATCGCGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_sub(mutator):
    """
    Simple substitution 3C>G.
    """
    mutator.substitution(3, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('ATGGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_del_sub_1(mutator):
    """
    Deletion and substitution directly adjecent to each other [2del;3C>G].

    See Trac #83.
    """
    mutator.deletion(2, 2)
    mutator.substitution(3, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AGGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_del_sub_2(mutator):
    """
    Deletion and substitution directly adjecent to each other [3del;2T>G].
    """
    mutator.deletion(3, 3)
    mutator.substitution(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AGGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_near_adjecent_del_sub_1(mutator):
    """
    Deletion and substitution almost adjecent to each other [2del;4G>T].
    """
    mutator.deletion(2, 2)
    mutator.substitution(4, 'T')
    assert unicode(mutator.mutated) == unicode(Seq('ACTATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_near_adjecent_del_sub_2(mutator):
    """
    Deletion and substitution almost adjecent to each other [4del;2T>G].
    """
    mutator.deletion(4, 4)
    mutator.substitution(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AGCATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedel_sub_1(mutator):
    """
    Large deletion and substitution directly adjecent to each other
    [2_6del;7C>T].
    """
    mutator.deletion(2, 6)
    mutator.substitution(7, 'T')
    assert unicode(mutator.mutated) == unicode(Seq('ATG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedel_sub_2(mutator):
    """
    Large deletion and substitution directly adjecent to each other
    [3_7del;2T>C].
    """
    mutator.deletion(3, 7)
    mutator.substitution(2, 'C')
    assert unicode(mutator.mutated) == unicode(Seq('ACG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_near_adjecent_largedel_sub_1(mutator):
    """
    Large deletion and substitution almost adjecent to each other [2_5del;7C>T].
    """
    mutator.deletion(2, 5)
    mutator.substitution(7, 'T')
    assert unicode(mutator.mutated) == unicode(Seq('ATTG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_near_adjecent_largedel_sub_2(mutator):
    """
    Large deletion and substitution almost adjecent to each other [4_7del;2T>C].
    """
    mutator.deletion(4, 7)
    mutator.substitution(2, 'C')
    assert unicode(mutator.mutated) == unicode(Seq('ACCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_del_ins_1(mutator):
    """
    Deletion and insertion adjecent to each other [2del;2_3insG].
    """
    mutator.deletion(2, 2)
    mutator.insertion(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AGCGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_del_ins_2(mutator):
    """
    Deletion and insertion adjecent to each other [3del;2_3insA].
    """
    mutator.deletion(3, 3)
    mutator.insertion(2, 'A')
    assert unicode(mutator.mutated) == unicode(Seq('ATAGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_near_adjectent_del_ins(mutator):
    """
    Deletion and insertion almost adjecent to each other [2del;3_4insG].
    """
    mutator.deletion(2, 2)
    mutator.insertion(3, 'T')
    assert unicode(mutator.mutated) == unicode(Seq('ACTGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_ins_sub_1(mutator):
    """
    Insertion and substitution directly adjecent to each other
    [2_3insA;3C>G].
    """
    mutator.insertion(2, 'A')
    mutator.substitution(3, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('ATAGGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_ins_sub_2(mutator):
    """
    Insertion and substitution directly adjecent to each other
    [2_3insA;2T>G].
    """
    mutator.insertion(2, 'A')
    mutator.substitution(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AGACGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_near_adjecent_ins_sub(mutator):
    """
    Insertion and substitution almost adjecent to each other
    [2_3insA;4C>T].
    """
    mutator.insertion(2, 'A')
    mutator.substitution(4, 'T')
    assert unicode(mutator.mutated) == unicode(Seq('ATACTATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largeins_sub_1(mutator):
    """
    Large insertion and substitution directly adjecent to each other
    [2_3insATCG;3C>G].
    """
    mutator.insertion(2, 'ATCG')
    mutator.substitution(3, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('ATATCGGGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largeins_sub_2(mutator):
    """
    Large insertion and substitution directly adjecent to each other
    [2_3insATCG;2T>G].
    """
    mutator.insertion(2, 'ATCG')
    mutator.substitution(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AGATCGCGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_near_adjecent_largeins_sub(mutator):
    """
    Large insertion and substitution almost adjecent to each other
    [2_3insATCG;4C>T].
    """
    mutator.insertion(2, 'ATCG')
    mutator.substitution(4, 'T')
    assert unicode(mutator.mutated) == unicode(Seq('ATATCGCTATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_del_del_1(mutator):
    """
    Deletion and deletion directly adjecent to each other [2del;3del].
    """
    mutator.deletion(2, 2)
    mutator.deletion(3, 3)
    assert unicode(mutator.mutated) == unicode(Seq('AGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_del_del_2(mutator):
    """
    Deletion and deletion directly adjecent to each other [3del;2del].
    """
    mutator.deletion(3, 3)
    mutator.deletion(2, 2)
    assert unicode(mutator.mutated) == unicode(Seq('AGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_delins_snp_1(mutator):
    """
    Delins and deletion directly adjecent to each other [2delinsA;3C>G].
    """
    mutator.delins(2, 2, 'A')
    mutator.substitution(3, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AAGGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_delins_snp_2(mutator):
    """
    Delins and deletion directly adjecent to each other [3delinsA;2T>G].
    """
    mutator.delins(3, 3, 'A')
    mutator.substitution(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AGAGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedelins_eq_snp_1(mutator):
    """
    Large delins and deletion directly adjecent to each other
    [2_6delinsAAAAA;7C>G].
    """
    mutator.delins(2, 6, 'AAAAA')
    mutator.substitution(7, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AAAAAAGG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedelins_min_snp_1(mutator):
    """
    Large delins (min) and deletion directly adjecent to each other
    [2_6delinsAAA;7C>G].
    """
    mutator.delins(2, 6, 'AAA')
    mutator.substitution(7, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AAAAGG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedelins_plus_snp_1(mutator):
    """
    Large delins (plus) and deletion directly adjecent to each other
    [2_6delinsAAAAAAA;7C>G].
    """
    mutator.delins(2, 6, 'AAAAAAA')
    mutator.substitution(7, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AAAAAAAAGG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedelins_eq_snp_2(mutator):
    """
    Large delins and deletion directly adjecent to each other
    [3_7delinsAAAAA;2T>G].
    """
    mutator.delins(3, 7, 'AAAAA')
    mutator.substitution(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AGAAAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedelins_min_snp_2(mutator):
    """
    Large delins (min) and deletion directly adjecent to each other
    [3_7delinsAAA;2T>G].
    """
    mutator.delins(3, 7, 'AAA')
    mutator.substitution(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AGAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedelins_plus_snp_2(mutator):
    """
    Large delins (plus) and deletion directly adjecent to each other
    [3_7delinsAAAAAAA;2T>G].
    """
    mutator.delins(3, 7, 'AAAAAAA')
    mutator.substitution(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AGAAAAAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_delins_del_1(mutator):
    """
    Delins and deletion directly adjecent to each other [2delinsA;3del].
    """
    mutator.delins(2, 2, 'A')
    mutator.deletion(3, 3)
    assert unicode(mutator.mutated) == unicode(Seq('AAGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_delins_del_2(mutator):
    """
    Delins and deletion directly adjecent to each other [3delinsA;2del].
    """
    mutator.delins(3, 3, 'A')
    mutator.deletion(2, 2)
    assert unicode(mutator.mutated) == unicode(Seq('AAGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedelins_eq_del_1(mutator):
    """
    Large delins and deletion directly adjecent to each other
    [2_6delinsAAAAA;7del].
    """
    mutator.delins(2, 6, 'AAAAA')
    mutator.deletion(7, 7)
    assert unicode(mutator.mutated) == unicode(Seq('AAAAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedelins_min_del_1(mutator):
    """
    Large delins (min) and deletion directly adjecent to each other
    [2_6delinsAAA;7del].
    """
    mutator.delins(2, 6, 'AAA')
    mutator.deletion(7, 7)
    assert unicode(mutator.mutated) == unicode(Seq('AAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedelins_plus_del_1(mutator):
    """
    Large delins (plus) and deletion directly adjecent to each other
    [2_6delinsAAAAAAA;7del].
    """
    mutator.delins(2, 6, 'AAAAAAA')
    mutator.deletion(7, 7)
    assert unicode(mutator.mutated) == unicode(Seq('AAAAAAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedelins_eq_del_2(mutator):
    """
    Large delins and deletion directly adjecent to each other
    [3_7delinsAAAAA;2del].
    """
    mutator.delins(3, 7, 'AAAAA')
    mutator.deletion(2, 2)
    assert unicode(mutator.mutated) == unicode(Seq('AAAAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedelins_min_del_2(mutator):
    """
    Large delins (min) and deletion directly adjecent to each other
    [3_7delinsAAA;2del].
    """
    mutator.delins(3, 7, 'AAA')
    mutator.deletion(2, 2)
    assert unicode(mutator.mutated) == unicode(Seq('AAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_largedelins_plus_del_2(mutator):
    """
    Large delins (plus) and deletion directly adjecent to each other
    [3_7delinsAAAAAAA;2del].
    """
    mutator.delins(3, 7, 'AAAAAAA')
    mutator.deletion(2, 2)
    assert unicode(mutator.mutated) == unicode(Seq('AAAAAAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_delins_ins_1(mutator):
    """
    Delins and insertion adjecent to each other [2delinsA;2_3insG].
    """
    mutator.delins(2, 2, 'A')
    mutator.insertion(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AAGCGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_delins_ins_2(mutator):
    """
    Delins and insertion adjecent to each other [3delinsA;2_3insG].
    """
    mutator.delins(3, 3, 'A')
    mutator.insertion(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('ATGAGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_largedelins_eq_ins_1(mutator):
    """
    Large delins and insertion adjecent to each other [2_6delinsAAAAA;6_7insG].
    """
    mutator.delins(2, 6, 'AAAAA')
    mutator.insertion(6, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AAAAAAGCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_largedelins_min_ins_1(mutator):
    """
    Large delins (min) and insertion adjecent to each other [2_6delinsAAA;6_7insG].
    """
    mutator.delins(2, 6, 'AAA')
    mutator.insertion(6, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AAAAGCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_largedelins_plus_ins_1(mutator):
    """
    Large delins (plus) and insertion adjecent to each other [2_6delinsAAAAAAA;6_7insG].
    """
    mutator.delins(2, 6, 'AAAAAAA')
    mutator.insertion(6, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('AAAAAAAAGCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_largedelins_eq_ins_2(mutator):
    """
    Large delins and insertion adjecent to each other [3_7delinsAAAAA;2_3insG].
    """
    mutator.delins(3, 7, 'AAAAA')
    mutator.insertion(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('ATGAAAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_largedelins_min_ins_2(mutator):
    """
    Large delins (min) and insertion adjecent to each other [3_7delinsAAA;2_3insG].
    """
    mutator.delins(3, 7, 'AAA')
    mutator.insertion(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('ATGAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_largedelins_plus_ins_2(mutator):
    """
    Large delins (plus) and insertion adjecent to each other [3_7delinsAAAAAAA;2_3insG].
    """
    mutator.delins(3, 7, 'AAAAAAA')
    mutator.insertion(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('ATGAAAAAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_delins_del_delins(mutator):
    """
    Delins (deletion) and delins (SNP) adjecent to each other [2_3delinsA;4delinsT].
    """
    mutator.delins(2, 3, 'A')
    mutator.delins(4, 4, 'T')
    assert unicode(mutator.mutated) == unicode(Seq('AATATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_largedelins_plus_delins_1(mutator):
    """
    Large delins (plus) and delins adjecent to each other [2_6delinsAAAAAAA;7delinsT].
    """
    mutator.delins(2, 6, 'AAAAAAA')
    mutator.delins(7, 7, 'T')
    assert unicode(mutator.mutated) == unicode(Seq('AAAAAAAATG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_largedelins_plus_delins_2(mutator):
    """
    Large delins (plus) and delins adjecent to each other [3_7delinsAAAAAAA;2delinsC].
    """
    mutator.delins(3, 7, 'AAAAAAA')
    mutator.delins(2, 2, 'C')
    assert unicode(mutator.mutated) == unicode(Seq('ACAAAAAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_largedelins_min_delins_1(mutator):
    """
    Large delins (min) and delins adjecent to each other [2_6delinsAAA;7delinsT].
    """
    mutator.delins(2, 6, 'AAA')
    mutator.delins(7, 7, 'T')
    assert unicode(mutator.mutated) == unicode(Seq('AAAATG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_largedelins_min_delins_2(mutator):
    """
    Large delins (min) and delins adjecent to each other [3_7delinsAAA;2delinsC].
    """
    mutator.delins(3, 7, 'AAA')
    mutator.delins(2, 2, 'C')
    assert unicode(mutator.mutated) == unicode(Seq('ACAAAG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_del_dup_1(mutator):
    """
    Deletion and duplication adjecent to each other [2del;3dup].
    """
    mutator.deletion(2, 2)
    mutator.duplication(3, 3)
    assert unicode(mutator.mutated) == unicode(Seq('ACCGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_del_dup_2(mutator):
    """
    Deletion and duplication adjecent to each other [3del;2dup].
    """
    mutator.deletion(3, 3)
    mutator.duplication(2, 2)
    assert unicode(mutator.mutated) == unicode(Seq('ATTGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_ins_dup_1(mutator):
    """
    Insertion and duplication adjecent to each other [2_3insG;3dup].
    """
    mutator.insertion(2, 'G')
    mutator.duplication(3, 3)
    assert unicode(mutator.mutated) == unicode(Seq('ATGCCGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_ins_dup_2(mutator):
    """
    Insertion and duplication adjecent to each other [2_3insG;2dup].
    """
    mutator.insertion(2, 'G')
    mutator.duplication(2, 2)
    assert unicode(mutator.mutated) == unicode(Seq('ATTGCGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_ins_ins_1(mutator):
    """
    Insertion and insertion adjecent to each other [2_3insG;3_4insA].
    """
    mutator.insertion(2, 'G')
    mutator.insertion(3, 'A')
    assert unicode(mutator.mutated) == unicode(Seq('ATGCAGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjectent_ins_ins_2(mutator):
    """
    Insertion and insertion adjecent to each other [3_4insA;2_3insG].
    """
    mutator.insertion(3, 'A')
    mutator.insertion(2, 'G')
    assert unicode(mutator.mutated) == unicode(Seq('ATGCAGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_ins_ins(mutator):
    """
    Insertion and insertion at same position [2_3insG;2_3insA].
    """
    mutator.insertion(2, 'G')
    mutator.insertion(2, 'A')
    assert unicode(mutator.mutated) in (unicode(Seq('ATGACGATCG')), unicode(Seq('ATAGCGATCG')))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_inv_inv_1(mutator):
    """
    Inversion and inversion directly adjecent to each other [2inv;3inv].
    """
    mutator.inversion(2, 2)
    mutator.inversion(3, 3)
    assert unicode(mutator.mutated) == unicode(Seq('AAGGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_inv_inv_2(mutator):
    """
    Inversion and inversion directly adjecent to each other [3inv;2inv].
    """
    mutator.inversion(3, 3)
    mutator.inversion(2, 2)
    assert unicode(mutator.mutated) == unicode(Seq('AAGGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_dup_dup_1(mutator):
    """
    Duplication and duplication directly adjecent to each other [2dup;3dup].
    """
    mutator.duplication(2, 2)
    mutator.duplication(3, 3)
    assert unicode(mutator.mutated) == unicode(Seq('ATTCCGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_dup_dup_2(mutator):
    """
    Duplication and duplication directly adjecent to each other [3dup;2dup].
    """
    mutator.duplication(3, 3)
    mutator.duplication(2, 2)
    assert unicode(mutator.mutated) == unicode(Seq('ATTCCGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_del_inv_1(mutator):
    """
    Deletion and inversion directly adjecent to each other [2del;3inv].
    """
    mutator.deletion(2, 2)
    mutator.inversion(3, 3)
    assert unicode(mutator.mutated) == unicode(Seq('AGGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_del_inv_2(mutator):
    """
    Deletion and inversion directly adjecent to each other [3del;2inv].
    """
    mutator.deletion(3, 3)
    mutator.inversion(2, 2)
    assert unicode(mutator.mutated) == unicode(Seq('AAGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_ins_inv_1(mutator):
    """
    Insertion and inversion directly adjecent to each other [2_3insG;3inv].
    """
    mutator.insertion(2, 'G')
    mutator.inversion(3, 3)
    assert unicode(mutator.mutated) == unicode(Seq('ATGGGATCG'))


@pytest.mark.parametrize('sequence', [Seq('ATCGATCG')])
def test_adjecent_ins_inv_2(mutator):
    """
    Insertion and inversion directly adjecent to each other [2_3insG;2inv].
    """
    mutator.insertion(2, 'G')
    mutator.inversion(2, 2)
    assert unicode(mutator.mutated) == unicode(Seq('AAGCGATCG'))
