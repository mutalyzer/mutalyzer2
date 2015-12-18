"""
Tests for the mutalyzer.variantchecker module.
"""


from __future__ import unicode_literals

import pytest

from mutalyzer.variantchecker import check_variant

from fixtures import with_references


# Todo: We had a test for checking a variant on a CONTIG RefSeq reference
#   (NG_005990.1), but instead we should have separate tests for the retriever
#   module, including a test for fetching a CONTIG RefSeq reference.


@pytest.fixture
def checker(output):
    def check(description):
        check_variant(description, output)
    return check


@with_references('AL449423.14')
def test_deletion_in_frame(output, checker):
    """
    Simple in-frame deletion should give a simple description on protein
    level.
    """
    checker('AL449423.14(CDKN2A_v001):c.161_163del')
    assert (output.getIndexedOutput('genomicDescription', 0) ==
            'AL449423.14:g.61937_61939del')
    assert 'AL449423.14(CDKN2A_v001):c.161_163del' \
           in output.getOutput('descriptions')
    assert 'AL449423.14(CDKN2A_i001):p.(Met54_Gly55delinsSer)' \
           in output.getOutput('protDescriptions')
    assert output.getOutput('newProtein')


@with_references('AL449423.14')
def test_insertion_in_frame(output, checker):
    """
    Simple in-frame insertion should give a simple description on protein
    level.
    """
    checker('AL449423.14(CDKN2A_v001):c.161_162insATC')
    assert (output.getIndexedOutput('genomicDescription', 0) ==
            'AL449423.14:g.61938_61939insGAT')
    assert 'AL449423.14(CDKN2A_v001):c.161_162insATC' \
           in output.getOutput('descriptions')
    assert 'AL449423.14(CDKN2A_i001):p.(Met54delinsIleSer)' \
           in output.getOutput('protDescriptions')
    assert output.getOutput('newProtein')


@with_references('AL449423.14')
def test_insertion_list_in_frame(output, checker):
    """
    Simple in-frame insertion of a list should give a simple description
    on protein level.
    """
    checker('AL449423.14(CDKN2A_v001):c.161_162ins[ATC]')
    assert (output.getIndexedOutput('genomicDescription', 0) ==
            'AL449423.14:g.61938_61939insGAT')
    assert 'AL449423.14(CDKN2A_v001):c.161_162insATC' \
           in output.getOutput('descriptions')
    assert 'AL449423.14(CDKN2A_i001):p.(Met54delinsIleSer)' \
           in output.getOutput('protDescriptions')
    assert output.getOutput('newProtein')


@with_references('AL449423.14')
def test_deletion_insertion_in_frame(output, checker):
    """
    Simple in-frame deletion/insertion should give a simple description on
    protein level.
    """
    check_variant('AL449423.14(CDKN2A_v001):c.161_162delinsATCCC',
                  output)
    assert output.getIndexedOutput('genomicDescription', 0) == 'AL449423.14:g.61938_61939delinsGGGAT'
    assert 'AL449423.14(CDKN2A_v001):c.161_162delinsATCCC' \
           in output.getOutput('descriptions')
    assert 'AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)' \
           in output.getOutput('protDescriptions')
    assert output.getOutput('newProtein')


@with_references('AL449423.14')
def test_deletion_insertion_list_in_frame(output, checker):
    """
    Simple in-frame deletion-insertion of a list should give a simple
    description on protein level.
    """
    check_variant('AL449423.14(CDKN2A_v001):c.161_162delins[ATCCC]',
                  output)
    assert output.getIndexedOutput('genomicDescription', 0) == 'AL449423.14:g.61938_61939delinsGGGAT'
    assert 'AL449423.14(CDKN2A_v001):c.161_162delinsATCCC' \
           in output.getOutput('descriptions')
    assert 'AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)' \
           in output.getOutput('protDescriptions')
    assert output.getOutput('newProtein')


@with_references('AL449423.14')
def test_deletion_insertion_in_frame_complete(output, checker):
    """
    Simple in-frame deletion/insertion should give a simple description on
    protein level, also with the optional deleted sequence argument.
    """
    check_variant('AL449423.14(CDKN2A_v001):c.161_162delTGinsATCCC',
                  output)
    assert output.getIndexedOutput('genomicDescription', 0) == 'AL449423.14:g.61938_61939delinsGGGAT'
    assert 'AL449423.14(CDKN2A_v001):c.161_162delinsATCCC' \
           in output.getOutput('descriptions')
    assert 'AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)' \
           in output.getOutput('protDescriptions')
    assert output.getOutput('newProtein')


@with_references('AL449423.14')
def test_deletion_insertion_list_in_frame_complete(output, checker):
    """
    Simple in-frame deletion-insertion of a list should give a simple
    description on protein level, also with the optional deleted sequence
    argument.
    """
    check_variant('AL449423.14(CDKN2A_v001):c.161_162delTGins[ATCCC]',
                  output)
    assert output.getIndexedOutput('genomicDescription', 0) == 'AL449423.14:g.61938_61939delinsGGGAT'
    assert 'AL449423.14(CDKN2A_v001):c.161_162delinsATCCC' \
           in output.getOutput('descriptions')
    assert 'AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)' \
           in output.getOutput('protDescriptions')
    assert output.getOutput('newProtein')


@with_references('NM_003002.2')
def test_est_warning_nm_est(output, checker):
    """
    Warning for EST positioning on NM reference.
    """
    checker('NM_003002.2:274del')
    west = output.getMessagesWithErrorCode('WEST')
    assert len(west) == 1


@with_references('NM_003002.2')
def test_no_est_warning_nm_c(output, checker):
    """
    No EST warning for c. positioning on NM reference.
    """
    checker('NM_003002.2:c.274del')
    west = output.getMessagesWithErrorCode('WEST')
    assert len(west) == 0


@with_references('NM_003002.2')
def test_no_est_warning_nm_n(output, checker):
    """
    No EST warning for n. positioning on NM reference.
    """
    checker('NM_003002.2:n.274del')
    west = output.getMessagesWithErrorCode('WEST')
    assert len(west) == 0


@with_references('NG_012772.1')
def test_est_warning_ng_est(output, checker):
    """
    Warning for EST positioning on NG reference.
    """
    checker('NG_012772.1:128del')
    west = output.getMessagesWithErrorCode('WEST')
    assert len(west) == 1


@with_references('NG_012772.1')
def test_no_est_warning_ng_g(output, checker):
    """
    No EST warning for g. positioning on NG reference.
    """
    checker('NG_012772.1:g.128del')
    west = output.getMessagesWithErrorCode('WEST')
    assert len(west) == 0


@with_references('AA010203.1')
def test_no_est_warning_est_est(output, checker):
    """
    No warning for EST positioning on EST reference.
    """
    checker('AA010203.1:54_55insG')
    west = output.getMessagesWithErrorCode('WEST')
    assert len(west) == 0


@with_references('NM_003002.2')
def test_roll(output, checker):
    """
    Just a variant where we should roll.
    """
    checker('NM_003002.2:c.273del')
    wroll = output.getMessagesWithErrorCode('WROLLFORWARD')
    assert len(wroll) > 0


@with_references('NM_003002.2')
def test_no_roll(output, checker):
    """
    Just a variant where we cannot roll.
    """
    checker('NM_003002.2:c.274del')
    wroll = output.getMessagesWithErrorCode('WROLLFORWARD')
    assert len(wroll) == 0


@with_references('NM_000088.3')
def test_no_roll_splice(output, checker):
    """
    Here we can roll but should not, because it is over a splice site.
    """
    checker('NM_000088.3:g.459del')
    wrollback = output.getMessagesWithErrorCode('IROLLBACK')
    assert len(wrollback) > 0
    wroll = output.getMessagesWithErrorCode('WROLLFORWARD')
    assert len(wroll) == 0


@with_references('NM_000088.3')
def test_partial_roll_splice(output, checker):
    """
    Here we can roll two positions, but should roll only one because
    otherwise it is over a splice site.
    """
    checker('NM_000088.3:g.494del')
    wrollback = output.getMessagesWithErrorCode('IROLLBACK')
    assert len(wrollback) > 0
    wroll = output.getMessagesWithErrorCode('WROLLFORWARD')
    assert len(wroll) > 0


@with_references('NM_000088.3')
def test_roll_after_splice(output, checker):
    """
    Here we can roll and should, we stay in the same exon.
    """
    checker('NM_000088.3:g.460del')
    wroll = output.getMessagesWithErrorCode('WROLLFORWARD')
    assert len(wroll) > 0


@with_references('AL449423.14')
def test_roll_both_ins(output, checker):
    """
    Insertion that rolls should not use the same inserted sequence in
    descriptions on forward and reverse strands.

    Here we have the following situation on the forward strand:

                            65470 (genomic)
                              |
      CGGTGCGTTGGGCAGCGCCCCCGCCTCCAGCAGCGCCCGCACCTCCTCTA

    Now, an insertion of TAC after 65470 should be rolled to an insertion
    of ACT after 65471:

      CGGTGCGTTGGGCAGCGCCCCCGCC --- TCCAGCAGCGCCCGCACCTCCTCTA
      CGGTGCGTTGGGCAGCGCCCCCGCC TAC TCCAGCAGCGCCCGCACCTCCTCTA  =>

      CGGTGCGTTGGGCAGCGCCCCCGCCT --- CCAGCAGCGCCCGCACCTCCTCTA
      CGGTGCGTTGGGCAGCGCCCCCGCCT ACT CCAGCAGCGCCCGCACCTCCTCTA

    However, in CDKN2A_v001 (on the reverse strand), this insertion should
    roll the other direction and the inserted sequence should be the reverse
    complement of CTA, which is TAG, and not that of ACT, which is AGT.

    The next test (test_roll_reverse_ins) tests the situation for an input
    of AL449423.14:g.65471_65472insACT, where only the reverse roll should
    be done.
    """
    checker('AL449423.14:g.65470_65471insTAC')
    assert 'AL449423.14(CDKN2A_v001):c.99_100insTAG' in output.getOutput('descriptions')
    assert 'AL449423.14:g.65471_65472insACT' == output.getIndexedOutput('genomicDescription', 0, '')
    assert len(output.getMessagesWithErrorCode('WROLLFORWARD')) == 1


@with_references('AL449423.14')
def test_roll_reverse_ins(output, checker):
    """
    Insertion that rolls on the reverse strand should not use the same
    inserted sequence in descriptions on forward and reverse strands.
    """
    checker('AL449423.14:g.65471_65472insACT')
    assert 'AL449423.14(CDKN2A_v001):c.99_100insTAG' in output.getOutput('descriptions')
    assert 'AL449423.14:g.65471_65472insACT' == output.getIndexedOutput('genomicDescription', 0, '')
    assert len(output.getMessagesWithErrorCode('WROLLFORWARD')) == 0


@with_references('AL449423.14')
def test_roll_message_forward(output, checker):
    """
    Roll warning message should only be shown for currently selected
    strand (forward).
    """
    checker('AL449423.14:g.65470_65471insTAC')
    assert len(output.getMessagesWithErrorCode('WROLLFORWARD')) == 1
    assert len(output.getMessagesWithErrorCode('WROLLREVERSE')) == 0


@with_references('AL449423.14')
def test_roll_message_reverse(output, checker):
    """
    Roll warning message should only be shown for currently selected
    strand (reverse).
    """
    checker('AL449423.14(CDKN2A_v001):c.98_99insGTA')
    assert len(output.getMessagesWithErrorCode('WROLLFORWARD')) == 0
    assert len(output.getMessagesWithErrorCode('WROLLREVERSE')) == 1


@with_references('NM_000143.3')
def test_ins_cds_start(output, checker):
    """
    Insertion on CDS start boundary should not be included in CDS.
    """
    checker('NM_000143.3:c.-1_1insCAT')
    assert output.getIndexedOutput("newProtein", 0) is None
    # Todo: Is this a good test?


@with_references('NM_000143.3')
def test_ins_cds_start_after(output, checker):
    """
    Insertion after CDS start boundary should be included in CDS.
    """
    checker('NM_000143.3:c.1_2insCAT')
    assert output.getIndexedOutput("newProtein", 0) == '?'
    # Todo: Is this a good test?


@with_references('NG_012772.1')
def test_del_splice_site(output, checker):
    """
    Deletion hitting one splice site should not do a protein prediction.
    """
    checker('NG_012772.1(BRCA2_v001):c.632-5_670del')
    assert len(output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
    assert output.getOutput('removedSpliceSites') == []
    # Todo: For now, the following is how to check if no protein
    # prediction is done.
    assert not output.getOutput('newProtein')


@with_references('NG_012772.1')
def test_del_exon(output, checker):
    """
    Deletion of an entire exon should be possible.
    """
    checker('NG_012772.1(BRCA2_v001):c.632-5_681+7del')
    assert len(output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
    assert output.getOutput('removedSpliceSites') == [2]
    # Todo: For now, the following is how to check if protein
    # prediction is done.
    assert output.getOutput('newProtein')


@with_references('NG_012772.1')
def test_del_exon_exact(output, checker):
    """
    Deletion of exactly an exon should be possible.
    """
    checker('NG_012772.1(BRCA2_v001):c.632_681del')
    assert len(output.getMessagesWithErrorCode('WOVERSPLICE')) == 0
    assert output.getOutput('removedSpliceSites') == [2]
    # Todo: For now, the following is how to check if protein
    # prediction is done.
    assert output.getOutput('newProtein')


@with_references('NG_012772.1')
def test_del_exon_in_frame(output, checker):
    """
    Deletion of an entire exon with length a triplicate should give a
    proteine product with just this deletion (and possibly substitutions
    directly before and after).

    NG_012772.1(BRCA2_v001):c.68-7_316+7del is such a variant, since
    positions 68 through 316 are exactly one exon and (316-68+1)/3 = 83.
    """
    checker('NG_012772.1(BRCA2_v001):c.68-7_316+7del')
    assert len(output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
    assert output.getOutput('removedSpliceSites') == [2]
    # Todo: For now, the following is how to check if protein
    # prediction is done.
    assert output.getOutput('newProtein')
    # Todo: assert that protein products indeed have only this difference.


@with_references('NG_012772.1')
def test_del_exons(output, checker):
    """
    Deletion of two entire exons should be possible.
    """
    checker('NG_012772.1(BRCA2_v001):c.632-5_793+7del')
    assert len(output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
    assert output.getOutput('removedSpliceSites') == [4]
    # Todo: For now, the following is how to check if protein
    # prediction is done.
    assert output.getOutput('newProtein')


@with_references('NG_012772.1')
def test_del_intron(output, checker):
    """
    Deletion of an entire intron should be possible (fusion of remaining
    exonic parts).
    """
    checker('NG_012772.1(BRCA2_v001):c.622_674del')
    assert len(output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
    assert output.getOutput('removedSpliceSites') == [2]
    # Todo: For now, the following is how to check if protein
    # prediction is done.
    assert output.getOutput('newProtein')


@with_references('NG_012772.1')
def test_del_intron_exact(output, checker):
    """
    Deletion of exactly an intron should be possible (fusion of flanking
    exons).
    """
    checker('NG_012772.1(BRCA2_v001):c.681+1_682-1del')
    assert output.getMessagesWithErrorCode('WOVERSPLICE') == []
    assert output.getOutput('removedSpliceSites') == [2]
    # Note: The protein prediction is done, but 'newProtein' is not set
    # because we have no change. So to check if the prediction is done, we
    # check if 'oldProtein' is set and to check if the prediction is
    # correct, we check if 'newProtein' is not set.
    assert output.getOutput('oldProtein')
    assert not output.getOutput('newProtein')


@with_references('NG_012772.1')
def test_del_intron_in_frame(output, checker):
    """
    Deletion of an entire intron should be possible (fusion of remaining
    exonic parts).
    """
    checker('NG_012772.1(BRCA2_v001):c.622_672del')
    assert len(output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
    assert output.getOutput('removedSpliceSites') == [2]
    # Todo: For now, the following is how to check if protein
    # prediction is done.
    assert output.getOutput('newProtein')
    # Todo: assert that protein products indeed have only this difference.


@with_references('NG_012772.1')
def test_del_exon_unknown_offsets(output, checker):
    """
    Deletion of an entire exon with unknown offsets should be possible.
    """
    checker('NG_012772.1(BRCA2_v001):c.632-?_681+?del')
    assert len(output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
    assert len(output.getMessagesWithErrorCode('IDELSPLICE')) > 0
    # Todo: For now, the following is how to check if protein
    # prediction is done.
    assert output.getOutput('newProtein')
    # Genomic positions should be centered in flanking introns and unsure.
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_012772.1:g.(17550_19725)del'
    assert 'NG_012772.1(BRCA2_v001):c.632-?_681+?del' \
           in output.getOutput('descriptions')
    assert 'NG_012772.1(BRCA2_i001):p.(Val211Glufs*10)' \
           in output.getOutput('protDescriptions')
    # Todo: .c notation should still be c.632-?_681+?del, but what about
    # other transcripts?


@with_references('NG_012772.1')
def test_del_exon_unknown_offsets_in_frame(output, checker):
    """
    Deletion of an entire exon with unknown offsets and length a
    triplicate should give a proteine product with just this deletion
    (and possibly substitutions directly before and after).

    NG_012772.1(BRCA2_v001):c.68-?_316+?del is such a variant, since
    positions 68 through 316 are exactly one exon and (316-68+1)/3 = 83.
    """
    checker('NG_012772.1(BRCA2_v001):c.68-?_316+?del')
    assert len(output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
    assert len(output.getMessagesWithErrorCode('IDELSPLICE')) > 0
    # Todo: For now, the following is how to check if protein
    # prediction is done.
    assert output.getOutput('newProtein')
    # Genomic positions should be centered in flanking introns and unsure.
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_012772.1:g.(7324_11720)del'
    assert 'NG_012772.1(BRCA2_v001):c.68-?_316+?del' \
           in output.getOutput('descriptions')
    # Todo: .c notation should still be c.632-?_681+?del, but what about
    # other transcripts?


@with_references('NG_012772.1')
def test_del_exon_unknown_offsets_composed(output, checker):
    """
    Deletion of an entire exon with unknown offsets and another composed
    variant with exact positioning should be possible.
    """
    check_variant('NG_012772.1(BRCA2_v001):c.[632-?_681+?del;681+4del]',
                  output)
    assert len(output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
    assert len(output.getMessagesWithErrorCode('IDELSPLICE')) > 0
    # Todo: For now, the following is how to check if protein
    # prediction is done.
    assert output.getOutput('newProtein')
    # Genomic positions should be centered in flanking introns and unsure.
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_012772.1:g.[(17550_19725)del;19017del]'
    assert 'NG_012772.1(BRCA2_v001):c.[632-?_681+?del;681+4del]' \
           in output.getOutput('descriptions')
    # Todo: .c notation should still be c.632-?_681+?del, but what about
    # other transcripts?


@with_references('AL449423.14')
def test_del_exon_unknown_offsets_reverse(output, checker):
    """
    Deletion of an entire exon with unknown offsets should be possible,
    also on the reverse strand.
    """
    check_variant('AL449423.14(CDKN2A_v001):c.151-?_457+?del',
                  output)
    assert len(output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
    assert len(output.getMessagesWithErrorCode('IDELSPLICE')) > 0
    # Todo: For now, the following is how to check if protein
    # prediction is done.
    assert output.getOutput('newProtein')
    # Genomic positions should be centered in flanking introns and unsure.
    assert output.getIndexedOutput('genomicDescription', 0) == 'AL449423.14:g.(60314_63683)del'
    assert 'AL449423.14(CDKN2A_v001):c.151-?_457+?del' \
           in output.getOutput('descriptions')
    # Todo: .c notation should still be c.632-?_681+?del, but what about
    # other transcripts?


@with_references('NM_000143.3')
def test_del_exon_transcript_reference(output, checker):
    """
    Deletion of entire exon on a transcript reference should remove the
    expected splice sites (only that of the deleted exon), and not those
    of the flanking exons (as would happen using the mechanism for genomic
    references).
    """
    # checker('NM_018723.3:c.758_890del')
    checker('NM_000143.3:c.739_904del')
    assert len(output.getMessagesWithErrorCode('WOVERSPLICE')) == 0
    assert output.getOutput('removedSpliceSites') == [2]
    # Todo: For now, the following is how to check if protein
    # prediction is done.
    assert output.getOutput('newProtein')


@with_references('NG_008939.1')
def test_ins_seq(output, checker):
    """
    Insertion of a sequence.
    """
    checker('NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_012337.1')
def test_ins_seq_reverse(output, checker):
    """
    Insertion of a sequence on reverse strand.
    """
    checker('NG_012337.1(TIMM8B_v001):c.12_13insGATC')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_012337.1:g.4911_4912insATCG'
    assert 'NG_012337.1(TIMM8B_v001):c.12_13insGATC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_ins_range(output, checker):
    """
    Insertion of a range.
    """
    checker('NG_008939.1:g.5207_5208ins4300_4320')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 0


@with_references('NG_008939.1')
def test_ins_range_inv(output, checker):
    """
    Insertion of an inverse range.
    """
    checker('NG_008939.1:g.5207_5208ins4300_4320inv')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGCCAGATAATGAGCACAGGAC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGCCAGATAATGAGCACAGGAC' \
           in output.getOutput('descriptions')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 0


@with_references('NG_008939.1')
def test_ins_seq_list(output, checker):
    """
    Insertion of a sequence as a list.
    """
    checker('NG_008939.1:g.5207_5208ins[GTCCTGTGCTCATTATCTGGC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_012337.1')
def test_ins_seq_list_reverse(output, checker):
    """
    Insertion of a sequence as a list on reverse strand.
    """
    checker('NG_012337.1(TIMM8B_v001):c.12_13ins[GATC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_012337.1:g.4911_4912insATCG'
    assert 'NG_012337.1(TIMM8B_v001):c.12_13insGATC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_ins_range_list(output, checker):
    """
    Insertion of a range as a list.
    """
    checker('NG_008939.1:g.5207_5208ins[4300_4320]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 0


@with_references('NG_008939.1')
def test_ins_range_inv_list(output, checker):
    """
    Insertion of an inverse range as a list.
    """
    checker('NG_008939.1:g.5207_5208ins[4300_4320inv]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGCCAGATAATGAGCACAGGAC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGCCAGATAATGAGCACAGGAC' \
           in output.getOutput('descriptions')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 0


@with_references('NG_008939.1')
def test_ins_seq_seq(output, checker):
    """
    Insertion of two sequences.
    """
    checker('NG_008939.1:g.5207_5208ins[GTCCTGTGCTC;ATTATCTGGC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_012337.1')
def test_ins_seq_seq_reverse(output, checker):
    """
    Insertion of two sequences on reverse strand.
    """
    checker('NG_012337.1(TIMM8B_v001):c.12_13ins[TTT;GATC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_012337.1:g.4911_4912insATCAAAG'
    assert 'NG_012337.1(TIMM8B_v001):c.12_13insTTTGATC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_ins_range_range(output, checker):
    """
    Insertion of two ranges.
    """
    checker('NG_008939.1:g.5207_5208ins[4300_4309;4310_4320]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 0


@with_references('NG_008939.1')
def test_ins_range_range_inv(output, checker):
    """
    Insertion of a range and an inverse range.
    """
    checker('NG_008939.1:g.5207_5208ins[4300_4309;4310_4320inv]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTGCCAGATAATG'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTGCCAGATAATG' \
           in output.getOutput('descriptions')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 0


@with_references('NG_008939.1')
def test_ins_seq_range(output, checker):
    """
    Insertion of a sequence and a range.
    """
    checker('NG_008939.1:g.5207_5208ins[GTCCTGTGCT;4310_4320]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_ins_seq_range_inv(output, checker):
    """
    Insertion of a sequence and an inverse range.
    """
    checker('NG_008939.1:g.5207_5208ins[GTCCTGTGCT;4310_4320inv]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTGCCAGATAATG'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTGCCAGATAATG' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_ins_range_seq(output, checker):
    """
    Insertion of a range and a sequence.
    """
    checker('NG_008939.1:g.5207_5208ins[4300_4309;CATTATCTGGC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_ins_range_inv_seq(output, checker):
    """
    Insertion of an inverse range and a sequence.
    """
    checker('NG_008939.1:g.5207_5208ins[4300_4309inv;CATTATCTGGC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insAGCACAGGACCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insAGCACAGGACCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_ins_seq_coding(output, checker):
    """
    Insertion of a sequence (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_ins_seq_list_coding(output, checker):
    """
    Insertion of a sequence as a list (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_157ins[GTCCTGTGCTCATTATCTGGC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_ins_seq_seq_coding(output, checker):
    """
    Insertion of two sequences (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_157ins[GTCCTGTGCTC;ATTATCTGGC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_ins_range_coding(output, checker):
    """
    Insertion of a range (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_157ins180_188')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 1


@with_references('NG_008939.1')
def test_ins_range_inv_coding(output, checker):
    """
    Insertion of an inverse range (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_157ins180_188inv')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 1


@with_references('NG_008939.1')
def test_ins_range_list_coding(output, checker):
    """
    Insertion of a range as a list (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_157ins[180_188]')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 1


@with_references('NG_008939.1')
def test_ins_range_inv_list_coding(output, checker):
    """
    Insertion of an inverse range as a list (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_157ins[180_188inv]')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 1


@with_references('NG_008939.1')
def test_delins_seq(output, checker):
    """
    Insertion-deletion of a sequence.
    """
    checker('NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_delins_range(output, checker):
    """
    Insertion-deletion of a range.
    """
    checker('NG_008939.1:g.5207_5212delins4300_4320')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 0


@with_references('NG_008939.1')
def test_delins_range_inv(output, checker):
    """
    Insertion-deletion of an inverse range.
    """
    checker('NG_008939.1:g.5207_5212delins4300_4320inv')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGCCAGATAATGAGCACAGGAC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGCCAGATAATGAGCACAGGAC' \
           in output.getOutput('descriptions')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 0


@with_references('NG_008939.1')
def test_delins_seq_list(output, checker):
    """
    Insertion-deletion of a sequence as a list.
    """
    checker('NG_008939.1:g.5207_5212delins[GTCCTGTGCTCATTATCTGGC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_delins_range_list(output, checker):
    """
    Insertion-deletion of a range as a list.
    """
    checker('NG_008939.1:g.5207_5212delins[4300_4320]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 0


@with_references('NG_008939.1')
def test_delins_range_inv_list(output, checker):
    """
    Insertion-deletion of an inverse range as a list.
    """
    checker('NG_008939.1:g.5207_5212delins[4300_4320inv]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGCCAGATAATGAGCACAGGAC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGCCAGATAATGAGCACAGGAC' \
           in output.getOutput('descriptions')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 0


@with_references('NG_008939.1')
def test_delins_seq_seq(output, checker):
    """
    Insertion-deletion of two sequences.
    """
    checker('NG_008939.1:g.5207_5212delins[GTCCTGTGCT;CATTATCTGGC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_delins_range_range(output, checker):
    """
    Insertion-deletion of two ranges.
    """
    checker('NG_008939.1:g.5207_5212delins[4300_4309;4310_4320]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 0


@with_references('NG_008939.1')
def test_delins_range_inv_range(output, checker):
    """
    Insertion-deletion of an inverse range and a range.

    Note that the delins is also shortened by one position here.
    """
    checker('NG_008939.1:g.5207_5212delins[4300_4309inv;4310_4320]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5208_5212delinsGCACAGGACCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.157_161delinsGCACAGGACCATTATCTGGC' \
           in output.getOutput('descriptions')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 0


@with_references('NG_008939.1')
def test_delins_seq_range(output, checker):
    """
    Insertion-deletion of a sequence and a range.
    """
    checker('NG_008939.1:g.5207_5212delins[GTCCTGTGCT;4310_4320]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_delins_seq_range_inv(output, checker):
    """
    Insertion-deletion of a sequence and an inverse range.

    Note that the delins is also shortened by one position here.
    """
    checker('NG_008939.1:g.5207_5212delins[GTCCTGTGCT;4310_4320inv]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5211delinsGTCCTGTGCTGCCAGATAAT'
    assert 'NG_008939.1(PCCB_v001):c.156_160delinsGTCCTGTGCTGCCAGATAAT' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_delins_range_seq(output, checker):
    """
    Insertion-deletion of a range and a sequence.
    """
    checker('NG_008939.1:g.5207_5212delins[4300_4309;CATTATCTGGC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_delins_range_inv_seq(output, checker):
    """
    Insertion-deletion of an inverse range and a sequence.

    Note that the delins is also shortened by one position here.
    """
    checker('NG_008939.1:g.5207_5212delins[4300_4309inv;CATTATCTGGC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5208_5212delinsGCACAGGACCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.157_161delinsGCACAGGACCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_delins_seq_coding(output, checker):
    """
    Insertion-deletion of a sequence (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_delins_seq_list_coding(output, checker):
    """
    Insertion-deletion of a sequence as a list (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_161delins[GTCCTGTGCTCATTATCTGGC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_delins_seq_seq_coding(output, checker):
    """
    Insertion-deletion of two sequences (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_161delins[GTCCTGTGCT;CATTATCTGGC]')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC'
    assert 'NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_delins_range_coding(output, checker):
    """
    Insertion-deletion of a range (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_161delins180_188')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 1


@with_references('NG_008939.1')
def test_delins_range_inv_coding(output, checker):
    """
    Insertion-deletion of an inverse range (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_161delins180_188inv')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 1


@with_references('NG_008939.1')
def test_delins_range_list_coding(output, checker):
    """
    Insertion-deletion of a range as a list (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_161delins[180_188]')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 1


@with_references('NG_008939.1')
def test_delins_range_inv_list_coding(output, checker):
    """
    Insertion-deletion of an inverse range as a list (coding).
    """
    checker('NG_008939.1(PCCB_v001):c.156_161delins[180_188inv]')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 1


def test_no_reference(output, checker):
    """
    Variant description without a reference.
    """
    checker('g.244355733del')
    assert len(output.getMessagesWithErrorCode('ENOREF')) == 1


@pytest.mark.usefixtures('hg19_transcript_mappings')
@with_references('NM_003002.2')
def test_chromosomal_positions(output, checker):
    """
    Variants on transcripts in c. notation should have chromosomal positions
    defined.
    """
    checker('NM_003002.2:c.274G>T')
    assert output.getIndexedOutput('rawVariantsChromosomal', 0) == ('chr11', '+', [('274G>T', (111959695, 111959695))])


@with_references('NM_002001.2')
def test_ex_notation(output, checker):
    """
    Variant description using EX notation should not crash but deletion of
    one exon should delete two splice sites.
    """
    checker('NM_002001.2:c.EX1del')
    assert len(output.getMessagesWithErrorCode('IDELSPLICE')) == 1


@with_references('LRG_1')
def test_lrg_reference(output, checker):
    """
    We should be able to use LRG reference sequence without error.
    """
    checker('LRG_1t1:c.266G>T')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicDescription', 0) == 'LRG_1:g.6855G>T'


@with_references('NM_002001.2')
def test_gi_reference_plain(output, checker):
    """
    Test reference sequence notation with GI number.
    """
    checker('31317229:c.6del')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicDescription', 0) == '31317229:n.105del'
    assert '31317229(FCER1A_v001):c.6del' \
           in output.getOutput('descriptions')


@with_references('NM_002001.2')
def test_gi_reference_prefix(output, checker):
    """
    Test reference sequence notation with GI number and prefix.
    """
    checker('GI31317229:c.6del')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicDescription', 0) == '31317229:n.105del'
    assert '31317229(FCER1A_v001):c.6del' \
           in output.getOutput('descriptions')


@with_references('NM_002001.2')
def test_gi_reference_prefix_colon(output, checker):
    """
    Test reference sequence notation with GI number and prefix with colon.
    """
    checker('GI:31317229:c.6del')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicDescription', 0) == '31317229:n.105del'
    assert '31317229(FCER1A_v001):c.6del' \
           in output.getOutput('descriptions')


@with_references('NM_002001.2')
def test_nop_nm(output, checker):
    """
    Variant on NM without effect should be described as '='.
    """
    checker('NM_002001.2:c.1_3delinsATG')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicDescription', 0) == 'NM_002001.2:n.='
    assert 'NM_002001.2(FCER1A_v001):c.=' \
           in output.getOutput('descriptions')


@with_references('DMD')
def test_nop_ud(output, references, checker):
    """
    Variant on UD without effect should be described as '='.
    """
    ud = references[0].accession
    checker(ud + ':g.5T>T')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicChromDescription', 0) == 'NC_000023.11:g.='
    assert output.getIndexedOutput('genomicDescription', 0) == ud + ':g.='
    assert ud + '(DMD_v001):c.=' in output.getOutput('descriptions')


@with_references('DPYD')
def test_ud_reverse_sequence(output, references, checker):
    """
    Variant on UD from reverse strand should have reverse complement
    sequence.
    """
    ud = references[0].accession
    checker(ud + '(DPYD_v1):c.85C>T')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicChromDescription', 0) == 'NC_000001.10:g.98348885G>A'
    assert output.getIndexedOutput('genomicDescription', 0) == ud + ':g.42731C>T'
    assert ud + '(DPYD_v001):c.85C>T' in output.getOutput('descriptions')


@with_references('MARK1')
def test_ud_forward_sequence(output, references, checker):
    """
    Variant on UD from forward strand should have forward sequence.
    """
    ud = references[0].accession
    checker(ud + '(MARK1_v001):c.400T>C')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicChromDescription', 0) == 'NC_000001.10:g.220773181T>C'
    assert output.getIndexedOutput('genomicDescription', 0) == ud + ':g.76614T>C'
    assert ud + '(MARK1_v001):c.400T>C' in output.getOutput('descriptions')


@with_references('chr9_reverse')
def test_ud_reverse_range(output, references, checker):
    """
    Variant on UD from reverse strand should have reversed range
    positions.
    """
    # This is just some slice on from the reverse strand of hg19 chr9.
    ud = references[0].accession
    checker(ud + ':g.10624_78132del')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicChromDescription', 0) == 'NC_000009.11:g.32928508_32996016del'
    assert output.getIndexedOutput('genomicDescription', 0) == ud + ':g.10624_78132del'


@with_references('MARK1')
def test_ud_forward_range(output, references, checker):
    """
    Variant on UD from forward strand should have forward range positions.
    """
    ud = references[0].accession
    checker(ud + '(MARK1_v001):c.400_415del')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicChromDescription', 0) == 'NC_000001.10:g.220773181_220773196del'
    assert output.getIndexedOutput('genomicDescription', 0) == ud + ':g.76614_76629del'


@with_references('chr9_reverse')
def test_ud_reverse_del_length(output, references, checker):
    """
    Variant on UD from reverse strand should have reversed range
    positions, but not reverse complement of first argument (it is not a
    sequence, but a length).
    """
    # This is just some slice on from the reverse strand of hg19 chr9.
    ud = references[0].accession
    checker(ud + ':g.10624_78132del67509')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicChromDescription', 0) == 'NC_000009.11:g.32928508_32996016del'
    assert output.getIndexedOutput('genomicDescription', 0) == ud + ':g.10624_78132del'


@with_references('DPYD')
def test_ud_reverse_roll(output, references, checker):
    """
    Variant on UD from reverse strand should roll the oposite direction.

    The situation is as follows:

                  G    A    A    A    T    T
            c.   102  103  104  105  106  107
            g.   748  749  750  751  752  753
        chr g.   868  867  866  865  864  863
    """
    ud = references[0].accession
    checker(ud + '(DPYD_v001):c.104del')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicChromDescription', 0) == 'NC_000001.10:g.98348867del'
    assert output.getIndexedOutput('genomicDescription', 0) == ud + ':g.42751del'
    assert ud + '(DPYD_v001):c.105del' in output.getOutput('descriptions')


@with_references('MARK1')
def test_ud_forward_roll(output, references, checker):
    """
    Variant on UD from forward strand should roll the same.

    The situation is as follows:

                  A    T    T    T    A
            c.   398  399  400  401  402
            g.   612  613  614  615  616
        chr g.   179  180  181  182  183
    """
    ud = references[0].accession
    checker(ud + '(MARK1_v001):c.400del')
    error_count, _, _ = output.Summary()
    assert error_count == 0
    assert output.getIndexedOutput('genomicChromDescription', 0) == 'NC_000001.10:g.220773182del'
    assert output.getIndexedOutput('genomicDescription', 0) == ud + ':g.76615del'
    assert ud + '(MARK1_v001):c.401del' in output.getOutput('descriptions')


@with_references('AL449423.14')
def test_deletion_with_sequence_forward_genomic(output, checker):
    """
    Specify the deleted sequence in a deletion.
    """
    checker('AL449423.14:g.65471_65472delTC')
    assert output.getIndexedOutput('genomicDescription', 0) == 'AL449423.14:g.65471_65472del'
    assert 'AL449423.14(CDKN2A_v001):c.98_99del' \
           in output.getOutput('descriptions')


@with_references('AL449423.14')
def test_deletion_with_length_forward_genomic(output, checker):
    """
    Specify the deleted sequence length in a deletion.
    """
    checker('AL449423.14:g.65471_65472del2')
    assert output.getIndexedOutput('genomicDescription', 0) == 'AL449423.14:g.65471_65472del'
    assert 'AL449423.14(CDKN2A_v001):c.98_99del' \
           in output.getOutput('descriptions')


@with_references('AL449423.14')
def test_deletion_with_sequence_reverse_coding(output, checker):
    """
    Specify the deleted sequence in a deletion on the reverse strand.
    """
    checker('AL449423.14(CDKN2A_v001):c.161_163delTGG')
    assert output.getIndexedOutput('genomicDescription', 0) == 'AL449423.14:g.61937_61939del'
    assert 'AL449423.14(CDKN2A_v001):c.161_163del' \
           in output.getOutput('descriptions')


@with_references('AL449423.14')
def test_deletion_with_length_reverse_coding(output, checker):
    """
    Specify the deleted sequence length in a deletion on the reverse strand.
    """
    checker('AL449423.14(CDKN2A_v001):c.161_163del3')
    assert output.getIndexedOutput('genomicDescription', 0) == 'AL449423.14:g.61937_61939del'
    assert 'AL449423.14(CDKN2A_v001):c.161_163del' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_deletion_with_sequence_reverse_ng_coding(output, checker):
    """
    Specify the deleted sequence in a deletion on the reverse strand
    using a genomic reference.
    """
    checker('NG_008939.1:c.155_157delAAC')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5206_5208del'
    assert 'NG_008939.1(PCCB_v001):c.155_157del' \
           in output.getOutput('descriptions')


@with_references('NG_008939.1')
def test_deletion_with_length_reverse_ng_coding(output, checker):
    """
    Specify the deleted sequence length in a deletion on the reverse strand
    using a genomic reference.
    """
    checker('NG_008939.1:c.155_157del3')
    assert output.getIndexedOutput('genomicDescription', 0) == 'NG_008939.1:g.5206_5208del'
    assert 'NG_008939.1(PCCB_v001):c.155_157del' \
           in output.getOutput('descriptions')


@with_references('AB026906.1')
def test_inversion(output, checker):
    """
    Inversion variant.
    """
    checker('AB026906.1:c.274_275inv')
    assert output.getIndexedOutput('genomicDescription', 0) == 'AB026906.1:g.7872_7873inv'
    assert 'AB026906.1(SDHD_v001):c.274_275inv' \
        in output.getOutput('descriptions')


@with_references('NM_000193.2')
def test_delins_with_length(output, checker):
    """
    Delins with explicit length of deleted sequence (bug #108).
    """
    checker('NM_000193.2:c.108_109del2insG')
    assert 'NM_000193.2(SHH_i001):p.(Lys38Serfs*2)' in output.getOutput('protDescriptions')


@with_references('NG_009105.1')
def test_protein_level_description(output, checker):
    """
    Currently protein level descriptions are not implemented.
    """
    checker('NG_009105.1(OPN1LW):p.=')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 1


@with_references('NP_064445.1')
def test_protein_reference(output, checker):
    """
    Currently protein references are not implemented.
    """
    checker('NP_064445.1:p.=')
    assert len(output.getMessagesWithErrorCode('ENOTIMPLEMENTED')) == 1


@with_references('AF230870.1')
def test_wnomrna_other(output, checker):
    """
    Warning for no mRNA field on other than currently selected transcript
    should give WNOMRNA_OTHER warning.
    """
    # Contains mtmC2 and mtmB2, both without mRNA
    checker('AF230870.1(mtmC2_v001):c.13del')
    wnomrna_other = output.getMessagesWithErrorCode('WNOMRNA_OTHER')
    assert len(wnomrna_other) == 1


@with_references('AF230870.1')
def test_wnomrna(output, checker):
    """
    Warning for no mRNA field on currently selected transcript should give
    WNOMRNA warning.
    """
    # Contains mtmC2 and mtmB2, both without mRNA
    checker('AF230870.1(mtmC2_v001):c.13del')
    wnomrna = output.getMessagesWithErrorCode('WNOMRNA')
    wnomrna_other = output.getMessagesWithErrorCode('WNOMRNA_OTHER')
    assert len(wnomrna) == 1
    assert len(wnomrna_other) == 1


@with_references('L41870.1')
def test_mrna_ref_adjacent_exons_warn(output, checker):
    """
    Warning for mRNA reference where exons are not adjacent.

    In L41870.1 exon 15 ends on 1558 and 16 starts on 1636.
    """
    checker('L41870.1:c.1del')
    w_exon_annotation = output.getMessagesWithErrorCode('WEXON_ANNOTATION')
    assert len(w_exon_annotation) == 1


@with_references('NM_003002.2')
def test_mrna_ref_adjacent_exons_no_warn(output, checker):
    """
    No warning for mRNA reference where exons are adjacent.
    """
    checker('NM_003002.2:c.1del')
    w_exon_annotation = output.getMessagesWithErrorCode('WEXON_ANNOTATION')
    assert len(w_exon_annotation) == 0


@with_references('NM_001199.3')
def test_fs_no_stop(output, checker):
    """
    Frame shift yielding no stop codon should be described with
    uncertainty of the stop codon.

    http://www.hgvs.org/mutnomen/FAQ.html#nostop
    """
    checker('NM_001199.3(BMP1):c.2188dup')
    assert 'NM_001199.3(BMP1_i001):p.(Gln730Profs*?)' in output.getOutput('protDescriptions')


@with_references('NM_000193.2')
def test_ext_no_stop(output, checker):
    """
    Extension yielding no stop codon should be described with
    uncertainty of the stop codon.

    http://www.hgvs.org/mutnomen/FAQ.html#nostop
    """
    checker('NM_000193.2:c.1388G>C')
    assert 'NM_000193.2(SHH_i001):p.(*463Serext*?)' in output.getOutput('protDescriptions')


@with_references('NM_000193.2')
def test_fs_ext_no_stop(output, checker):
    """
    Extension yielding no stop codon should be described with
    uncertainty of the stop codon.

    http://www.hgvs.org/mutnomen/FAQ.html#nostop
    """
    checker('NM_000193.2:c.1388_1389insC')
    assert 'NM_000193.2(SHH_i001):p.(*463Cysext*?)' in output.getOutput('protDescriptions')


@with_references('AB026906.1')
def test_synonymous_p_is(output, checker):
    """
    Synonymous mutation should yield a p.(=) description.
    """
    checker('AB026906.1:c.276C>T')
    assert 'AB026906.1(SDHD_i001):p.(=)' in output.getOutput('protDescriptions')
    assert not output.getOutput('newProteinFancy')


@with_references('NM_024426.4')
def test_synonymous_p_is_alt_start(output, checker):
    """
    Synonymous mutation should yield a p.(=) description, also with an
    alternative start codon.
    """
    checker('NM_024426.4:c.1107A>G')
    assert 'NM_024426.4(WT1_i001):p.(=)' in output.getOutput('protDescriptions')
    assert not output.getOutput('newProteinFancy')
    waltstart = output.getMessagesWithErrorCode('WALTSTART')
    assert len(waltstart) == 1
    assert output.getOutput('oldProtein')[0].startswith('M')
    assert not output.getOutput('newProtein')
    assert not output.getOutput('altStart')
    assert not output.getOutput('altProteinFancy')


@with_references('AB026906.1')
def test_start_codon(output, checker):
    """
    Mutation of start codon should yield a p.? description.
    """
    checker('AB026906.1:c.1A>G')
    assert 'AB026906.1(SDHD_i001):p.?' in output.getOutput('protDescriptions')
    wstart = output.getMessagesWithErrorCode('WSTART')
    assert len(wstart) == 1
    assert output.getOutput('newProtein')[0] == '?'
    waltstart = output.getMessagesWithErrorCode('WALTSTART')
    assert len(waltstart) == 0
    assert not output.getOutput('altStart')


@with_references('NM_024426.4')
def test_start_codon_alt_start(output, checker):
    """
    Mutation of start codon should yield a p.? description, also with an
    alternative start codon.
    """
    checker('NM_024426.4:c.1C>G')
    assert 'NM_024426.4(WT1_i001):p.?' in output.getOutput('protDescriptions')
    west = output.getMessagesWithErrorCode('WSTART')
    assert len(west) == 1
    assert output.getOutput('newProtein')[0] == '?'
    waltstart = output.getMessagesWithErrorCode('WALTSTART')
    assert len(waltstart) == 1
    assert not output.getOutput('altStart')


@with_references('AB026906.1')
def test_start_codon_yield_start_p_is(output, checker):
    """
    Silent mutation creating new start codon should yield a p.?
    description. The visualisation should also render the case for the new
    start codon.
    """
    checker('AB026906.1:c.1A>T')  # yields TTG start codon
    assert 'AB026906.1(SDHD_i001):p.?' in output.getOutput('protDescriptions')
    wstart = output.getMessagesWithErrorCode('WSTART')
    assert len(wstart) == 1
    assert output.getOutput('newProtein')[0] == '?'
    waltstart = output.getMessagesWithErrorCode('WALTSTART')
    assert len(waltstart) == 0
    assert output.getOutput('oldProtein')[0].startswith('M')
    assert 'TTG' in output.getOutput('altStart')
    assert not output.getOutput('altProteinFancy')


@with_references('NM_024426.4')
def test_start_codon_alt_start_yield_start_p_is(output, checker):
    """
    Silent mutation creating new start codon should yield a p.?
    description, also with an alternative start codon. The visualisation
    should also render the case for the new start codon.
    """
    checker('NM_024426.4:c.1C>A')  # yields ATG start codon
    assert 'NM_024426.4(WT1_i001):p.?' in output.getOutput('protDescriptions')
    west = output.getMessagesWithErrorCode('WSTART')
    assert len(west) == 1
    assert output.getOutput('newProtein')[0] == '?'
    waltstart = output.getMessagesWithErrorCode('WALTSTART')
    assert len(waltstart) == 1
    assert output.getOutput('oldProtein')[0].startswith('M')
    assert 'ATG' in output.getOutput('altStart')
    assert not output.getOutput('altProteinFancy')


@with_references('AB026906.1')
def test_start_codon_yield_start(output, checker):
    """
    Mutation creating new start codon should yield a p.? description. The
    visualisation should also render the case for the new start codon.
    """
    checker('AB026906.1:c.1_4delinsTTGA')  # yields TTG start codon
    assert 'AB026906.1(SDHD_i001):p.?' in output.getOutput('protDescriptions')
    wstart = output.getMessagesWithErrorCode('WSTART')
    assert len(wstart) == 1
    assert output.getOutput('newProtein')[0] == '?'
    waltstart = output.getMessagesWithErrorCode('WALTSTART')
    assert len(waltstart) == 0
    assert 'TTG' in output.getOutput('altStart')
    assert output.getOutput('altProtein')[0].startswith('M')


@with_references('NM_024426.4')
def test_start_codon_alt_start_yield_start(output, checker):
    """
    Mutation creating new start codon should yield a p.? description, also
    with an alternative start codon. The visualisation should also render
    the new start codon.
    """
    checker('NM_024426.4:c.1_4delinsATGA')  # yields ATG start codon
    assert 'NM_024426.4(WT1_i001):p.?' in output.getOutput('protDescriptions')
    west = output.getMessagesWithErrorCode('WSTART')
    assert len(west) == 1
    assert output.getOutput('newProtein')[0] == '?'
    waltstart = output.getMessagesWithErrorCode('WALTSTART')
    assert len(waltstart) == 1
    assert output.getOutput('oldProtein')[0].startswith('M')
    assert 'ATG' in output.getOutput('altStart')
    assert output.getOutput('altProtein')[0].startswith('M')


@with_references('AB026906.1')
def test_legend_mrna_by_construction(output, checker):
    """
    Transcript created from CDS by construction should be in the legend.
    """
    checker('AB026906.1:g.7872G>T')
    assert output.getOutput('legends') == [
        ['SDHD_v001', None, None, None, 'construction'],
        ['SDHD_i001', 'BAA81889.1', None, 'small subunit of cytochrome b of succinate dehydrogenase', 'construction']
    ]
