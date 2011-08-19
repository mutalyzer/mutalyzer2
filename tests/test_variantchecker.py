"""
Tests for the variantchecker module.
"""


#import logging; logging.basicConfig()
import re
import os
import random
from nose.tools import *
from Bio.Seq import Seq

import mutalyzer
from mutalyzer.config import Config
from mutalyzer.output import Output
from mutalyzer.variantchecker import check_variant


class TestVariantchecker():
    """
    Test the variantchecker module.
    """
    def setUp(self):
        """
        Initialize test variantchecker module.
        """
        self.config = Config()
        self.output = Output(__file__, self.config.Output)

    def test_deletion_in_frame(self):
        """
        Simple in-frame deletion should give a simple description on protein
        level.
        """
        check_variant('AL449423.14(CDKN2A_v001):c.161_163del',
                      self.config, self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'AL449423.14:g.61937_61939del')
        assert 'AL449423.14(CDKN2A_v001):c.161_163del' \
               in self.output.getOutput('descriptions')
        assert 'AL449423.14(CDKN2A_i001):p.(Met54_Gly55delinsSer)' \
               in self.output.getOutput('protDescriptions')
        assert self.output.getOutput('newprotein')

    def test_insertion_in_frame(self):
        """
        Simple in-frame insertion should give a simple description on protein
        level.
        """
        check_variant('AL449423.14(CDKN2A_v001):c.161_162insATC',
                      self.config, self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'AL449423.14:g.61938_61939insGAT')
        assert 'AL449423.14(CDKN2A_v001):c.161_162insATC' \
               in self.output.getOutput('descriptions')
        assert 'AL449423.14(CDKN2A_i001):p.(Met54delinsIleSer)' \
               in self.output.getOutput('protDescriptions')
        assert self.output.getOutput('newprotein')

    def test_deletion_insertion_in_frame(self):
        """
        Simple in-frame deletion/insertion should give a simple description on
        protein level.
        """
        check_variant('AL449423.14(CDKN2A_v001):c.161_162delinsATCCC',
                      self.config, self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'AL449423.14:g.61938_61939delinsGGGAT')
        assert 'AL449423.14(CDKN2A_v001):c.161_162delinsATCCC' \
               in self.output.getOutput('descriptions')
        assert 'AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)' \
               in self.output.getOutput('protDescriptions')
        assert self.output.getOutput('newprotein')

    def test_deletion_insertion_in_frame_complete(self):
        """
        Simple in-frame deletion/insertion should give a simple description on
        protein level, also with the optional deleted sequence argument.
        """
        check_variant('AL449423.14(CDKN2A_v001):c.161_162delTGinsATCCC',
                      self.config, self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'AL449423.14:g.61938_61939delinsGGGAT')
        assert 'AL449423.14(CDKN2A_v001):c.161_162delinsATCCC' \
               in self.output.getOutput('descriptions')
        assert 'AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)' \
               in self.output.getOutput('protDescriptions')
        assert self.output.getOutput('newprotein')

    def test_roll(self):
        """
        Just a variant where we should roll.
        """
        check_variant('NM_003002.2:c.273del', self.config, self.output)
        wroll = self.output.getMessagesWithErrorCode('WROLLFORWARD')
        assert len(wroll) > 0

    def test_no_roll(self):
        """
        Just a variant where we cannot roll.
        """
        check_variant('NM_003002.2:c.274del', self.config, self.output)
        wroll = self.output.getMessagesWithErrorCode('WROLLFORWARD')
        assert_equal(len(wroll), 0)

    def test_no_roll_splice(self):
        """
        Here we can roll but should not, because it is over a splice site.
        """
        check_variant('NM_000088.3:g.459del', self.config, self.output)
        wrollback = self.output.getMessagesWithErrorCode('IROLLBACK')
        assert len(wrollback) > 0
        wroll = self.output.getMessagesWithErrorCode('WROLLFORWARD')
        assert_equal(len(wroll), 0)

    def test_partial_roll_splice(self):
        """
        Here we can roll two positions, but should roll only one because
        otherwise it is over a splice site.
        """
        check_variant('NM_000088.3:g.494del', self.config, self.output)
        wrollback = self.output.getMessagesWithErrorCode('IROLLBACK')
        assert len(wrollback) > 0
        wroll = self.output.getMessagesWithErrorCode('WROLLFORWARD')
        assert len(wroll) > 0

    def test_roll_after_splice(self):
        """
        Here we can roll and should, we stay in the same exon.
        """
        check_variant('NM_000088.3:g.460del', self.config, self.output)
        wroll = self.output.getMessagesWithErrorCode('WROLLFORWARD')
        assert len(wroll) > 0

    def test_roll_both_ins(self):
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
        check_variant('AL449423.14:g.65470_65471insTAC', self.config, self.output)
        assert 'AL449423.14(CDKN2A_v001):c.99_100insTAG' in self.output.getOutput('descriptions')
        assert_equal ('AL449423.14:g.65471_65472insACT', self.output.getIndexedOutput('genomicDescription', 0, ''))
        assert_equal(len(self.output.getMessagesWithErrorCode('WROLLFORWARD')), 1)

    def test_roll_reverse_ins(self):
        """
        Insertion that rolls on the reverse strand should not use the same
        inserted sequence in descriptions on forward and reverse strands.
        """
        check_variant('AL449423.14:g.65471_65472insACT', self.config, self.output)
        assert 'AL449423.14(CDKN2A_v001):c.99_100insTAG' in self.output.getOutput('descriptions')
        assert_equal ('AL449423.14:g.65471_65472insACT', self.output.getIndexedOutput('genomicDescription', 0, ''))
        assert_equal(len(self.output.getMessagesWithErrorCode('WROLLFORWARD')), 0)

    def test_roll_message_forward(self):
        """
        Roll warning message should only be shown for currently selected
        strand (forward).
        """
        check_variant('AL449423.14:g.65470_65471insTAC', self.config, self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('WROLLFORWARD')), 1)
        assert_equal(len(self.output.getMessagesWithErrorCode('WROLLREVERSE')), 0)

    def test_roll_message_reverse(self):
        """
        Roll warning message should only be shown for currently selected
        strand (reverse).
        """
        check_variant('AL449423.14(CDKN2A_v001):c.98_99insGTA', self.config, self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('WROLLFORWARD')), 0)
        assert_equal(len(self.output.getMessagesWithErrorCode('WROLLREVERSE')), 1)

    def test_ins_cds_start(self):
        """
        Insertion on CDS start boundary should not be included in CDS.
        """
        check_variant('NM_000143.3:c.-1_1insCAT', self.config, self.output)
        assert_equal(self.output.getIndexedOutput("newprotein", 0), None)
        # Todo: Is this a good test?

    def test_ins_cds_start_after(self):
        """
        Insertion after CDS start boundary should be included in CDS.
        """
        check_variant('NM_000143.3:c.1_2insCAT', self.config, self.output)
        assert_equal(self.output.getIndexedOutput("newprotein", 0), '?')
        # Todo: Is this a good test?

    def test_del_splice_site(self):
        """
        Deletion hitting one splice site should not do a protein prediction.
        """
        check_variant('NG_012772.1(BRCA2_v001):c.632-5_670del',
                      self.config, self.output)
        assert len(self.output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
        assert_equal(self.output.getOutput('removedSpliceSites'), [])
        # Todo: For now, the following is how to check if no protein
        # prediction is done.
        assert not self.output.getOutput('newprotein')

    def test_del_exon(self):
        """
        Deletion of an entire exon should be possible.
        """
        check_variant('NG_012772.1(BRCA2_v001):c.632-5_681+7del',
                      self.config, self.output)
        assert len(self.output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
        assert_equal(self.output.getOutput('removedSpliceSites'), [2])
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')

    def test_del_exon_exact(self):
        """
        Deletion of exactly an exon should be possible.
        """
        check_variant('NG_012772.1(BRCA2_v001):c.632_681del',
                      self.config, self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('WOVERSPLICE')), 0)
        assert_equal(self.output.getOutput('removedSpliceSites'), [2])
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')

    def test_del_exon_in_frame(self):
        """
        Deletion of an entire exon with length a triplicate should give a
        proteine product with just this deletion (and possibly substitutions
        directly before and after).

        NG_012772.1(BRCA2_v001):c.68-7_316+7del is such a variant, since
        positions 68 through 316 are exactly one exon and (316-68+1)/3 = 83.
        """
        check_variant('NG_012772.1(BRCA2_v001):c.68-7_316+7del',
                      self.config, self.output)
        assert len(self.output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
        assert_equal(self.output.getOutput('removedSpliceSites'), [2])
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')
        # Todo: assert that protein products indeed have only this difference.

    def test_del_exons(self):
        """
        Deletion of two entire exons should be possible.
        """
        check_variant('NG_012772.1(BRCA2_v001):c.632-5_793+7del',
                      self.config, self.output)
        assert len(self.output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
        assert_equal(self.output.getOutput('removedSpliceSites'), [4])
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')

    def test_del_intron(self):
        """
        Deletion of an entire intron should be possible (fusion of remaining
        exonic parts).
        """
        check_variant('NG_012772.1(BRCA2_v001):c.622_674del',
                      self.config, self.output)
        assert len(self.output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
        assert_equal(self.output.getOutput('removedSpliceSites'), [2])
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')

    def test_del_intron_exact(self):
        """
        Deletion of exactly an intron should be possible (fusion of flanking
        exons).
        """
        check_variant('NG_012772.1(BRCA2_v001):c.681+1_682-1del',
                      self.config, self.output)
        assert_equal(self.output.getMessagesWithErrorCode('WOVERSPLICE'), [])
        assert_equal(self.output.getOutput('removedSpliceSites'), [2])
        # Note: The protein prediction is done, but 'newprotein' is not set
        # because we have no change. So to check if the prediction is done, we
        # check if 'oldprotein' is set and to check if the prediction is
        # correct, we check if 'newprotein' is not set.
        assert self.output.getOutput('oldprotein')
        assert not self.output.getOutput('newprotein')

    def test_del_intron_in_frame(self):
        """
        Deletion of an entire intron should be possible (fusion of remaining
        exonic parts).
        """
        check_variant('NG_012772.1(BRCA2_v001):c.622_672del',
                      self.config, self.output)
        assert len(self.output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
        assert_equal(self.output.getOutput('removedSpliceSites'), [2])
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')
        # Todo: assert that protein products indeed have only this difference.

    def test_del_exon_unknown_offsets(self):
        """
        Deletion of an entire exon with unknown offsets should be possible.
        """
        check_variant('NG_012772.1(BRCA2_v001):c.632-?_681+?del',
                      self.config, self.output)
        assert len(self.output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
        assert len(self.output.getMessagesWithErrorCode('IDELSPLICE')) > 0
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')
        # Genomic positions should be centered in flanking introns and unsure.
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'NG_012772.1:g.(17550_19725)del')
        assert 'NG_012772.1(BRCA2_v001):c.632-?_681+?del' \
               in self.output.getOutput('descriptions')
        assert 'NG_012772.1(BRCA2_i001):p.(Val211Glufs*10)' \
               in self.output.getOutput('protDescriptions')
        # Todo: .c notation should still be c.632-?_681+?del, but what about
        # other transcripts?

    def test_del_exon_unknown_offsets_in_frame(self):
        """
        Deletion of an entire exon with unknown offsets and length a
        triplicate should give a proteine product with just this deletion
        (and possibly substitutions directly before and after).

        NG_012772.1(BRCA2_v001):c.68-?_316+?del is such a variant, since
        positions 68 through 316 are exactly one exon and (316-68+1)/3 = 83.
        """
        check_variant('NG_012772.1(BRCA2_v001):c.68-?_316+?del',
                      self.config, self.output)
        assert len(self.output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
        assert len(self.output.getMessagesWithErrorCode('IDELSPLICE')) > 0
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')
        # Genomic positions should be centered in flanking introns and unsure.
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'NG_012772.1:g.(7324_11720)del')
        assert 'NG_012772.1(BRCA2_v001):c.68-?_316+?del' \
               in self.output.getOutput('descriptions')
        # Todo: .c notation should still be c.632-?_681+?del, but what about
        # other transcripts?

    def test_del_exon_unknown_offsets_composed(self):
        """
        Deletion of an entire exon with unknown offsets and another composed
        variant with exact positioning should be possible.
        """
        check_variant('NG_012772.1(BRCA2_v001):c.[632-?_681+?del;681+4del]',
                      self.config, self.output)
        assert len(self.output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
        assert len(self.output.getMessagesWithErrorCode('IDELSPLICE')) > 0
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')
        # Genomic positions should be centered in flanking introns and unsure.
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'NG_012772.1:g.[(17550_19725)del;19017del]')
        assert 'NG_012772.1(BRCA2_v001):c.[632-?_681+?del;681+4del]' \
               in self.output.getOutput('descriptions')
        # Todo: .c notation should still be c.632-?_681+?del, but what about
        # other transcripts?

    def test_del_exon_unknown_offsets_reverse(self):
        """
        Deletion of an entire exon with unknown offsets should be possible,
        also on the reverse strand.
        """
        check_variant('AL449423.14(CDKN2A_v001):c.151-?_457+?del',
                      self.config, self.output)
        assert len(self.output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
        assert len(self.output.getMessagesWithErrorCode('IDELSPLICE')) > 0
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')
        # Genomic positions should be centered in flanking introns and unsure.
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'AL449423.14:g.(60314_63683)del')
        assert 'AL449423.14(CDKN2A_v001):c.151-?_457+?del' \
               in self.output.getOutput('descriptions')
        # Todo: .c notation should still be c.632-?_681+?del, but what about
        # other transcripts?

    def test_del_exon_transcript_reference(self):
        """
        Deletion of entire exon on a transcript reference should remove the
        expected splice sites (only that of the deleted exon), and not those
        of the flanking exons (as would happen using the mechanism for genomic
        references).
        """
        check_variant('NM_018723.3:c.758_890del', self.config, self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('WOVERSPLICE')), 0)
        assert_equal(self.output.getOutput('removedSpliceSites'), [2])
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')

    def test_ins_range(self):
        """
        Insertion of a range is not implemented yet.
        """
        check_variant('AB026906.1:c.274_275ins262_268', self.config, self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('ENOTIMPLEMENTED')), 1)

    def test_delins_range(self):
        """
        Deletion/insertion of a range is not implemented yet.
        """
        check_variant('AB026906.1:c.274delins262_268', self.config, self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('ENOTIMPLEMENTED')), 1)
