"""
Tests for the variantchecker module.
"""


#import logging; logging.basicConfig()
from nose.tools import *

from mutalyzer.output import Output
from mutalyzer.Db import Cache
from mutalyzer.Retriever import GenBankRetriever
from mutalyzer.variantchecker import check_variant


class TestVariantchecker():
    """
    Test the variantchecker module.
    """
    def setUp(self):
        """
        Initialize test variantchecker module.
        """
        self.output = Output(__file__)
        self.cache_database = Cache()
        self.retriever = GenBankRetriever(self.output, self.cache_database)

    def _slice(self, chromosome, start, stop, orientation):
        """
        Get a UD slice.

        Orientation: 1 for forward, 2 for reverse.
        """
        return self.retriever.retrieveslice(chromosome, start, stop, orientation)

    def _slice_gene(self, gene, organism='human', upstream=5000, downstream=2000):
        """
        Get a UD slice for a gene.
        """
        return self.retriever.retrievegene(gene, organism, upstream, downstream)

    def _load_record(self, identifier):
        """
        Load a record in the database and cache.
        """
        return self.retriever.loadrecord(identifier)

    def test_deletion_in_frame(self):
        """
        Simple in-frame deletion should give a simple description on protein
        level.
        """
        check_variant('AL449423.14(CDKN2A_v001):c.161_163del', self.output)
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
        check_variant('AL449423.14(CDKN2A_v001):c.161_162insATC', self.output)
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
                      self.output)
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
                      self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'AL449423.14:g.61938_61939delinsGGGAT')
        assert 'AL449423.14(CDKN2A_v001):c.161_162delinsATCCC' \
               in self.output.getOutput('descriptions')
        assert 'AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)' \
               in self.output.getOutput('protDescriptions')
        assert self.output.getOutput('newprotein')

    def test_est_warning_nm_est(self):
        """
        Warning for EST positioning on NM reference.
        """
        check_variant('NM_003002.2:274del', self.output)
        west = self.output.getMessagesWithErrorCode('WEST')
        assert len(west) == 1

    def test_no_est_warning_nm_c(self):
        """
        No EST warning for c. positioning on NM reference.
        """
        check_variant('NM_003002.2:c.274del', self.output)
        west = self.output.getMessagesWithErrorCode('WEST')
        assert len(west) == 0

    def test_no_est_warning_nm_n(self):
        """
        No EST warning for n. positioning on NM reference.
        """
        check_variant('NM_003002.2:n.274del', self.output)
        west = self.output.getMessagesWithErrorCode('WEST')
        assert len(west) == 0

    def test_est_warning_ng_est(self):
        """
        Warning for EST positioning on NG reference.
        """
        check_variant('NG_012772.1:128del', self.output)
        west = self.output.getMessagesWithErrorCode('WEST')
        assert len(west) == 1

    def test_no_est_warning_ng_g(self):
        """
        No EST warning for g. positioning on NG reference.
        """
        check_variant('NG_012772.1:g.128del', self.output)
        west = self.output.getMessagesWithErrorCode('WEST')
        assert len(west) == 0

    def test_no_est_warning_est_est(self):
        """
        No warning for EST positioning on EST reference.
        """
        check_variant('AA010203.1:54_55insG', self.output)
        west = self.output.getMessagesWithErrorCode('WEST')
        assert len(west) == 0

    def test_roll(self):
        """
        Just a variant where we should roll.
        """
        check_variant('NM_003002.2:c.273del', self.output)
        wroll = self.output.getMessagesWithErrorCode('WROLLFORWARD')
        assert len(wroll) > 0

    def test_no_roll(self):
        """
        Just a variant where we cannot roll.
        """
        check_variant('NM_003002.2:c.274del', self.output)
        wroll = self.output.getMessagesWithErrorCode('WROLLFORWARD')
        assert_equal(len(wroll), 0)

    def test_no_roll_splice(self):
        """
        Here we can roll but should not, because it is over a splice site.
        """
        check_variant('NM_000088.3:g.459del', self.output)
        wrollback = self.output.getMessagesWithErrorCode('IROLLBACK')
        assert len(wrollback) > 0
        wroll = self.output.getMessagesWithErrorCode('WROLLFORWARD')
        assert_equal(len(wroll), 0)

    def test_partial_roll_splice(self):
        """
        Here we can roll two positions, but should roll only one because
        otherwise it is over a splice site.
        """
        check_variant('NM_000088.3:g.494del', self.output)
        wrollback = self.output.getMessagesWithErrorCode('IROLLBACK')
        assert len(wrollback) > 0
        wroll = self.output.getMessagesWithErrorCode('WROLLFORWARD')
        assert len(wroll) > 0

    def test_roll_after_splice(self):
        """
        Here we can roll and should, we stay in the same exon.
        """
        check_variant('NM_000088.3:g.460del', self.output)
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
        check_variant('AL449423.14:g.65470_65471insTAC', self.output)
        assert 'AL449423.14(CDKN2A_v001):c.99_100insTAG' in self.output.getOutput('descriptions')
        assert_equal ('AL449423.14:g.65471_65472insACT', self.output.getIndexedOutput('genomicDescription', 0, ''))
        assert_equal(len(self.output.getMessagesWithErrorCode('WROLLFORWARD')), 1)

    def test_roll_reverse_ins(self):
        """
        Insertion that rolls on the reverse strand should not use the same
        inserted sequence in descriptions on forward and reverse strands.
        """
        check_variant('AL449423.14:g.65471_65472insACT', self.output)
        assert 'AL449423.14(CDKN2A_v001):c.99_100insTAG' in self.output.getOutput('descriptions')
        assert_equal ('AL449423.14:g.65471_65472insACT', self.output.getIndexedOutput('genomicDescription', 0, ''))
        assert_equal(len(self.output.getMessagesWithErrorCode('WROLLFORWARD')), 0)

    def test_roll_message_forward(self):
        """
        Roll warning message should only be shown for currently selected
        strand (forward).
        """
        check_variant('AL449423.14:g.65470_65471insTAC', self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('WROLLFORWARD')), 1)
        assert_equal(len(self.output.getMessagesWithErrorCode('WROLLREVERSE')), 0)

    def test_roll_message_reverse(self):
        """
        Roll warning message should only be shown for currently selected
        strand (reverse).
        """
        check_variant('AL449423.14(CDKN2A_v001):c.98_99insGTA', self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('WROLLFORWARD')), 0)
        assert_equal(len(self.output.getMessagesWithErrorCode('WROLLREVERSE')), 1)

    def test_ins_cds_start(self):
        """
        Insertion on CDS start boundary should not be included in CDS.
        """
        check_variant('NM_000143.3:c.-1_1insCAT', self.output)
        assert_equal(self.output.getIndexedOutput("newprotein", 0), None)
        # Todo: Is this a good test?

    def test_ins_cds_start_after(self):
        """
        Insertion after CDS start boundary should be included in CDS.
        """
        check_variant('NM_000143.3:c.1_2insCAT', self.output)
        assert_equal(self.output.getIndexedOutput("newprotein", 0), '?')
        # Todo: Is this a good test?

    def test_del_splice_site(self):
        """
        Deletion hitting one splice site should not do a protein prediction.
        """
        check_variant('NG_012772.1(BRCA2_v001):c.632-5_670del', self.output)
        assert len(self.output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
        assert_equal(self.output.getOutput('removedSpliceSites'), [])
        # Todo: For now, the following is how to check if no protein
        # prediction is done.
        assert not self.output.getOutput('newprotein')

    def test_del_exon(self):
        """
        Deletion of an entire exon should be possible.
        """
        check_variant('NG_012772.1(BRCA2_v001):c.632-5_681+7del', self.output)
        assert len(self.output.getMessagesWithErrorCode('WOVERSPLICE')) > 0
        assert_equal(self.output.getOutput('removedSpliceSites'), [2])
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')

    def test_del_exon_exact(self):
        """
        Deletion of exactly an exon should be possible.
        """
        check_variant('NG_012772.1(BRCA2_v001):c.632_681del', self.output)
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
        check_variant('NG_012772.1(BRCA2_v001):c.68-7_316+7del', self.output)
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
        check_variant('NG_012772.1(BRCA2_v001):c.632-5_793+7del', self.output)
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
        check_variant('NG_012772.1(BRCA2_v001):c.622_674del', self.output)
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
        check_variant('NG_012772.1(BRCA2_v001):c.681+1_682-1del', self.output)
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
        check_variant('NG_012772.1(BRCA2_v001):c.622_672del', self.output)
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
        check_variant('NG_012772.1(BRCA2_v001):c.632-?_681+?del', self.output)
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
        check_variant('NG_012772.1(BRCA2_v001):c.68-?_316+?del', self.output)
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
                      self.output)
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
                      self.output)
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
        check_variant('NM_018723.3:c.758_890del', self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('WOVERSPLICE')), 0)
        assert_equal(self.output.getOutput('removedSpliceSites'), [2])
        # Todo: For now, the following is how to check if protein
        # prediction is done.
        assert self.output.getOutput('newprotein')

    def test_ins_range(self):
        """
        Insertion of a range is not implemented yet.
        """
        check_variant('AB026906.1:c.274_275ins262_268', self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('ENOTIMPLEMENTED')), 1)

    def test_delins_range(self):
        """
        Deletion/insertion of a range is not implemented yet.
        """
        check_variant('AB026906.1:c.274delins262_268', self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('ENOTIMPLEMENTED')), 1)

    def test_contig_reference(self):
        """
        Variant description on a CONTIG RefSeq reference.
        """
        check_variant('NG_005990.1:g.1del', self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'NG_005990.1:g.1del')

    def test_no_reference(self):
        """
        Variant description without a reference.
        """
        check_variant('g.244355733del', self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('ENOREF')), 1)

    def test_chromosomal_positions(self):
        """
        Variants on transcripts in c. notation should have chromosomal positions
        defined.
        """
        check_variant('NM_003002.2:c.274G>T', self.output)
        assert_equal(self.output.getIndexedOutput('rawVariantsChromosomal', 0),
                     ('chr11', '+', [('274G>T', (111959695, 111959695))]))

    def test_ex_notation(self):
        """
        Variant description using EX notation should not crash but deletion of
        one exon should delete two splice sites.
        """
        check_variant('NM_002001.2:c.EX1del', self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('IDELSPLICE')), 1)

    def test_lrg_reference(self):
        """
        We should be able to use LRG reference sequence without error.
        """
        check_variant('LRG_1t1:c.266G>T', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'LRG_1:g.6855G>T')

    def test_lrg_reference_new(self):
        """
        We should be able to use new LRG reference sequence without error.

        Note that all LRG sequences are now in a new format and essentially
        this test is no different from the previous, except that LRG_218 was
        not yet in our cache which makes it easier to test the new format.
        """
        check_variant('LRG_218:c.1786_1788delAAT', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)

    def test_non_numeric_locus_tag_ending(self):
        """
        Locus tag in NC_002128 does not end in an underscore and three digits
        but we should not crash on it.
        """
        check_variant('NC_002128(tagA):c.3del', self.output)

    def test_gi_reference_plain(self):
        """
        Test reference sequence notation with GI number.
        """
        assert self._load_record('NM_002001.2')  # Make sure it's in our database
        check_variant('31317229:c.6del', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     '31317229:n.105del')
        assert '31317229(FCER1A_v001):c.6del' \
               in self.output.getOutput('descriptions')

    def test_gi_reference_prefix(self):
        """
        Test reference sequence notation with GI number and prefix.
        """
        assert self._load_record('NM_002001.2')  # Make sure it's in our database
        check_variant('GI31317229:c.6del', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     '31317229:n.105del')
        assert '31317229(FCER1A_v001):c.6del' \
               in self.output.getOutput('descriptions')

    def test_gi_reference_prefix_colon(self):
        """
        Test reference sequence notation with GI number and prefix with colon.
        """
        assert self._load_record('NM_002001.2')  # Make sure it's in our database
        check_variant('GI:31317229:c.6del', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     '31317229:n.105del')
        assert '31317229(FCER1A_v001):c.6del' \
               in self.output.getOutput('descriptions')

    def test_nop_nm(self):
        """
        Variant on NM without effect should be described as '='.
        """
        check_variant('NM_002001.2:c.1_3delinsATG', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'NM_002001.2:n.=')
        assert 'NM_002001.2(FCER1A_v001):c.=' \
               in self.output.getOutput('descriptions')

    def test_nop_ud(self):
        """
        Variant on UD without effect should be described as '='.
        """
        ud = self._slice_gene('DMD')
        check_variant(ud + ':g.5T>T', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicChromDescription', 0),
                     'NC_000023.10:g.=')
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     ud + ':g.=')
        assert ud + '(DMD_v001):c.=' \
               in self.output.getOutput('descriptions')

    def test_ud_reverse_sequence(self):
        """
        Variant on UD from reverse strand should have reverse complement
        sequence.
        """
        ud = self._slice_gene('DPYD')
        check_variant(ud + '(DPYD_v1):c.85C>T', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicChromDescription', 0),
                     'NC_000001.10:g.98348885G>A')
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     ud + ':g.42731C>T')
        assert ud + '(DPYD_v001):c.85C>T' \
               in self.output.getOutput('descriptions')

    def test_ud_forward_sequence(self):
        """
        Variant on UD from forward strand should have forward sequence.
        """
        ud = self._slice_gene('MARK1')
        check_variant(ud + '(MARK1_v001):c.400T>C', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicChromDescription', 0),
                     'NC_000001.10:g.220773181T>C')
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     ud + ':g.76614T>C')
        assert ud + '(MARK1_v001):c.400T>C' \
               in self.output.getOutput('descriptions')

    def test_ud_reverse_range(self):
        """
        Variant on UD from reverse strand should have reversed range
        positions.
        """
        ud = self._slice('NC_000009.11', 32922603, 33006639, 2)
        check_variant(ud + ':g.10624_78132del', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicChromDescription', 0),
                     'NC_000009.11:g.32928508_32996016del')
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     ud + ':g.10624_78132del')

    def test_ud_forward_range(self):
        """
        Variant on UD from forward strand should have forward range positions.
        """
        ud = self._slice_gene('MARK1')
        check_variant(ud + '(MARK1_v001):c.400_415del', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicChromDescription', 0),
                     'NC_000001.10:g.220773181_220773196del')
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     ud + ':g.76614_76629del')

    def test_ud_reverse_del_length(self):
        """
        Variant on UD from reverse strand should have reversed range
        positions, but not reverse complement of first argument (it is not a
        sequence, but a length).
        """
        ud = self._slice('NC_000009.11', 32922603, 33006639, 2)
        check_variant(ud + ':g.10624_78132del67509', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicChromDescription', 0),
                     'NC_000009.11:g.32928508_32996016del')
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     ud + ':g.10624_78132del')

    def test_ud_reverse_roll(self):
        """
        Variant on UD from reverse strand should roll the oposite direction.

        The situation is as follows:

                      G    A    A    A    T    T
                c.   102  103  104  105  106  107
                g.   748  749  750  751  752  753
            chr g.   868  867  866  865  864  863
        """
        ud = self._slice_gene('DPYD')
        check_variant(ud + '(DPYD_v001):c.104del', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicChromDescription', 0),
                     'NC_000001.10:g.98348867del')
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     ud + ':g.42751del')
        assert ud + '(DPYD_v001):c.105del' \
               in self.output.getOutput('descriptions')

    def test_ud_forward_roll(self):
        """
        Variant on UD from forward strand should roll the same.

        The situation is as follows:

                      A    T    T    T    A
                c.   398  399  400  401  402
                g.   612  613  614  615  616
            chr g.   179  180  181  182  183
        """
        ud = self._slice_gene('MARK1')
        check_variant(ud + '(MARK1_v001):c.400del', self.output)
        error_count, _, _ = self.output.Summary()
        assert_equal(error_count, 0)
        assert_equal(self.output.getIndexedOutput('genomicChromDescription', 0),
                     'NC_000001.10:g.220773182del')
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     ud + ':g.76615del')
        assert ud + '(MARK1_v001):c.401del' \
               in self.output.getOutput('descriptions')

    def test_deletion_with_sequence_forward_genomic(self):
        """
        Specify the deleted sequence in a deletion.
        """
        check_variant('AL449423.14:g.65471_65472delTC', self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'AL449423.14:g.65471_65472del')
        assert 'AL449423.14(CDKN2A_v001):c.98_99del' \
               in self.output.getOutput('descriptions')

    def test_deletion_with_length_forward_genomic(self):
        """
        Specify the deleted sequence length in a deletion.
        """
        check_variant('AL449423.14:g.65471_65472del2', self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'AL449423.14:g.65471_65472del')
        assert 'AL449423.14(CDKN2A_v001):c.98_99del' \
               in self.output.getOutput('descriptions')

    def test_deletion_with_sequence_reverse_coding(self):
        """
        Specify the deleted sequence in a deletion on the reverse strand.
        """
        check_variant('AL449423.14(CDKN2A_v001):c.161_163delTGG', self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'AL449423.14:g.61937_61939del')
        assert 'AL449423.14(CDKN2A_v001):c.161_163del' \
               in self.output.getOutput('descriptions')

    def test_deletion_with_length_reverse_coding(self):
        """
        Specify the deleted sequence length in a deletion on the reverse strand.
        """
        check_variant('AL449423.14(CDKN2A_v001):c.161_163del3', self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'AL449423.14:g.61937_61939del')
        assert 'AL449423.14(CDKN2A_v001):c.161_163del' \
               in self.output.getOutput('descriptions')

    def test_deletion_with_sequence_reverse_ng_coding(self):
        """
        Specify the deleted sequence in a deletion on the reverse strand
        using a genomic reference.
        """
        check_variant('NG_008939.1:c.155_157delAAC', self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'NG_008939.1:g.5206_5208del')
        assert 'NG_008939.1(PCCB_v001):c.155_157del' \
               in self.output.getOutput('descriptions')

    def test_deletion_with_length_reverse_ng_coding(self):
        """
        Specify the deleted sequence length in a deletion on the reverse strand
        using a genomic reference.
        """
        check_variant('NG_008939.1:c.155_157del3', self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'NG_008939.1:g.5206_5208del')
        assert 'NG_008939.1(PCCB_v001):c.155_157del' \
               in self.output.getOutput('descriptions')

    def test_inversion(self):
        """
        Inversion variant.
        """
        check_variant('AB026906.1:c.274_275inv', self.output)
        assert_equal(self.output.getIndexedOutput('genomicDescription', 0),
                     'AB026906.1:g.7872_7873inv')
        assert 'AB026906.1(SDHD_v001):c.274_275inv' \
            in self.output.getOutput('descriptions')

    def test_delins_with_length(self):
        """
        Delins with explicit length of deleted sequence (bug #108).
        """
        check_variant('NM_000193.2:c.108_109del2insG', self.output)
        assert 'NM_000193.2(SHH_i001):p.(Lys38Serfs*2)' in self.output.getOutput('protDescriptions')

    def test_protein_level_description(self):
        """
        Currently protein level descriptions are not implemented.
        """
        check_variant('NG_009105.1(OPN1LW):p.=', self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('ENOTIMPLEMENTED')), 1)

    def test_protein_reference(self):
        """
        Currently protein references are not implemented.
        """
        check_variant('NP_064445.1:p.=', self.output)
        assert_equal(len(self.output.getMessagesWithErrorCode('ENOTIMPLEMENTED')), 1)

    def test_wnomrna_other(self):
        """
        Warning for no mRNA field on other than currently selected transcript
        should give WNOMRNA_OTHER warning.
        """
        ud = self._slice_gene('A1BG') # Contains ZNF497 (v1 and v2) with no mRNA
        check_variant(ud + '(A1BG_v001):c.13del', self.output)
        wnomrna_other = self.output.getMessagesWithErrorCode('WNOMRNA_OTHER')
        assert len(wnomrna_other) == 2

    def test_wnomrna(self):
        """
        Warning for no mRNA field on currently selected transcript should give
        WNOMRNA warning.
        """
        ud = self._slice_gene('A1BG') # Contains ZNF497 (v1 and v2) with no mRNA
        check_variant(ud + '(ZNF497_v001):c.13del', self.output)
        wnomrna = self.output.getMessagesWithErrorCode('WNOMRNA')
        wnomrna_other = self.output.getMessagesWithErrorCode('WNOMRNA_OTHER')
        assert len(wnomrna) == len(wnomrna_other) == 1

    def test_mrna_ref_adjacent_exons_warn(self):
        """
        Warning for mRNA reference where exons are not adjacent.

        In L41870.1 exon 15 ends on 1558 and 16 starts on 1636.
        """
        check_variant('L41870.1:c.1del', self.output)
        w_exon_annotation = self.output.getMessagesWithErrorCode('WEXON_ANNOTATION')
        assert len(w_exon_annotation) == 1

    def test_mrna_ref_adjacent_exons_no_warn(self):
        """
        No warning for mRNA reference where exons are adjacent.
        """
        check_variant('NM_133378.3:c.1del', self.output)
        w_exon_annotation = self.output.getMessagesWithErrorCode('WEXON_ANNOTATION')
        assert len(w_exon_annotation) == 0
