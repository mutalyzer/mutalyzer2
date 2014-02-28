"""
Tests for the mapping module.
"""


#import logging; logging.basicConfig()
from nose.tools import *
from sqlalchemy import or_

from mutalyzer.db.models import Assembly
from mutalyzer.output import Output
from mutalyzer.mapping import Converter

from fixtures import database, hg19, hg19_transcript_mappings
from utils import MutalyzerTest


class TestConverter(MutalyzerTest):
    """
    Test the Converter class.
    """
    fixtures = (database, hg19, hg19_transcript_mappings)

    def setup(self):
        super(TestConverter, self).setup()
        self.output = Output(__file__)

    def _converter(self, assembly_name_or_alias):
        """
        Create a Converter instance for a given genome assembly.
        """
        assembly = Assembly.query \
            .filter(or_(Assembly.name == assembly_name_or_alias,
                        Assembly.alias == assembly_name_or_alias)) \
            .one()
        return Converter(assembly, self.output)

    def test_converter(self):
        """
        Simple test.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_003002.2:c.274G>T')
        assert_equal(genomic, 'NC_000011.9:g.111959695G>T')
        coding = converter.chrom2c(genomic, 'list')
        assert 'NM_003002.2:c.274G>T' in coding
        # Fix for r536: disable the -u and +d convention.
        #assert 'NR_028383.1:c.1-u2173C>A' in coding
        assert 'NR_028383.1:n.-2173C>A' in coding

    def test_converter_non_coding(self):
        """
        Test with variant on non-coding transcript.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NR_028383.1:n.-2173C>A')
        assert_equal(genomic, 'NC_000011.9:g.111959695G>T')
        coding = converter.chrom2c(genomic, 'list')
        assert 'NM_003002.2:c.274G>T' in coding
        # Fix for r536: disable the -u and +d convention.
        #assert 'NR_028383.1:c.1-u2173C>A' in coding
        assert 'NR_028383.1:n.-2173C>A' in coding

    def test_converter_compound(self):
        """
        Test with compound variant.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_003002.2:c.[274G>T;278A>G]')
        assert_equal(genomic, 'NC_000011.9:g.[111959695G>T;111959699A>G]')
        coding = converter.chrom2c(genomic, 'list')
        assert 'NM_003002.2:c.[274G>T;278A>G]' in coding
        assert 'NR_028383.1:n.[-2173C>A;-2177T>C]' in coding

    def test_hla_cluster(self):
        """
        Convert to primary assembly.

        Transcript NM_000500.5 is mapped to different chromosome locations,
        but we like to just see the primary assembly mapping to chromosome 6.

        See also bug #58.
        """
        # Todo: This test is bogus now that we use a fixture that has just the
        #   mapping to chromosome 6. However, I think we only get this mapping
        #   from our current source (NCBI seq_gene.md) anyway, so I'm not sure
        #   where we got the other mappings from in the past (but haven't
        #   investigated really).
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_000500.5:c.92C>T')
        assert_equal(genomic, 'NC_000006.11:g.32006291C>T')
        coding = converter.chrom2c(genomic, 'list')
        assert 'NM_000500.5:c.92C>T' in coding

    def test_converter_del_length_reverse(self):
        """
        Position converter on deletion (denoted by length) on transcripts
        located on the reverse strand.
        """
        converter = self._converter('hg19')
        coding = converter.chrom2c('NC_000022.10:g.51016285_51017117del123456789', 'list')
        # Fix for r536: disable the -u and +d convention.
        #assert 'NM_001145134.1:c.-138-u21_60del123456789' in coding
        #assert 'NR_021492.1:c.1-u5170_1-u4338del123456789' in coding
        assert 'NM_001145134.1:c.-159_60del123456789' in coding
        assert 'NR_021492.1:n.-5170_-4338del123456789' in coding

    def test_S_Venkata_Suresh_Kumar(self):
        """
        Test for correct mapping information on genes where CDS start or stop
        is exactly on the border of an exon.

        Bug reported February 24, 2012 by S Venkata Suresh Kumar.
        """
        converter = self._converter('hg19')
        coding = converter.chrom2c('NC_000001.10:g.115259837_115259837delT', 'list')
        assert 'NM_001007553.1:c.3863delA' not in coding
        assert 'NM_001007553.1:c.*953delA' in coding
        assert 'NM_001130523.1:c.*953delA' in coding

    def test_S_Venkata_Suresh_Kumar_more(self):
        """
        Another test for correct mapping information on genes where CDS start
        or stop is exactly on the border of an exon.

        Bug reported March 21, 2012 by S Venkata Suresh Kumar.
        """
        converter = self._converter('hg19')
        coding = converter.chrom2c('NC_000001.10:g.160012314_160012329del16', 'list')
        assert 'NM_002241.4:c.-27250-7_-27242del16' not in coding
        assert 'NM_002241.4:c.1-7_9del16' in coding

    def test_range_order_forward_correct(self):
        """
        Just a normal position converter call, both directions.  See Trac #95.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_003002.2:c.-1_274del')
        assert_equal(genomic, 'NC_000011.9:g.111957631_111959695del')
        coding = converter.chrom2c(genomic, 'list')
        assert 'NM_003002.2:c.-1_274del' in coding

    def test_range_order_forward_incorrect_c2chrom(self):
        """
        Incorrect order of a range on the forward strand. See Trac #95.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_003002.2:c.274_-1del')
        assert_equal(genomic, None)
        erange = self.output.getMessagesWithErrorCode('ERANGE')
        assert_equal(len(erange), 1)

    def test_range_order_reverse_correct(self):
        """
        Just a normal position converter call on the reverse strand, both
        directions. See Trac #95.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_001162505.1:c.-1_40del')
        assert_equal(genomic, 'NC_000020.10:g.48770135_48770175del')
        coding = converter.chrom2c(genomic, 'list')
        assert 'NM_001162505.1:c.-1_40del' in coding

    def test_range_order_reverse_incorrect_c2chrom(self):
        """
        Incorrect order of a range on the reverse strand. See Trac #95.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_001162505.1:c.40_-1del')
        assert_equal(genomic, None)
        erange = self.output.getMessagesWithErrorCode('ERANGE')
        assert_equal(len(erange), 1)

    def test_range_order_incorrect_chrom2c(self):
        """
        Incorrect order of a chromosomal range. See Trac #95.
        """
        converter = self._converter('hg19')
        coding = converter.chrom2c('NC_000011.9:g.111959695_111957631del', 'list')
        assert_equal(coding, None)
        erange = self.output.getMessagesWithErrorCode('ERANGE')
        assert_equal(len(erange), 1)

    def test_delins_large_ins_c2chrom(self):
        """
        Delins with multi-base insertion c. to chrom.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_003002.2:c.274delinsTAAA')
        assert_equal(genomic, 'NC_000011.9:g.111959695delinsTAAA')
        coding = converter.chrom2c(genomic, 'list')
        assert 'NM_003002.2:c.274delinsTAAA' in coding

    def test_delins_large_ins_explicit_c2chrom(self):
        """
        Delins with multi-base insertion and explicit deleted sequence c. to chrom.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_003002.2:c.274delGinsTAAA')
        assert_equal(genomic, 'NC_000011.9:g.111959695delinsTAAA')
        coding = converter.chrom2c(genomic, 'list')
        assert 'NM_003002.2:c.274delinsTAAA' in coding

    def test_delins_large_ins_chrom2c(self):
        """
        Delins with multi-base insertion chrom to c.
        """
        converter = self._converter('hg19')
        coding = converter.chrom2c('NC_000011.9:g.111959695delinsTAAA', 'list')
        assert 'NM_003002.2:c.274delinsTAAA' in coding

    def test_delins_large_ins_explicit_chrom2c(self):
        """
        Delins with multi-base insertion and explicit deleted sequence chrom to c.
        """
        converter = self._converter('hg19')
        coding = converter.chrom2c('NC_000011.9:g.111959695delGinsTAAA', 'list')
        assert 'NM_003002.2:c.274delinsTAAA' in coding

    def test_chrm_chrom2c(self):
        """
        Mitochondrial m. to c.
        """
        converter = self._converter('hg19')
        coding = converter.chrom2c('NC_012920.1:m.12030del', 'list')
        assert 'NC_012920.1(ND4_v001):c.1271del' in coding

    def test_chrm_name_chrom2c(self):
        """
        Mitochondrial m. (by chromosome name) to c.
        """
        converter = self._converter('hg19')
        variant = converter.correctChrVariant('chrM:m.12030del')
        coding = converter.chrom2c(variant, 'list')
        assert 'NC_012920.1(ND4_v001):c.1271del' in coding

    def test_chrm_c2chrom(self):
        """
        Mitochondrial c. to m.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NC_012920.1(ND4_v001):c.1271del')
        assert_equal(genomic, 'NC_012920.1:m.12030del')

    def test_nm_without_selector_chrom2c(self):
        """
        NM reference without transcript selection c. to g.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_017780.2:c.109A>T')
        assert_equal(genomic, 'NC_000008.10:g.61654100A>T')

    def test_nm_with_selector_chrom2c(self):
        """
        NM reference with transcript selection c. to g.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_017780.2(CHD7_v001):c.109A>T')
        assert_equal(genomic, 'NC_000008.10:g.61654100A>T')

    def test_nm_c2chrom_no_selector(self):
        """
        To NM reference should never result in transcript selection.
        """
        converter = self._converter('hg19')
        variant = converter.correctChrVariant('NC_000008.10:g.61654100A>T')
        coding = converter.chrom2c(variant, 'list')
        assert 'NM_017780.2:c.109A>T' in coding

    def test_incorrect_selector_c2chrom(self):
        """
        Incorrect selector.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_017780.2(CHD8):c.109A>T')
        erange = self.output.getMessagesWithErrorCode('EACCNOTINDB')
        assert_equal(len(erange), 1)

    def test_incorrect_selector_version_c2chrom(self):
        """
        Incorrect selector version.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_017780.2(CHD7_v002):c.109A>T')
        erange = self.output.getMessagesWithErrorCode('EACCNOTINDB')
        assert_equal(len(erange), 1)

    def test_no_selector_version_c2chrom(self):
        """
        Selector but no selector version.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_017780.2(CHD7):c.109A>T')
        assert_equal(genomic, 'NC_000008.10:g.61654100A>T')

    def test_incorrect_selector_no_selector_version_c2chrom(self):
        """
        Incorrect selector, no selector version.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_017780.2(CHD8):c.109A>T')
        erange = self.output.getMessagesWithErrorCode('EACCNOTINDB')
        assert_equal(len(erange), 1)

    def test_ins_seq_chrom2c(self):
        """
        Insertion of a sequence (chrom2c).
        """
        converter = self._converter('hg19')
        coding = converter.chrom2c('NC_000011.9:g.111957482_111957483insGAT', 'list')
        assert 'NM_003002.2:c.-150_-149insGAT' in coding
        assert 'NM_012459.2:c.10_11insATC' in coding

    def test_ins_seq_seq(self):
        """
        Insertion of two sequences (chrom2c).
        """
        converter = self._converter('hg19')
        coding = converter.chrom2c('NC_000011.9:g.111957482_111957483ins[GAT;AAA]', 'list')
        assert 'NM_003002.2:c.-150_-149ins[GAT;AAA]' in coding
        assert 'NM_012459.2:c.10_11ins[TTT;ATC]' in coding

    def test_ins_seq_c2chrom_reverse(self):
        """
        Insertion of a sequence on reverse strand (c2chrom).
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_012459.2:c.10_11insATC')
        assert_equal(genomic, 'NC_000011.9:g.111957482_111957483insGAT')

    def test_ins_seq_seq_c2chrom_reverse(self):
        """
        Insertion of two sequences on reverse strand (c2chrom).
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_012459.2:c.10_11ins[TTT;ATC]')
        assert_equal(genomic, 'NC_000011.9:g.111957482_111957483ins[GAT;AAA]')
