"""
Tests for the mapping module.
"""


#import logging; logging.basicConfig()
from nose.tools import *

from mutalyzer.output import Output
from mutalyzer.mapping import Converter


class TestConverter():
    """
    Test the Converter class.
    """
    def setUp(self):
        """
        Initialize test converter module.
        """
        self.output = Output(__file__)

    def _converter(self, build):
        """
        Create a Converter instance for a given build.
        """
        return Converter(build, self.output)

    def test_converter(self):
        """
        Simple test.
        """
        converter = self._converter('hg19')
        genomic = converter.c2chrom('NM_003002.2:c.274G>T')
        assert_equal(genomic, 'NC_000011.9:g.111959695G>T')
        coding = converter.chrom2c(genomic, 'list')
        assert 'NM_003002.2:c.274G>T' in coding

    def test_hla_cluster(self):
        """
        Convert to primary assembly.

        Transcript NM_000500.5 is mapped to different chromosome locations,
        but we like to just see the primary assembly mapping to chromosome 6.

        See also bug #58.
        """
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
        assert 'NM_001145134.1:c.-138-u21_60del123456789' in coding
        assert 'NR_021492.1:c.1-u5170_1-u4338del123456789' in coding
