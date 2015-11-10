"""
Tests for the mutalyzer.mapping module.
"""


from __future__ import unicode_literals

import codecs
import os

import pytest

from mutalyzer.db.models import TranscriptMapping
from mutalyzer import mapping


pytestmark = pytest.mark.usefixtures('hg19_transcript_mappings')


@pytest.fixture
def converter(output, hg19):
    return mapping.Converter(hg19, output)


def test_converter(converter):
    """
    Simple test.
    """
    genomic = converter.c2chrom('NM_003002.2:c.274G>T')
    assert genomic == 'NC_000011.9:g.111959695G>T'
    coding = converter.chrom2c(genomic, 'list')
    assert 'NM_003002.2:c.274G>T' in coding
    # Fix for r536: disable the -u and +d convention.
    # assert 'NR_028383.1:c.1-u2173C>A' in coding
    assert 'NR_028383.1:n.-2173C>A' in coding


def test_converter_non_coding(converter):
    """
    Test with variant on non-coding transcript.
    """
    genomic = converter.c2chrom('NR_028383.1:n.-2173C>A')
    assert genomic == 'NC_000011.9:g.111959695G>T'
    coding = converter.chrom2c(genomic, 'list')
    assert 'NM_003002.2:c.274G>T' in coding
    # Fix for r536: disable the -u and +d convention.
    # assert 'NR_028383.1:c.1-u2173C>A' in coding
    assert 'NR_028383.1:n.-2173C>A' in coding


def test_converter_compound(converter):
    """
    Test with compound variant.
    """
    genomic = converter.c2chrom('NM_003002.2:c.[274G>T;278A>G]')
    assert genomic == 'NC_000011.9:g.[111959695G>T;111959699A>G]'
    coding = converter.chrom2c(genomic, 'list')
    assert 'NM_003002.2:c.[274G>T;278A>G]' in coding
    assert 'NR_028383.1:n.[-2173C>A;-2177T>C]' in coding


def test_hla_cluster(converter):
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
    genomic = converter.c2chrom('NM_000500.5:c.92C>T')
    assert genomic == 'NC_000006.11:g.32006291C>T'
    coding = converter.chrom2c(genomic, 'list')
    assert 'NM_000500.5:c.92C>T' in coding


def test_converter_del_length_reverse(converter):
    """
    Position converter on deletion (denoted by length) on transcripts
    located on the reverse strand.
    """
    coding = converter.chrom2c(
        'NC_000022.10:g.51016285_51017117del123456789', 'list')
    # Fix for r536: disable the -u and +d convention.
    # assert 'NM_001145134.1:c.-138-u21_60del123456789' in coding
    # assert 'NR_021492.1:c.1-u5170_1-u4338del123456789' in coding
    assert 'NM_001145134.1:c.-159_60del123456789' in coding
    assert 'NR_021492.1:n.-5170_-4338del123456789' in coding


def test_S_Venkata_Suresh_Kumar(converter):
    """
    Test for correct mapping information on genes where CDS start or stop
    is exactly on the border of an exon.

    Bug reported February 24, 2012 by S Venkata Suresh Kumar.
    """
    coding = converter.chrom2c(
        'NC_000001.10:g.115259837_115259837delT', 'list')
    assert 'NM_001007553.1:c.3863delA' not in coding
    assert 'NM_001007553.1:c.*953delA' in coding
    assert 'NM_001130523.1:c.*953delA' in coding


def test_S_Venkata_Suresh_Kumar_more(converter):
    """
    Another test for correct mapping information on genes where CDS start
    or stop is exactly on the border of an exon.

    Bug reported March 21, 2012 by S Venkata Suresh Kumar.
    """
    coding = converter.chrom2c(
        'NC_000001.10:g.160012314_160012329del16', 'list')
    assert 'NM_002241.4:c.-27250-7_-27242del16' not in coding
    assert 'NM_002241.4:c.1-7_9del16' in coding


def test_range_order_forward_correct(converter):
    """
    Just a normal position converter call, both directions.  See Trac #95.
    """
    genomic = converter.c2chrom('NM_003002.2:c.-1_274del')
    assert genomic == 'NC_000011.9:g.111957631_111959695del'
    coding = converter.chrom2c(genomic, 'list')
    assert 'NM_003002.2:c.-1_274del' in coding


def test_range_order_forward_incorrect_c2chrom(output, converter):
    """
    Incorrect order of a range on the forward strand. See Trac #95.
    """
    genomic = converter.c2chrom('NM_003002.2:c.274_-1del')
    assert genomic is None
    erange = output.getMessagesWithErrorCode('ERANGE')
    assert len(erange) == 1


def test_range_order_reverse_correct(converter):
    """
    Just a normal position converter call on the reverse strand, both
    directions. See Trac #95.
    """
    genomic = converter.c2chrom('NM_001162505.1:c.-1_40del')
    assert genomic == 'NC_000020.10:g.48770135_48770175del'
    coding = converter.chrom2c(genomic, 'list')
    assert 'NM_001162505.1:c.-1_40del' in coding


def test_range_order_reverse_incorrect_c2chrom(output, converter):
    """
    Incorrect order of a range on the reverse strand. See Trac #95.
    """
    genomic = converter.c2chrom('NM_001162505.1:c.40_-1del')
    assert genomic is None
    erange = output.getMessagesWithErrorCode('ERANGE')
    assert len(erange) == 1


def test_range_order_incorrect_chrom2c(output, converter):
    """
    Incorrect order of a chromosomal range. See Trac #95.
    """
    coding = converter.chrom2c('NC_000011.9:g.111959695_111957631del', 'list')
    assert coding is None
    erange = output.getMessagesWithErrorCode('ERANGE')
    assert len(erange) == 1


def test_delins_large_ins_c2chrom(converter):
    """
    Delins with multi-base insertion c. to chrom.
    """
    genomic = converter.c2chrom('NM_003002.2:c.274delinsTAAA')
    assert genomic == 'NC_000011.9:g.111959695delinsTAAA'
    coding = converter.chrom2c(genomic, 'list')
    assert 'NM_003002.2:c.274delinsTAAA' in coding


def test_delins_large_ins_explicit_c2chrom(converter):
    """
    Delins with multi-base insertion and explicit deleted sequence c. to chrom.
    """
    genomic = converter.c2chrom('NM_003002.2:c.274delGinsTAAA')
    assert genomic == 'NC_000011.9:g.111959695delinsTAAA'
    coding = converter.chrom2c(genomic, 'list')
    assert 'NM_003002.2:c.274delinsTAAA' in coding


def test_delins_large_ins_chrom2c(converter):
    """
    Delins with multi-base insertion chrom to c.
    """
    coding = converter.chrom2c('NC_000011.9:g.111959695delinsTAAA', 'list')
    assert 'NM_003002.2:c.274delinsTAAA' in coding


def test_delins_large_ins_explicit_chrom2c(converter):
    """
    Delins with multi-base insertion and explicit deleted sequence chrom to c.
    """
    coding = converter.chrom2c('NC_000011.9:g.111959695delGinsTAAA', 'list')
    assert 'NM_003002.2:c.274delinsTAAA' in coding


def test_chrm_chrom2c(converter):
    """
    Mitochondrial m. to c.
    """
    coding = converter.chrom2c('NC_012920.1:m.12030del', 'list')
    assert 'NC_012920.1(ND4_v001):c.1271del' in coding


def test_chrm_name_chrom2c(converter):
    """
    Mitochondrial m. (by chromosome name) to c.
    """
    variant = converter.correctChrVariant('chrM:m.12030del')
    coding = converter.chrom2c(variant, 'list')
    assert 'NC_012920.1(ND4_v001):c.1271del' in coding


def test_chrm_c2chrom(converter):
    """
    Mitochondrial c. to m.
    """
    genomic = converter.c2chrom('NC_012920.1(ND4_v001):c.1271del')
    assert genomic == 'NC_012920.1:m.12030del'


def test_nm_without_selector_chrom2c(converter):
    """
    NM reference without transcript selection c. to g.
    """
    genomic = converter.c2chrom('NM_017780.2:c.109A>T')
    assert genomic == 'NC_000008.10:g.61654100A>T'


def test_nm_with_selector_chrom2c(converter):
    """
    NM reference with transcript selection c. to g.
    """
    genomic = converter.c2chrom('NM_017780.2(CHD7_v001):c.109A>T')
    assert genomic == 'NC_000008.10:g.61654100A>T'


def test_nm_c2chrom_no_selector(converter):
    """
    To NM reference should never result in transcript selection.
    """
    variant = converter.correctChrVariant('NC_000008.10:g.61654100A>T')
    coding = converter.chrom2c(variant, 'list')
    assert 'NM_017780.2:c.109A>T' in coding


def test_incorrect_selector_c2chrom(output, converter):
    """
    Incorrect selector.
    """
    converter.c2chrom('NM_017780.2(CHD8):c.109A>T')
    erange = output.getMessagesWithErrorCode('EACCNOTINDB')
    assert len(erange) == 1


def test_incorrect_selector_version_c2chrom(output, converter):
    """
    Incorrect selector version.
    """
    converter.c2chrom('NM_017780.2(CHD7_v002):c.109A>T')
    erange = output.getMessagesWithErrorCode('EACCNOTINDB')
    assert len(erange) == 1


def test_no_selector_version_c2chrom(converter):
    """
    Selector but no selector version.
    """
    genomic = converter.c2chrom('NM_017780.2(CHD7):c.109A>T')
    assert genomic == 'NC_000008.10:g.61654100A>T'


def test_incorrect_selector_no_selector_version_c2chrom(output, converter):
    """
    Incorrect selector, no selector version.
    """
    converter.c2chrom('NM_017780.2(CHD8):c.109A>T')
    erange = output.getMessagesWithErrorCode('EACCNOTINDB')
    assert len(erange) == 1


def test_ins_seq_chrom2c(converter):
    """
    Insertion of a sequence (chrom2c).
    """
    coding = converter.chrom2c(
        'NC_000011.9:g.111957482_111957483insGAT', 'list')
    assert 'NM_003002.2:c.-150_-149insGAT' in coding
    assert 'NM_012459.2:c.10_11insATC' in coding


def test_ins_seq_seq(converter):
    """
    Insertion of two sequences (chrom2c).
    """
    coding = converter.chrom2c(
        'NC_000011.9:g.111957482_111957483ins[GAT;AAA]', 'list')
    assert 'NM_003002.2:c.-150_-149ins[GAT;AAA]' in coding
    assert 'NM_012459.2:c.10_11ins[TTT;ATC]' in coding


def test_ins_seq_c2chrom_reverse(converter):
    """
    Insertion of a sequence on reverse strand (c2chrom).
    """
    genomic = converter.c2chrom('NM_012459.2:c.10_11insATC')
    assert genomic == 'NC_000011.9:g.111957482_111957483insGAT'


def test_ins_seq_seq_c2chrom_reverse(converter):
    """
    Insertion of two sequences on reverse strand (c2chrom).
    """
    genomic = converter.c2chrom('NM_012459.2:c.10_11ins[TTT;ATC]')
    assert genomic == 'NC_000011.9:g.111957482_111957483ins[GAT;AAA]'


def test_import_mapview(hg19):
    original_count = TranscriptMapping.query.count()

    group_label = 'GRCh37.p13-Primary Assembly'

    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data',
                        'hg19.chr11.111771755-112247252.seq_gene.sorted.md')
    mapview = codecs.open(path, encoding='utf-8')
    mapview_count = sum(1 for line in mapview
                        if line.split('\t')[12] == group_label
                        and line.split('\t')[11] == 'RNA')
    mapview.seek(0)

    mapping.import_from_mapview_file(hg19, mapview, group_label)

    # Two transcripts were already in, the rest is new:
    # - NR_028383.1
    # - NM_012459.2
    assert TranscriptMapping.query.count() == original_count + mapview_count - 2

    # No changes here.
    unchanged = TranscriptMapping.query.filter_by(accession='NM_012459').one()
    assert unchanged.start == 111955524
    assert unchanged.stop == 111957522
    assert unchanged.exon_starts == [111955524, 111957364]
    assert unchanged.exon_stops == [111956186, 111957522]
    assert unchanged.cds == (111956019, 111957492)

    # We made some artificial changes to the mapview file here.
    updated = TranscriptMapping.query.filter_by(accession='NR_028383').one()
    assert updated.start == 111955524
    assert updated.stop == 111957525
    assert updated.exon_starts == [111955524, 111956700, 111957364]
    assert updated.exon_stops == [111956180, 111957034, 111957525]

    # This is a new entry.
    new = TranscriptMapping.query.filter_by(accession='NM_000317').one()
    assert new.version == 2
    assert new.start == 112097088
    assert new.stop == 112104696
    assert new.exon_starts == [112097088, 112099317, 112100931, 112101349,
                               112103886, 112104155]
    assert new.exon_stops == [112097249, 112099396, 112100953, 112101405,
                              112103956, 112104696]
    assert new.cds == (112097167, 112104278)
    assert new.gene == 'PTS'
    assert new.orientation == 'forward'
    assert new.reference_type == 'refseq'
    assert new.source == 'ncbi'
