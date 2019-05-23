"""
Tests for the mutalyzer.parsers.genbank module.
"""


from __future__ import unicode_literals

import os

import pytest

from mutalyzer.parsers.genbank import GBparser

from fixtures import with_references


@pytest.fixture
def parser():
    return GBparser()


@pytest.mark.parametrize('products,expected', [
    (['a b c d e', 'a b C D e', 'a b c d e'], (2, 1)),
    (['a b c d e', 'a b C d e', 'a B c d e'], (1, 2)),
    (['a c d a', 'a b a', 'a a', 'a'], (1, 1)),
    ([''], (-1, -1)),
    (['', ''], (-1, -1)),
    (['a', 'a'], (-1, -1)),
    (['a', 'b'], (0, 0)),
    (['a b c', 'a b c'], (-1, -1)),
    (['a b c d a b', 'a b'], (2, 2))
])
def test_product_lists_mismatch(parser, products, expected):
    """
    Test finding mismatches in some product lists.
    """
    assert parser._find_mismatch(products) == expected


@with_references('AB026906.1')
def test_include_cds_without_mrna(settings, references, parser):
    """
    Annotated CDS without mRNA feature should be included since Mutalyzer can
    construct the RNA from the CDS.
    """
    # Contains one gene with only a CDS annotated, no mRNA.
    accession = references[0].accession
    filename = os.path.join(settings.CACHE_DIR, '%s.gb.bz2' % accession)
    record = parser.create_record(filename)
    assert record.geneList[0].transcriptList[0].name == '001'


@with_references('A1BG')
def test_only_complete_mrna_included(settings, references, parser):
    """
    Incomplete transcripts from the reference file should be ignored.
    """
    # Contains A1BG (two complete transcripts) and A1BG-AS1, ZNF497,
    # LOC100419840 (no complete transcripts).
    accession = references[0].accession
    filename = os.path.join(settings.CACHE_DIR, '%s.gb.bz2' % accession)
    record = parser.create_record(filename)
    assert [g.name for g in record.geneList] == ['A1BG']
    assert len(record.geneList[0].transcriptList) == 2


@with_references('PIK3R2')
def test_complete_and_incomplete_mrna(settings, references, parser):
    """
    Incomplete transcripts from the reference file should be ignored, but the
    gene should be included if it contains another complete transcript.
    """
    # Contains MAST3 without complete transcripts and PIK3R2 with one complete
    # and one incomplete transcript.
    accession = references[0].accession
    filename = os.path.join(settings.CACHE_DIR, '%s.gb.bz2' % accession)
    record = parser.create_record(filename)
    assert [g.name for g in record.geneList] == ['PIK3R2']
    assert len(record.geneList[0].transcriptList) == 1


@with_references('ADAC')
def test_no_version(settings, references, parser):
    """
    Genbank file without 'version' field, so BioPython record.id is the
    accession number without version. Our parser used to crash on that.

    This genbank file was contributed by Gerard Schaafsma (original
    source unknown).
    """
    accession = references[0].accession
    genbank_filename = os.path.join(settings.CACHE_DIR,
                                    '%s.gb.bz2' % accession)
    parser.create_record(genbank_filename)


@with_references('UD_144413132067')
def test_transcript_to_protein_link(settings, references, parser):
    """
    Checks if the transcript to protein linking was performed in the right manner.
    Added for `issue #430 <https://github.com/mutalyzer/mutalyzer/issues/430>`_.
    """
    accession = references[0].accession
    filename = os.path.join(settings.CACHE_DIR,'%s.gb.bz2' % accession)
    record = parser.create_record(filename)

    correct_links = {
        'CDKN1A': {
            'NM_001220777.1': 'NP_001207706.1',
            'NM_078467.2': 'NP_510867.1',
            'NM_001291549.1': 'NP_001278478.1',
            'NM_000389.4': 'NP_000380.1',
            'NM_001220778.1': 'NP_001207707.1'
            },
        'PANDAR': {
            'NR_109836.1': None
            }
        }

    real_links = {}
    for gene in record.geneList:
        current_gene_links = {}
        for transcript in gene.transcriptList:
            transcriptID = transcript.transcriptID
            proteinID = transcript.proteinID
            current_gene_links[transcriptID] = proteinID
        real_links[gene.name] = current_gene_links

    assert correct_links == real_links


@with_references('UD_150167851083')
def test_mrna_and_cds_addition_to_genes(settings, references, parser):
    """
    Checking if the mRNA and CDS features are correctly added to the genes.
    Added for `issue #468 <https://github.com/mutalyzer/mutalyzer/issues/468>`_.
    """
    accession = references[0].accession
    filename = os.path.join(settings.CACHE_DIR,'%s.gb.bz2' % accession)
    record = parser.create_record(filename)

    for gene in record.geneList:
        for transcript in gene.transcriptList:
            print transcript.transcriptID
            print transcript.proteinID

    assert [g.name for g in record.geneList] == ['UGT1A3', 'UGT1A1']
    assert [t.transcriptID for t in record.geneList[0].transcriptList] == \
           ['XR_241239.1']
    assert [t.transcriptID for t in record.geneList[1].transcriptList] == \
           ['NM_000463.2']
    assert [t.proteinID for t in record.geneList[0].transcriptList] == \
           [None]
    assert [t.proteinID for t in record.geneList[1].transcriptList] == \
           ['NP_000454.1']
