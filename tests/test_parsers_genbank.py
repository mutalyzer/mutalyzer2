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


@with_references('A1BG')
def test_only_complete_genes_included(settings, references, parser):
    """
    Incomplete genes from the reference file should be ignored.
    """
    # contains A1BG (complete) and A1BG-AS1, ZNF497, LOC100419840
    # (incomplete).
    accession = references[0].accession
    filename = os.path.join(settings.CACHE_DIR, '%s.gb.bz2' % accession)
    record = parser.create_record(filename)
    assert [g.name for g in record.geneList] == ['A1BG']

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
