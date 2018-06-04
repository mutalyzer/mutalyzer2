"""
Tests for the mutalyzer.parsers.lrg module.
"""


from __future__ import unicode_literals

import os
import bz2

from mutalyzer.parsers.lrg import create_record

from fixtures import with_references


@with_references('LRG_1')
def test_lrg_basic(settings, references):
    """
    """
    accession = references[0].accession
    filename = os.path.join(settings.CACHE_DIR, '%s.xml.bz2' % accession)
    file_handle = bz2.BZ2File(filename, 'r')
    record = create_record(file_handle.read())
    file_handle.close()

    assert [g.name for g in record.geneList] == ['COL1A1']
    assert record.geneList[0].transcriptList[0].name == '1'
    assert len(record.geneList[0].transcriptList) == 1
    assert record.geneList[0].transcriptList[0].CDS is not None
    assert record.organism == 'Homo sapiens'


@with_references('LRG_24')
def test_lrg_multiple_transcripts(settings, references):
    """
    """
    accession = references[0].accession
    filename = os.path.join(settings.CACHE_DIR, '%s.xml.bz2' % accession)
    file_handle = bz2.BZ2File(filename, 'r')
    record = create_record(file_handle.read())
    file_handle.close()

    assert len(record.geneList[0].transcriptList) == 2


@with_references('LRG_163')
def test_lrg_no_coding_sequence(settings, references):
    """
    """
    accession = references[0].accession
    filename = os.path.join(settings.CACHE_DIR, '%s.xml.bz2' % accession)
    file_handle = bz2.BZ2File(filename, 'r')
    record = create_record(file_handle.read())
    file_handle.close()

    assert len(record.geneList[0].transcriptList) == 1
    assert record.geneList[0].transcriptList[0].CDS is None
