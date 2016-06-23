"""
Tests for the mutalyzer.Scheduler module.
"""


from __future__ import unicode_literals

import bz2
import os
import io

import pytest
from Bio import Entrez
from mock import patch

from mutalyzer.config import settings
from mutalyzer.db.models import BatchJob
from mutalyzer import File
from mutalyzer import output
from mutalyzer import Scheduler

from fixtures import with_references


pytestmark = pytest.mark.usefixtures('db')


def _batch_job(batch_file, expected, job_type, argument=None):
    file_instance = File.File(output.Output('test'))
    scheduler = Scheduler.Scheduler()

    job, columns = file_instance.parseBatchFile(batch_file)
    result_id = scheduler.addJob('test@test.test', job, columns,
                                 job_type, argument=argument)

    batch_job = BatchJob.query.filter_by(result_id=result_id).one()

    left = batch_job.batch_queue_items.count()
    assert left == len(expected)

    scheduler.process()

    left = batch_job.batch_queue_items.count()
    assert left == 0

    filename = 'batch-job-%s.txt' % result_id
    result = io.open(os.path.join(settings.CACHE_DIR, filename),
                     encoding='utf-8')

    next(result)  # Header.
    assert expected == [line.strip().split('\t') for line in result]


def _batch_job_plain_text(variants, expected, job_type, argument=None):
    batch_file = io.BytesIO(('\n'.join(variants) + '\n').encode('utf-8'))
    _batch_job(batch_file, expected, job_type, argument=argument)


def test_syntax_checker():
    """
    Simple syntax checker batch job.
    """
    variants = ['AB026906.1:c.274G>T',
                'AL449423.14(CDKN2A_v002):c.5_400del']
    expected = [['AB026906.1:c.274G>T',
                 'OK'],
                ['AL449423.14(CDKN2A_v002):c.5_400del',
                 'OK']]
    _batch_job_plain_text(variants, expected, 'syntax-checker')


def test_snp_converter():
    """
    Simple SNP converter batch job.
    """
    snps = ['rs9919552']
    expected = [['rs9919552',
                 '|'.join(['NC_000011.9:g.111959625C>T',
                           'NG_012337.2:g.7055C>T',
                           'NG_033145.1:g.2898G>A',
                           'NM_001276503.1:c.169+928C>T',
                           'NM_001276504.1:c.87C>T',
                           'NM_001276506.1:c.204C>T',
                           'NM_003002.3:c.204C>T',
                           'NP_001263433.1:p.Ser29=',
                           'NP_001263435.1:p.Ser68=',
                           'NP_002993.1:p.Ser68=',
                           'NR_077060.1:n.288C>T',
                           'NT_033899.8:g.15522041C>T'])]]

    # Patch Bio.Entrez.efetch to return dbSNP record for rs9919552.
    def mock_efetch(*args, **kwargs):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data',
                            'rs9919552.xml.bz2')
        return bz2.BZ2File(path)

    with patch.object(Entrez, 'efetch', mock_efetch):
        _batch_job_plain_text(snps, expected, 'snp-converter')


def test_large_input():
    """
    Simple batch job with large input.
    """
    variants = ['chr13:g.114503915delCACCTGCGGGAGGTGAGGGGCGCTGGGGACCCCCG'
                'TATCTACACCTGCGGGAGGTGAGGGGCGCTGGGGACCCCTATATCTACACCTGAG'
                'GGAGGTGinsTGCCTGCGGGAGGTGAGGGGCGCTGGGGACCCCCGTATCTACACC'
                'TGCGGGAGGTGAGGGGCGCTGGGGACCCCTATATCTACACCTGAGGGAGGTG']

    expected = [['InputFields: chr13:g.114503915delCACCTGCGGGAGGTGAGGGGC'
                 'GCTGGGGACCCCCGTATCTACACCTGCGGGAGGTGAGGGGCGCTGGGGACCCCT'
                 'ATATCTACACCTGAGGGAGGTGinsTGCCTGCGGGAGGTGAGGGGCGCTGGGGA'
                 'CCCCCGTATCTACACCTGCGGGAGGTGAGGG...',
                 '(Scheduler): Entry could not be formatted correctly, '
                 'check batch input file help for details']]
    _batch_job_plain_text(variants, expected, 'syntax-checker')


@with_references('AB026906.1', 'NM_000059.3')
def test_name_checker():
    """
    Simple name checker batch job.
    """
    variants = ['AB026906.1:c.274G>T',
                'NM_000059.3:c.670G>T']
    expected = [['AB026906.1:c.274G>T',
                 '(GenRecord): No mRNA field found for gene SDHD, '
                 'transcript variant 001 in record, constructing it from '
                 'CDS. Please note that descriptions exceeding CDS '
                 'boundaries are invalid.',
                 'AB026906.1',
                 'SDHD_v001',
                 'c.274G>T',
                 'g.7872G>T',
                 'c.274G>T',
                 'p.(Asp92Tyr)',
                 'SDHD_v001:c.274G>T',
                 'SDHD_v001:p.(Asp92Tyr)',
                 '',
                 '',
                 'BAA81889.1',
                 'AB026906.1(SDHD_v001):c.274G>T',
                 'AB026906.1(SDHD_i001):p.(Asp92Tyr)',
                 'CviQI,RsaI',
                 'BccI'],
                ['NM_000059.3:c.670G>T',
                 '',
                 'NM_000059.3',
                 'BRCA2_v001',
                 'c.670G>T',
                 'n.897G>T',
                 'c.670G>T',
                 'p.(Asp224Tyr)',
                 'BRCA2_v001:c.670G>T',
                 'BRCA2_v001:p.(Asp224Tyr)',
                 '',
                 'NM_000059.3',
                 'NP_000050.2',
                 'NM_000059.3(BRCA2_v001):c.670G>T',
                 'NM_000059.3(BRCA2_i001):p.(Asp224Tyr)',
                 '',
                 'BspHI,CviAII,FatI,Hpy188III,NlaIII']]
    _batch_job_plain_text(variants, expected, 'name-checker')


def test_name_checker_altered():
    """
    Name checker job with altered entries.
    """
    variants = ['NM_000059:c.670dup',
                'NM_000059:c.670G>T',
                'NM_000059.3:c.670G>T']
    expected = [['NM_000059:c.670dup',
                 '|'.join(['(Retriever): No version number is given, '
                           'using NM_000059.3. Please use this number to '
                           'reduce downloading overhead.',
                           '(Scheduler): All further occurrences of '
                           'NM_000059 will be substituted by '
                           'NM_000059.3']),
                 'NM_000059',
                 'BRCA2_v001',
                 'c.670dup',
                 'n.897dup',
                 'c.670dup',
                 'p.(Asp224Glyfs*5)',
                 'BRCA2_v001:c.670dup',
                 'BRCA2_v001:p.(Asp224Glyfs*5)',
                 '',
                 'NM_000059.3',
                 'NP_000050.2',
                 'NM_000059(BRCA2_v001):c.670dup',
                 'NM_000059(BRCA2_i001):p.(Asp224Glyfs*5)',
                 'BciVI',
                 'BspHI,Hpy188III'],
                ['NM_000059.3:c.670G>T',
                 '(Scheduler): Entry altered before execution',
                 'NM_000059.3',
                 'BRCA2_v001',
                 'c.670G>T',
                 'n.897G>T',
                 'c.670G>T',
                 'p.(Asp224Tyr)',
                 'BRCA2_v001:c.670G>T',
                 'BRCA2_v001:p.(Asp224Tyr)',
                 '',
                 'NM_000059.3',
                 'NP_000050.2',
                 'NM_000059.3(BRCA2_v001):c.670G>T',
                 'NM_000059.3(BRCA2_i001):p.(Asp224Tyr)',
                 '',
                 'BspHI,CviAII,FatI,Hpy188III,NlaIII'],
                ['NM_000059.3:c.670G>T',
                 '',
                 'NM_000059.3',
                 'BRCA2_v001',
                 'c.670G>T',
                 'n.897G>T',
                 'c.670G>T',
                 'p.(Asp224Tyr)',
                 'BRCA2_v001:c.670G>T',
                 'BRCA2_v001:p.(Asp224Tyr)',
                 '',
                 'NM_000059.3',
                 'NP_000050.2',
                 'NM_000059.3(BRCA2_v001):c.670G>T',
                 'NM_000059.3(BRCA2_i001):p.(Asp224Tyr)',
                 '',
                 'BspHI,CviAII,FatI,Hpy188III,NlaIII']]

    # Patch GenBankRetriever.fetch to return the contents of NM_000059.3
    # for NM_000059.
    def mock_efetch(*args, **kwargs):
        if kwargs.get('id') != 'NM_000059':
            return Entrez.efetch(*args, **kwargs)
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data',
                            'NM_000059.3.gb.bz2')
        return bz2.BZ2File(path)

    with patch.object(Entrez, 'efetch', mock_efetch):
        _batch_job_plain_text(variants, expected, 'name-checker')


@with_references('NM_000059.3')
def test_name_checker_skipped():
    """
    Name checker job with skipped entries.
    """
    variants = ['NM_1234567890.3:c.670G>T',
                'NM_1234567890.3:c.570G>T',
                'NM_000059.3:c.670G>T']
    expected = [['NM_1234567890.3:c.670G>T',
                 '(Retriever): Could not retrieve NM_1234567890.3.|'
                 '(Scheduler): All further occurrences with '
                 '\'NM_1234567890.3\' will be skipped'],
                ['NM_1234567890.3:c.570G>T',
                 '(Scheduler): Skipping entry'],
                ['NM_000059.3:c.670G>T',
                 '',
                 'NM_000059.3',
                 'BRCA2_v001',
                 'c.670G>T',
                 'n.897G>T',
                 'c.670G>T',
                 'p.(Asp224Tyr)',
                 'BRCA2_v001:c.670G>T',
                 'BRCA2_v001:p.(Asp224Tyr)',
                 '',
                 'NM_000059.3',
                 'NP_000050.2',
                 'NM_000059.3(BRCA2_v001):c.670G>T',
                 'NM_000059.3(BRCA2_i001):p.(Asp224Tyr)',
                 '',
                 'BspHI,CviAII,FatI,Hpy188III,NlaIII']]

    # Patch GenBankRetriever.fetch to fail on NM_1234567890.3.
    def mock_efetch(*args, **kwargs):
        if kwargs.get('id') != 'NM_1234567890.3':
            return Entrez.efetch(*args, **kwargs)
        raise IOError()

    with patch.object(Entrez, 'efetch', mock_efetch):
        _batch_job_plain_text(variants, expected, 'name-checker')


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_position_converter():
    """
    Simple position converter batch job.
    """
    variants = ['chr11:g.111959695G>T']
    expected = [['chr11:g.111959695G>T',
                 '',
                 'NC_000011.9:g.111959695G>T',
                 'NM_003002.2:c.274G>T',
                 'NM_012459.2:c.-2203C>A',
                 'NR_028383.1:n.-2173C>A']]
    _batch_job_plain_text(variants, expected, 'position-converter', 'hg19')


def test_ods_file():
    """
    OpenDocument Spreadsheet input for batch job.
    """
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'data',
                        'batch_input.ods')
    batch_file = open(path, 'rb')
    expected = [['AB026906.1:c.274G>T',
                 'OK'],
                ['AL449423.14(CDKN2A_v002):c.5_400del',
                 'OK']]

    _batch_job(batch_file, expected, 'syntax-checker')


def test_sxc_file():
    """
    OpenOffice.org 1.x Calc spreadsheet input for batch job.
    """
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'data',
                        'batch_input.sxc')
    batch_file = open(path, 'rb')
    expected = [['AB026906.1:c.274G>T',
                 'OK'],
                ['AL449423.14(CDKN2A_v002):c.5_400del',
                 'OK']]

    _batch_job(batch_file, expected, 'syntax-checker')


def test_xls_file():
    """
    Microsoft Excel 97/2000/XP/2003 input for batch job.
    """
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'data',
                        'batch_input.xls')
    batch_file = open(path, 'rb')
    expected = [['AB026906.1:c.274G>T',
                 'OK'],
                ['AL449423.14(CDKN2A_v002):c.5_400del',
                 'OK']]

    _batch_job(batch_file, expected, 'syntax-checker')


def test_xlsx_file():
    """
    Office Open XML Spreadsheet input for batch job.
    """
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'data',
                        'batch_input.xlsx')
    batch_file = open(path, 'rb')
    expected = [['AB026906.1:c.274G>T',
                 'OK'],
                ['AL449423.14(CDKN2A_v002):c.5_400del',
                 'OK']]

    _batch_job(batch_file, expected, 'syntax-checker')


def test_invalid_zip_file():
    """
    Random zip file input for batch job (invalid).
    """
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'data',
                        'image.zip')
    batch_file = open(path, 'rb')

    file_instance = File.File(output.Output('test'))
    job, columns = file_instance.parseBatchFile(batch_file)
    assert job is None


def test_unicode_input():
    """
    Simple input with some non-ASCII unicode characters.
    """
    variants = ['\u2026AB026906.1:c.274G>T',
                '\u2026AL449423.14(CDKN2A_v002):c.5_400del']
    expected = [['\u2026AB026906.1:c.274G>T',
                 '(grammar): Expected W:(0123...) (at char 0), (line:1, col:1)'],
                ['\u2026AL449423.14(CDKN2A_v002):c.5_400del',
                 '(grammar): Expected W:(0123...) (at char 0), (line:1, col:1)']]
    _batch_job_plain_text(variants, expected, 'syntax-checker')


def test_windows_1252_input():
    """
    Simple input encoded as WINDOWS-1252.
    """
    variants = ['AB026906.1:c.274G>T',
                # Encoded as WINDOWS-1252, the following is not valid UTF8.
                'NM_000052.4:c.2407\u20132A>G',
                'AL449423.14(CDKN2A_v002):c.5_400del']
    batch_file = io.BytesIO(('\n'.join(variants) + '\n').encode('WINDOWS-1252'))
    expected = [['AB026906.1:c.274G>T',
                 'OK'],
                ['NM_000052.4:c.2407\u20132A>G',
                 '(grammar): Expected W:(acgt...) (at char 18), (line:1, col:19)'],
                ['AL449423.14(CDKN2A_v002):c.5_400del',
                 'OK']]

    _batch_job(batch_file, expected, 'syntax-checker')
