"""
Tests for the Scheduler module.
"""


import bz2
import os
import StringIO

#import logging; logging.basicConfig()
from Bio import Entrez
from mock import patch
from nose.tools import *

from mutalyzer.config import settings
from mutalyzer.db.models import BatchJob, BatchQueueItem
from mutalyzer import File
from mutalyzer import output
from mutalyzer import Scheduler

from fixtures import database, cache
from utils import MutalyzerTest
from utils import fix


class TestScheduler(MutalyzerTest):
    """
    Test the Scheduler class.
    """
    fixtures = (database, )

    @staticmethod
    def _batch_job(variants, expected, job_type, argument=None):
        file_instance = File.File(output.Output('test'))
        scheduler = Scheduler.Scheduler()

        batch_file = StringIO.StringIO('\n'.join(variants) + '\n')
        job, columns = file_instance.parseBatchFile(batch_file)
        result_id = scheduler.addJob('test@test.test', job, columns,
                                     job_type)

        batch_job = BatchJob.query.filter_by(result_id=result_id).one()

        left = batch_job.batch_queue_items.count()
        assert_equal(left, len(variants))

        scheduler.process()

        left = batch_job.batch_queue_items.count()
        assert_equal(left, 0)

        filename = 'batch-job-%s.txt' % result_id
        result = open(os.path.join(settings.CACHE_DIR, filename))

        next(result) # Header.
        assert_equal(expected,
                     [line.strip().split('\t') for line in result])

    def test_syntax_checker(self):
        """
        Simple syntax checker batch job.
        """
        variants = ['AB026906.1:c.274G>T',
                    'AL449423.14(CDKN2A_v002):c.5_400del']
        expected = [['AB026906.1:c.274G>T',
                     'OK'],
                    ['AL449423.14(CDKN2A_v002):c.5_400del',
                     'OK']]
        self._batch_job(variants, expected, 'syntax-checker')

    @fix(cache('AB026906.1', 'NM_000059.3'))
    def test_name_checker(self):
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
        self._batch_job(variants, expected, 'name-checker')

    def test_name_checker_altered(self):
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
            self._batch_job(variants, expected, 'name-checker')

    @fix(cache('NM_000059.3'))
    def test_name_checker_skipped(self):
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
            self._batch_job(variants, expected, 'name-checker')
