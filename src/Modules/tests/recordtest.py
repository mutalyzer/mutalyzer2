"""
recordtest.py contains
    TestRecord - a BaseClass for testing GenRecord.Record instances
"""
import unittest, types

class TestRecord(unittest.TestCase):

    def setUp(self):
        # This method should be redefined by a class that is derived from it
        # and it should initialize the self.record to a GenRecord.Record
        self.record = None
        raise ImportError("This method should be redefined by the Child")

    def test_moltype(self):
        self.assertEqual(self.record.molType, 'g')

    def test_sequence(self):
        self._test_if_seq(self.record.seq)

    def test_genes(self):
        self.assertTrue(len(self.record.geneList)>0)
        for gene in self.record.geneList:
            pass


    def test_transcripts(self):
        pass

    def _test_if_seq(self, seq):
        self.assertTrue(len(seq)>0)
        self.assertTrue(set(seq) - set('ACTG') == set())

    def _test_if_type(self, arr, typs):
        if type(typs) == types.TypeType:
            typs = (typs,)
        for el in arr:
            self.assertTrue(isinstance(el, typs))


if __name__ == "__main__":
    # This file should be imported 
    pass


