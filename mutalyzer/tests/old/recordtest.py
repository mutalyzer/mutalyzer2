"""
recordtest.py contains
    TestRecord - a BaseClass for testing GenRecord.Record instances
"""
import unittest, types
from Modules import GenRecord        #test class-types

class TestRecord(unittest.TestCase):

    def setUp(self):
        # This method should be redefined by a class that is derived from it
        # and it should initialize the self.record to a GenRecord.Record
        self.record = None
        raise ImportError("This method should be redefined by the Child")

    def test_moltype(self):
        self.assertTrue(self.record.molType in ['g','c']) #Extend options

    def test_sequence(self):
        self._test_if_seq(self.record.seq)

    def test_mapping(self):
        pass #Overwrite in lrg_test

    def test_organelle(self):
        self.assertTrue(isinstance(self.record.organelle, (types.ListType,
            types.NoneType)))

    def test_source(self):
        self.assertTrue(isinstance(self.record.source, GenRecord.Gene))

    def test_description(self):
        self.assertTrue(isinstance(self.record.description, types.StringType))

    def test_genes(self):
        self.assertTrue(len(self.record.geneList)>0)
        for gene in self.record.geneList:
            self.assertTrue(isinstance(gene.name, types.StringType))
            self.assertTrue(gene.orientation in (-1, 1))
            self._test_if_loc(gene.location)
            self.assertTrue(isinstance(gene.longName, types.StringType))


    def test_transcripts(self):
        for gene in self.record.geneList:
            self.assertTrue(len(gene.transcriptList) > 0)
            for tc in gene.transcriptList:
                strings = (tc.name, tc.description, tc.proteinDescription,
                        tc.transLongName, tc.protLongName)
                self._test_if_type(strings, types.StringType)
                self._test_if_loc(tc.location)
                self._test_if_type((tc.txTable,), int)
                self.assertTrue(tc.molType in "cn")
                plists = (tc.mRNA, tc.CDS, tc.exon)
                for plist in plists:
                    self.assertTrue(isinstance(plist,
                        (types.NoneType, GenRecord.PList)))

                #self.assertTrue(any(map(isinstance, 


    def _test_if_loc(self, loc):
        self.assertTrue(len(loc) == 2)
        self._test_if_type(loc, (int, long))
        self.assertTrue(loc[0] < loc[1])

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


