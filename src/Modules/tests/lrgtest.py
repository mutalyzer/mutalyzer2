import sys, os, unittest, types

#make it possible to import the Modules
sys.path.append("../..")

##############################################################################
from Modules import LRGparser
from recordtest import TestRecord

class TestRecordFromLRG(TestRecord):
    def setUp(self):
        LRGparser.DEBUG = False
        f = open("lrgtest_files/LRG_1.xml")
        data = f.read()
        f.close()
        #self.record = LRGparser.createLrgRecord_new(data)
        self.record = LRGparser.createLrgRecord(data)

    def test_moltype(self): #overwrite parent
        self.assertEqual(self.record.molType, 'g')

    def test_mapping(self):
        mapp = self.record.mapping
        loc = mapp["chr_location"]
        self.assertTrue(len(loc) == 2)
        numbers = loc
        numbers.extend(mapp["lrg_location"])
        numbers.append(mapp["strand"])
        self._test_if_type(numbers, (int, long))


    def test_mapping_diffs(self):
        mapp = self.record.mapping
        for diff in mapp["diffs"]:
            numbers = [diff["lrg_start"], diff["lrg_end"],
                       diff["start"], diff["end"]]
            self._test_if_type(numbers, (int, long))
            self.assertTrue(diff["type"] in\
                    ["mismatch", "genomic_ins", "lrg_ins"])
            gs = diff.get("genomic_sequence")
            ls = diff.get("lrg_sequence")
            if gs: self._test_if_seq(gs)
            if ls: self._test_if_seq(ls)
            self.assertTrue(gs or ls)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestRecordFromLRG)
    unittest.TextTestRunner(verbosity=2).run(suite)
