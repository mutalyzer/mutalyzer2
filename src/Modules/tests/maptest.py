import sys, os, unittest, types

#make it possible to import the Modules
os.chdir("../../..")
sys.path.append("src")

import warnings
warnings.filterwarnings("ignore")

##############################################################################
from Modules import Output
from Modules import Config
from Modules import Mapper

DEBUG = False

class TestConverter(unittest.TestCase):
    def setUp(self):
        C = Config.Config()
        self.O = Output.Output(__file__, C.Output)
        self.conv = Mapper.Converter("hg19", C, self.O)

    def tearDown(self):
        if DEBUG:
            print "\nDEBUG OUTPUT"
            print self.O.getMessages()

    def test_c2chromSingle(self):
        ret = self.conv.c2chrom("NM_000051.3:c.274G>T")
        self.assertEqual(ret, "NC_000011.9:g.108099991G>T")

    def test_c2chromSingle(self):
        ret = self.conv.c2chrom("NM_000051.3:c.274_300delinsAA")
        print ret

    def test_c2chromMulti(self):
        ret = self.conv.c2chrom("NM_000051.3:c.[274G>T;273A>G]")
        self.assertEqual(ret, None)

    def test_c2chromMapping(self):
        ret = self.conv.mainMapping("NM_000051.3","c.[274G>T;273A>G]")
        self.assertEqual(ret.errorcode, 1)

    def test_chrom2cSingle(self):
        ret = self.conv.chrom2c("NC_000011.9:g.108099991G>T")
        self.assertEqual(ret, ["NM_000051.3(ATM):c.274G>T"])

    def test_chrom2cSingle2(self):
        ret = self.conv.chrom2c("NC_000011.9:g.108099991_108099999delinsAA")
        print ret


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestConverter)
    unittest.TextTestRunner(verbosity=2).run(suite)
