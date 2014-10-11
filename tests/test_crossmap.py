"""
Tests for the Crossmap module.
"""


from __future__ import unicode_literals

#import logging; logging.basicConfig()

from mutalyzer.Crossmap import Crossmap

from utils import MutalyzerTest


class TestCrossmap(MutalyzerTest):
    """
    Test the Crossmap class.
    """
    def test_splice_sites(self):
        """
        Check whether the gene on the forward strand has the right splice
        sites in c. notation.
        """
        rna = [5002, 5125, 27745, 27939, 58661, 58762, 74680, 74767, 103409,
               103528, 119465, 119537, 144687, 144810, 148418, 149215]
        cds = [27925, 74736]
        cm = Crossmap(rna, cds, 1)
        assert (cm._Crossmap__crossmapping ==
                [-304, -181, -180, 15, 16, 117, 118, 205, 206, 325, 326, 398,
                 399, 522, 523, 1320])

    def test_splice_sites_reverse(self):
        """
        Check whether the gene on the reverse strand has the right splice
        sites in c. notation.
        """
        rna = [2000, 2797, 6405, 6528, 31678, 31750, 47687, 47806, 76448,
               76535, 92453, 92554, 123276, 123470, 146090, 146213]
        cds = [76479, 123290]
        cm = Crossmap(rna, cds, -1)
        assert (cm._Crossmap__crossmapping ==
                [1320, 523, 522, 399, 398, 326, 325, 206, 205, 118, 117, 16,
                 15, -180, -181, -304])

    def test_g2x(self):
        """
        Do some g. to c. conversion checking for the gene on the forward
        strand.
        """
        rna = [5002, 5125, 27745, 27939, 58661, 58762, 74680, 74767, 103409,
               103528, 119465, 119537, 144687, 144810, 148418, 149215]
        cds = [27925, 74736]
        cm = Crossmap(rna, cds, 1)
        # Fix for r536: disable the -u and +d convention.
        #assert cm.tuple2string(cm.g2x(5001)) == '-304-u1'
        assert cm.tuple2string(cm.g2x(5001)) == '-305'
        assert cm.tuple2string(cm.g2x(5124)) == '-182'
        assert cm.tuple2string(cm.g2x(5126)) == '-181+1'
        assert cm.tuple2string(cm.g2x(27924)) == '-1'
        assert cm.tuple2string(cm.g2x(27925)) == '1'
        assert cm.tuple2string(cm.g2x(58660)) == '16-1'
        assert cm.tuple2string(cm.g2x(74736)) == '174'
        assert cm.tuple2string(cm.g2x(74737)) == '*1'
        assert cm.tuple2string(cm.g2x(103408)) == '*32-1'
        assert cm.tuple2string(cm.g2x(103410)) == '*33'
        # Fix for r536: disable the -u and +d convention.
        #assert cm.tuple2string(cm.g2x(149216)) == '*1146+d1'
        assert cm.tuple2string(cm.g2x(149216)) == '*1147'

    def test_g2x_reverse(self):
        """
        Do some g. to c. conversion checking for the gene on the reverse
        strand.
        """
        rna = [2000, 2797, 6405, 6528, 31678, 31750, 47687, 47806, 76448,
               76535, 92453, 92554, 123276, 123470, 146090, 146213]
        cds = [76479, 123290]
        cm = Crossmap(rna, cds, -1)
        # Fix for r536: disable the -u and +d convention.
        #assert cm.tuple2string(cm.g2x(146214)) == '-304-u1'
        assert cm.tuple2string(cm.g2x(146214)) == '-305'
        assert cm.tuple2string(cm.g2x(146091)) == '-182'
        assert cm.tuple2string(cm.g2x(146089)) == '-181+1'
        assert cm.tuple2string(cm.g2x(123291)) == '-1'
        assert cm.tuple2string(cm.g2x(123290)) == '1'
        assert cm.tuple2string(cm.g2x(92555)) == '16-1'
        assert cm.tuple2string(cm.g2x(76479)) == '174'
        assert cm.tuple2string(cm.g2x(76478)) == '*1'
        assert cm.tuple2string(cm.g2x(47807)) == '*32-1'
        assert cm.tuple2string(cm.g2x(47805)) == '*33'
        # Fix for r536: disable the -u and +d convention.
        #assert cm.tuple2string(cm.g2x(1999)) == '*1146+d1'
        assert cm.tuple2string(cm.g2x(1999)) == '*1147'

    def test_x2g_more(self):
        """
        Do some c. to g. conversion checking for the gene on the forward
        strand.
        """
        rna = [5002, 5125, 27745, 27939, 58661, 58762, 74680, 74767, 103409,
               103528, 119465, 119537, 144687, 144810, 148418, 149215]
        cds = [27925, 74736]
        cm = Crossmap(rna, cds, 1)
        assert cm.x2g(-304, -1) == 5001
        assert cm.x2g(-182, 0) == 5124
        assert cm.x2g(-181, 1) ==  5126
        assert cm.x2g(-1, 0) == 27924
        assert cm.x2g(1, 0) == 27925
        assert cm.x2g(16, -1) == 58660
        assert cm.x2g(174, 0) == 74736
        assert cm.x2g(cm.main2int('*1'), 0) == 74737
        assert cm.x2g(cm.main2int('*32'), -1) == 103408
        assert cm.x2g(cm.main2int('*33'), 0) == 103410
        assert cm.x2g(cm.main2int('*1146'), 1) == 149216

    def test_x2g_more_reverse(self):
        """
        Do some c. to g. conversion checking for the gene on the reverse
        strand.
        """
        rna = [2000, 2797, 6405, 6528, 31678, 31750, 47687, 47806, 76448,
               76535, 92453, 92554, 123276, 123470, 146090, 146213]
        cds = [76479, 123290]
        cm = Crossmap(rna, cds, -1)
        assert cm.x2g(-304, -1) == 146214
        assert cm.x2g(-182, 0) == 146091
        assert cm.x2g(-181, 1) ==  146089
        assert cm.x2g(-1, 0) == 123291
        assert cm.x2g(1, 0) == 123290
        assert cm.x2g(16, -1) == 92555
        assert cm.x2g(174, 0) == 76479
        assert cm.x2g(cm.main2int('*1'), 0) == 76478
        assert cm.x2g(cm.main2int('*32'), -1) == 47807
        assert cm.x2g(cm.main2int('*33'), 0) == 47805
        assert cm.x2g(cm.main2int('*1146'), 1) == 1999

    def test_g2x_missing_exons(self):
        """
        Hypothetical gene, missing the first exon and the last two exons
        should have the same crossmapping on the shared part.
        """
        rna1 = [5002, 5125, 27745, 27939, 58661, 58762, 74680, 74767, 103409,
               103528, 119465, 119537, 144687, 144810, 148418, 149215]
        rna2 = [27745, 27939, 58661, 58762, 74680, 74767, 103409, 103528,
                119465, 119537]
        cds = [27925, 74736]
        cm1 = Crossmap(rna1, cds, 1)
        cm2 = Crossmap(rna2, cds, 1)
        assert cm1.g2x(27925) == cm2.g2x(27925)

    def test_g2x_missing_exons_reverse(self):
        """
        Hypothetical gene on the reverse strand, missing the first exon and
        the last two exons should have the same crossmapping on the shared
        part.
        """
        rna1 = [2000, 2797, 6405, 6528, 31678, 31750, 47687, 47806, 76448,
               76535, 92453, 92554, 123276, 123470, 146090, 146213]
        rna2 = [31678, 31750, 47687, 47806, 76448, 76535, 92453, 92554,
                123276, 123470]
        cds = [76479, 123290]
        cm1 = Crossmap(rna1, cds, -1)
        cm2 = Crossmap(rna2, cds, -1)
        assert cm1.g2x(123290) == cm2.g2x(123290)

    def test_splice_sites_noncoding(self):
        """
        Check whether the gene on the forward strand has the right splice
        sites in n. notation.
        """
        rna = [5002, 5125, 27745, 27939, 58661, 58762, 74680, 74767, 103409,
               103528, 119465, 119537, 144687, 144810, 148418, 149215]
        cm = Crossmap(rna, [], 1)
        assert (cm._Crossmap__crossmapping ==
                [1, 124, 125, 319, 320, 421, 422, 509, 510, 629, 630, 702,
                 703, 826, 827, 1624])

    def test_splice_sites_noncoding_reverse(self):
        """
        Check whether the gene on the reverse strand has the right splice
        sites in n. notation.
        """
        rna = [2000, 2797, 6405, 6528, 31678, 31750, 47687, 47806, 76448,
               76535, 92453, 92554, 123276, 123470, 146090, 146213]
        cm = Crossmap(rna, [], -1)
        assert (cm._Crossmap__crossmapping ==
                [1624, 827, 826, 703, 702, 630, 629, 510, 509, 422, 421, 320,
                 319, 125, 124, 1])

    def test_g2x_noncoding(self):
        """
        Do some g. to n. conversion checking for the gene on the forward
        strand.
        """
        rna = [5002, 5125, 27745, 27939, 58661, 58762, 74680, 74767, 103409,
               103528, 119465, 119537, 144687, 144810, 148418, 149215]
        cm = Crossmap(rna, [], 1)
        # Fix for r536: disable the -u and +d convention.
        #assert cm.tuple2string(cm.g2x(5001)) == '1-u1'
        assert cm.tuple2string(cm.g2x(5001)) == '-1'
        assert cm.tuple2string(cm.g2x(5002)) == '1'
        assert cm.tuple2string(cm.g2x(5126)) == '124+1'
        # Fix for r536: disable the -u and +d convention.
        #assert cm.tuple2string(cm.g2x(149216)) == '1624+d1'
        assert cm.tuple2string(cm.g2x(149216)) == '*1'

    def test_g2x_noncoding_reverse(self):
        """
        Do some g. to n. conversion checking for the gene on the reverse
        strand.
        """
        rna = [2000, 2797, 6405, 6528, 31678, 31750, 47687, 47806, 76448,
               76535, 92453, 92554, 123276, 123470, 146090, 146213]
        cm = Crossmap(rna, [], -1)
        # Fix for r536: disable the -u and +d convention.
        #assert cm.tuple2string(cm.g2x(146214)) == '1-u1'
        assert cm.tuple2string(cm.g2x(146214)) == '-1'
        assert cm.tuple2string(cm.g2x(146213)) == '1'
        assert cm.tuple2string(cm.g2x(146089)) == '124+1'
        # Fix for r536: disable the -u and +d convention.
        #assert cm.tuple2string(cm.g2x(1999)) == '1624+d1'
        assert cm.tuple2string(cm.g2x(1999)) == '*1'

    def test_x2g_noncoding(self):
        """
        Do some n. to g. conversion checking for the gene on the forward
        strand.
        """
        rna = [5002, 5125, 27745, 27939, 58661, 58762, 74680, 74767, 103409,
               103528, 119465, 119537, 144687, 144810, 148418, 149215]
        cm = Crossmap(rna, [], 1)
        assert cm.x2g(1, -1) == 5001
        assert cm.x2g(1, 0) == 5002
        assert cm.x2g(124, 1) ==  5126
        assert cm.x2g(1624, 1) == 149216

    def test_x2g_noncoding_reverse(self):
        """
        Do some n. to g. conversion checking for the gene on the reverse
        strand.
        """
        rna = [2000, 2797, 6405, 6528, 31678, 31750, 47687, 47806, 76448,
               76535, 92453, 92554, 123276, 123470, 146090, 146213]
        cm = Crossmap(rna, [], -1)
        assert cm.x2g(1, -1) == 146214
        assert cm.x2g(1, 0) == 146213
        assert cm.x2g(124, 1) ==  146089
        assert cm.x2g(1624, 1) == 1999

    def test_cds_one_exon(self):
        """
        Test a gene that has a CDS that lies entirely in one exon.
        """
        rna = [1, 80, 81, 3719]
        cds = [162, 2123]
        cm = Crossmap(rna, cds, 1)
        assert cm._Crossmap__crossmapping == [-161, -82, -81, 3558]
        assert cm.x2g(1, 0) == 162
        assert cm.tuple2string(cm.g2x(2123)) == '1962'
        assert cm.tuple2string(cm.g2x(2124)) == '*1'

    def test_cds_start_on_splice_site(self):
        """
        Test a gene that has a CDS that starts on an exon splice site.
        """
        rna = [23755059, 23755214, 23777833, 23778028, 23808749, 23808851,
               23824768, 23824856, 23853497, 23853617, 23869553, 23869626,
               23894775, 23894899, 23898506, 23899304]
        cds = [23777833, 23898680]
        cm = Crossmap(rna, cds, 1)
        assert (cm._Crossmap__crossmapping ==
                [-156, -1, 1, 196, 197, 299, 300, 388, 389, 509, 510, 583,
                 584, 708, 709, 1507])
        assert cm.x2g(1, 0) == 23777833
        # Fix for r536: disable the -u and +d convention.
        #assert cm.tuple2string(cm.g2x(2123)) == '-156-u23752936'
        #assert cm.tuple2string(cm.g2x(2124)) == '-156-u23752935'
        assert cm.tuple2string(cm.g2x(2123)) == '-23753092'
        assert cm.tuple2string(cm.g2x(2124)) == '-23753091'

    def test_cds_start_on_splice_site_reverse(self):
        """
        Test a gene on the reverse strand that has a CDS that starts on an
        exon splice site.
        """
        rna = [23777833, 23778028, 23808749, 23808851, 23824768, 23824856,
               23853497, 23853617, 23869553, 23869626, 23894775, 23894899,
               23898506, 23899304]
        cds = [23755214, 23778028]
        cm = Crossmap(rna, cds, -1)
        assert (cm._Crossmap__crossmapping ==
                [196, 1, -1, -103, -104, -192, -193, -313, -314, -387, -388,
                 -512, -513, -1311])

    def test_cds_start_on_splice_site_other(self):
        """
        Test a gene that has a CDS that starts on an other exon splice site.
        """
        rna = [23755059, 23755214, 23777833, 23778028, 23808749, 23808851,
               23824768, 23824856, 23853497, 23853617, 23869553, 23869626,
               23894775, 23894899, 23898506, 23899304]
        cds = [23755214, 23898680]
        cm = Crossmap(rna, cds, 1)
        assert (cm._Crossmap__crossmapping ==
                [-155, 1, 2, 197, 198, 300, 301, 389, 390, 510, 511, 584, 585,
                 709, 710, 1508])

    def test_cds_start_on_splice_site_other_reverse(self):
        """
        Test a gene on the reverse strand that has a CDS that starts on an
        other exon splice site.
        """
        rna = [23777833, 23778028, 23808749, 23808851, 23824768, 23824856,
               23853497, 23853617, 23869553, 23869626, 23894775, 23894899,
               23898506, 23899304]
        cds = [23755214, 23808749]
        cm = Crossmap(rna, cds, -1)
        assert (cm._Crossmap__crossmapping ==
                [197, 2, 1, -102, -103, -191, -192, -312, -313, -386, -387,
                 -511, -512, -1310])

    def test_cds_start_on_transcript_start(self):
        """
        Upstream correction (forward) for CDS start at the start of
        transcript.
        """
        rna = [23777833, 23778028, 23808749, 23808851, 23824768, 23824856,
               23853497, 23853617, 23869553, 23869626, 23894775, 23894899,
               23898506, 23899304]
        cds = [23777833, 23899304]
        cm = Crossmap(rna, cds, 1)
        assert cm.x2g(-1, 0) == cm.x2g(1, -1)

    def test_cds_start_on_transcript_start_reverse(self):
        """
        Upstream correction (reverse) for CDS start at the start of
        transcript.
        """
        rna = [23777833, 23778028, 23808749, 23808851, 23824768, 23824856,
               23853497, 23853617, 23869553, 23869626, 23894775, 23894899,
               23898506, 23899304]
        cds = [23777833, 23899304]
        cm = Crossmap(rna, cds, -1)
        assert cm.x2g(-1, 0) == cm.x2g(1, -1)

    def test_cds_is_exon(self):
        """
        Gene where CDS is exactly one exon.
        """
        rna = [27745, 27939, 58661, 58762, 74680, 74767]
        cds = [58661, 58762]
        cm = Crossmap(rna, cds, 1)
        assert cm._Crossmap__crossmapping == [-195, -1, 1, 102, 103, 190]

    def test_cds_is_exon_reverse(self):
        """
        Gene on the reverse strand where CDS is exactly one exon.
        """
        rna = [27745, 27939, 58661, 58762, 74680, 74767]
        cds = [58661, 58762]
        cm = Crossmap(rna, cds, -1)
        assert cm._Crossmap__crossmapping == [297, 103, 102, 1, -1, -88]
