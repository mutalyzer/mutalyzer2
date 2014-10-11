"""
Tests for the mutalyzer.grammar module.
"""


from __future__ import unicode_literals

#import logging; logging.basicConfig()
import os

import mutalyzer
from mutalyzer.grammar import Grammar
from mutalyzer.output import Output

from utils import MutalyzerTest


class TestGrammar(MutalyzerTest):
    """
    Test the mytalyzer.grammar module.
    """
    def setup(self):
        super(TestGrammar, self).setup()
        self.output = Output(__file__)
        self.grammar = Grammar(self.output)

    def _parse(self, description):
        """
        Parse a variant description.
        """
        self.grammar.parse(description)
        assert self.output.getOutput('parseError') == []

    def test_some_variants(self):
        """
        Some example variants.
        """
        self._parse('NM_002001.2:c.[12del]')
        self._parse('NM_002001.2:c.[(12del)]')
        self._parse('NM_002001.2:c.[(12del)?]')
        self._parse('NM_002001.2:c.[(12del);(12del)]')
        self._parse('NM_002001.2:c.[(12del;12del)]')
        self._parse('NM_002001.2:c.[((12del)?;12del)?]')

    def test_compound_insertion(self):
        """
        Some some compound insertions.
        """
        self._parse('NM_002001.2:c.15_16insA')
        self._parse('NM_002001.2:c.15_16insATC')
        self._parse('NM_002001.2:c.15_16ins[A]')
        self._parse('NM_002001.2:c.15_16ins[ATC]')
        self._parse('NM_002001.2:c.15_16ins28_39')
        self._parse('NM_002001.2:c.15_16ins[28_39]')
        self._parse('NM_002001.2:c.15_16ins[28_39;A]')
        self._parse('NM_002001.2:c.15_16ins[28_39;ATC]')
        self._parse('NM_002001.2:c.15_16ins[28_39;A;ATC]')
        self._parse('NM_002001.2:c.15_16ins28_39inv')
        self._parse('NM_002001.2:c.15_16ins[28_39inv]')
        self._parse('NM_002001.2:c.15_16ins[28_39inv;A]')
        self._parse('NM_002001.2:c.15_16ins[28_39inv;ATC]')
        self._parse('NM_002001.2:c.15_16ins[28_39inv;A;ATC]')

    def test_compound_delins(self):
        """
        Some some compound deletion-insertions.
        """
        self._parse('NM_002001.2:c.12_17delinsA')
        self._parse('NM_002001.2:c.12_17delinsATC')
        self._parse('NM_002001.2:c.12_17delins[A]')
        self._parse('NM_002001.2:c.12_17delins[ATC]')
        self._parse('NM_002001.2:c.12_17delins28_39')
        self._parse('NM_002001.2:c.12_17delins[28_39]')
        self._parse('NM_002001.2:c.12_17delins[28_39;A]')
        self._parse('NM_002001.2:c.12_17delins[28_39;ATC]')
        self._parse('NM_002001.2:c.12_17delins[28_39;A;ATC]')
        self._parse('NM_002001.2:c.12_17delins28_39inv')
        self._parse('NM_002001.2:c.12_17delins[28_39inv]')
        self._parse('NM_002001.2:c.12_17delins[28_39inv;A]')
        self._parse('NM_002001.2:c.12_17delins[28_39inv;ATC]')
        self._parse('NM_002001.2:c.12_17delins[28_39inv;A;ATC]')

    def test_protein_variants(self):
        """
        Some protein variants.
        """
        self._parse('NG_009105.1(OPN1LW):p.=')
        self._parse('NG_009105.1(OPN1LW):p.?')
        self._parse('NM_000076.2(CDKN1C):p.0')
        self._parse('NM_000076.2(CDKN1C):p.0?')
        self._parse('NG_009105.1(OPN1LW):p.(=)')
        self._parse('NM_000076.2(CDKN1C):p.(Ala123del)')
        self._parse('NM_000076.2(CDKN1C):p.(Ala123_Leu126del)')
        self._parse('NM_000076.2(CDKN1C):p.(Ala123_Leu126delinsVal)')
        self._parse('NM_000076.2(CDKN1C):p.Ala123del')
        self._parse('NM_000076.2(CDKN1C):p.Ala123_Leu126del')
        self._parse('NM_000076.2(CDKN1C):p.Ala123_Leu126delinsVal')
        self._parse('NM_000076.2(CDKN1C):p.Ala123_*317delinsVal')
        self._parse('NM_000076.2(CDKN1C):p.Ala123_X317delinsVal')
        self._parse('NM_000076.2(CDKN1C):p.Ala123delinsVal')
        self._parse('NM_000076.2(CDKN1C):p.Ala123delinsValPro')
        self._parse('NM_000076.2(CDKN1C):p.Ala123delinsVP')
        self._parse('NM_000076.2(CDKN1C):p.Ala123fs')
        self._parse('NM_000076.2(CDKN1C_i001):p.(Glu124Serfs*148)')
        self._parse('NM_000076.2(CDKN1C_i001):p.(Glu124SerfsX148)')
        self._parse('NM_000076.2(CDKN1C_i001):p.(E124Sfs*148)')
        self._parse('NM_000076.2(CDKN1C_i001):p.(E124SfsX148)')
        self._parse('NG_009105.1(OPN1LW):p.Met1Leu')
        self._parse('NP_064445.1(OPN1LW):p.Met1?')
        self._parse('NP_064445.1(OPN1LW):p.M1?')
        self._parse('NP_064445.1:p.Gln16del')
        self._parse('NP_064445.1:p.Gln16dup')
        self._parse('NP_064445.1:p.Gln3del')
        self._parse('NP_064445.1:p.Q16del')
        self._parse('NP_064445.1:p.Q16dup')
        self._parse('NP_064445.1:p.Q16*')
        self._parse('NP_064445.1:p.Q16X')
        self._parse('NG_009105.1:p.Gln3Leu')
        self._parse('NG_009105.1(OPN1LW):p.Gln3Leu')
        self._parse('NG_009105.1(OPN1LW_i1):p.Gln3Leu')
        self._parse('NG_009105.1(OPN1LW_v1):p.Gln3Leu')
        self._parse('NG_009105.1(OPN1LW):p.Gln3_Gln4insLeu')
        self._parse('NG_009105.1(OPN1LW):p.Gln3_Gln4insGln')
        self._parse('NG_009105.1(OPN1LW):p.Gln3_Gln4dup')
        self._parse('NG_009105.1(OPN1LW):p.Q3_Q4insQ')
        self._parse('NG_009105.1(OPN1LW):p.Q3_Q4insQQ')
        self._parse('NG_009105.1(OPN1LW):p.Q3_Q4dup')
        self._parse('NG_009105.1(OPN1LW):p.Gln3_Leu7del')
        self._parse('NG_009105.1(OPN1LW):p.Gln3_Leu7delinsValLeu')
        self._parse('NG_009105.1(OPN1LW):p.Gln3_Leu7delinsValPro')
        self._parse('NG_009105.1(OPN1LW):p.Gln3_Leu7delinsGlnGlnTrpSerLeu')
        self._parse('NG_009105.1(OPN1LW):p.Q3_L7delinsGlnGlnTrpSerLeu')
        self._parse('NG_009105.1(OPN1LW):p.Gln3_Leu7delinsQQWSL')
        #self._parse('NG_009105.1(OPN1LW):p.Met1AlaextMet-1')
        #self._parse('NG_009105.1(OPN1LW):p.M1AextM-1')
        #self._parse('NG_009105.1(OPN1LW):p.Gln3_Leu7[3]')
        self._parse('NG_009105.1(OPN1LW):p.Gln3_Leu7(1_6)')
        self._parse('NG_009105.1(OPN1LW):p.Gln3Leu')
        self._parse('NG_009105.1(OPN1LW):p.Gln3Leu')
        #self._parse('NM_000076.2(CDKN1C_i001):p.(*317Trpext*3)')
        self._parse('NM_000076.2(CDKN1C_i001):p.(*317TrpextX3)')
        #self._parse('NM_000076.2(CDKN1C_i001):p.(*317Cysext*1)')
        self._parse('NM_000076.2(CDKN1C_i001):p.(*317CysextX1)')
        #self._parse('NM_000076.2(CDKN1C_i001):p.(*317Cext*1)')
        self._parse('NM_000076.2(CDKN1C_i001):p.(*317CextX1)')
        #self._parse('t(X;17)(DMD:p.Met1_Val1506; SGCA:p.Val250_*387)')

    def test_minus_in_gene_symbol(self):
        """
        Gene symbol is allowed to contain a minus character.
        """
        self._parse('UD_132464528477(KRTAP2-4_v001):c.100del')
