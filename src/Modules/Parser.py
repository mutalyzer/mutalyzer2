#!/usr/bin/python

"""
    Module for parting a variant described using the HGVS nomenclature.

    A context-free parser is defined here, the nomenclature rules are specified
    in BNF, which is used (with some minor modifications) as source of this
    module.

    Public classes:
        Nomenclatureparser ; Parse an input string.
"""

from pyparsing import *

class Nomenclatureparser() :
    """
        Parse an input string.

        Private variables:
            __output ; The output object.

        Public variables:
            All variables defined below, they are all context-free grammar 
            rules.
        
        Special methods:
            __init__() ; Initialise the class and enable packrat parsing.

        Public methods:
            parse(input) ; Parse the input string and return a parse tree.
    """
    
    # New:
    # Nest -> `{' SimpleAlleleVarSet `}'
    SimpleAlleleVarSet = Forward()
    Nest = Suppress('{') + Group(SimpleAlleleVarSet)("Nest") + Suppress('}')

    # New:
    # Name -> ([a-z][a-Z][0-9])*
    Name = Word(alphanums, min = 1)

    # Nt -> `A' | `C' | `G' | `T'
    # and
    # Nt -> `a' | `c' | `g' | `u'
    # Changed to:
    # Nt -> `a' | `c' | `g' | `t' | `u' | `r' | `y' | `k' | 
    #       `m' | `s' | `w' | `b' | `d' | `h' | `v' | `i' | 
    #       `n' | `A' | `C' | `G' | `T' | `U' 
    Nt = Word("acgtuACGTU", exact = 1)

    # New:
    NtString = Combine(OneOrMore(Nt))

    # Number -> [0-9]+
    #Number = Combine(OneOrMore((nums)))
    Number = Word(nums)

    # RefSeqAcc  -> [A-Z][A-Z0-9_]+ (`.' [0-9]+)?
    # GeneSymbol -> [A-Z][A-Z0-9_]+
    # Changed to:
    # TransVar   -> `_v' Number
    # ProtIso    -> `_i' Number
    # GeneSymbol -> `(' Name (TransVar | ProtVar)? `)'
    # GI         -> (`GI' | `GI:')? Number
    # Version    -> `.' Number
    # AccNo      -> ([a-Z] Number `_')+ Version? 
    # UD         -> `UD_' [a-Z]+ (`_' Number)+
    TransVar = Suppress("_v") + Number("TransVar")
    ProtIso = Suppress("_i") + Number("ProtIso")
    GeneSymbol = Suppress('(') + Group(Name("GeneSymbol") + \
                 Optional(TransVar ^ ProtIso))("Gene") + Suppress(')')
    GI = Suppress(Optional("GI") ^ Optional("GI:")) + Number("RefSeqAcc")
    Version = Suppress('.') + Number("Version")
    AccNo = NotAny("LRG_") + \
            Combine(Word(alphas + '_') + Number)("RefSeqAcc") + \
            Optional(Version)
    UD = Combine("UD_" + Word(alphas) + OneOrMore('_' + Number))("RefSeqAcc")

    # LRGTranscriptID -> `t' [0-9]+
    # LRGProteinID    -> `p' [0-9]+
    # LRG             -> `LRG' [0-9]+ (`_' (LRGTranscriptID | LRGProteinID))?
    LRGTranscriptID = Suppress('t') + Number("LRGTranscriptID")
    LRGProteinID = Suppress('p') + Number("LRGProteinID")
    LRG = Combine("LRG_" + Number)("LrgAcc") + Optional(LRGTranscriptID ^ \
                                                        LRGProteinID)

    # RefSeqAcc  -> (GI | AccNo | UD | LRG) (`(' GeneSymbol `)')?
    GenBankRef = (GI ^ AccNo ^ UD) + Optional(GeneSymbol)
    RefSeqAcc = GenBankRef ^ LRG

    # Chrom -> `1'..`22' | `X' | `Y'
    # Changed to:
    # Chrom -> Name
    Chrom = Name("Chrom")

    # AbsLoc -> Number
    #AbsLoc = Number("AbsLoc")

    # Offset -> (`+' | `-') (Number | `?')
    # Changed to:
    # Offset -> (`+' | `-') (`u' | `d')? (Number | `?')
    Offset = Word("+-", exact = 1)("OffSgn") + \
             Optional(Word("ud", exact = 1))("OffOpt") + \
             (Number ^ '?')("Offset")

    # PtLoc -> AbsLoc | `-' AbsLoc Offset? | AbsLoc Offset | `*' AbsLoc Offset?
    # Changed to:
    # RealPtLoc -> ((`-' | `*')? Number Offset?) | `?'
    # IVSLoc -> `IVS' Number (`+' | `-') Number
    # PtLoc -> ((`-' | `*')? Number Offset?) | `?'
    RealPtLoc = Group((Optional(Word("-*", exact = 1))("MainSgn") + \
                Number("Main") + Optional(Offset)) ^ '?')
    IVSLoc = Group(Suppress("IVS") + Number("IVSNumber") + \
             Word("+-", exact = 1)("OffSgn") + Number("Offset"))("IVSLoc")
    PtLoc = IVSLoc ^ RealPtLoc

    # Ref -> ((RefSeqAcc | GeneSymbol) `:')? (`c.' | `g.' | `m.')
    # Changed to:
    # RefType -> (`c' | `g' | `m' | `n' | `r') `.'
    # Ref -> ((RefSeqAcc | GeneSymbol) `:')? RefType?
    RefType = Word("cgmnr", exact = 1)("RefType") + Suppress('.')
    Ref = Optional((RefSeqAcc ^ GeneSymbol) + Suppress(':')) + Optional(RefType)
    RefOne = RefSeqAcc + Suppress(':') + Optional(RefType)

    # Extent -> PtLoc `_' (`o'? (RefSeqAcc | GeneSymbol) `:')? PtLoc
    # Changed to:
    # Extent -> PtLoc `_' (`o'? (RefSeqAcc | GeneSymbol) `:')? RefType? PtLoc
    RealExtent = Group(PtLoc("PtLoc"))("StartLoc") + \
                 Suppress('_') + Group(Optional(Group(Optional('o') + \
                 (RefSeqAcc ^ GeneSymbol) + Suppress(':') + \
                 Optional(RefType)))("OptRef") + \
                 PtLoc("PtLoc"))("EndLoc")
    EXLoc = Group(Suppress("EX") + Number("EXNumberStart") + \
            Optional(Suppress('-') + Number("EXNumberStop")))("EXLoc")
    Extent = RealExtent ^ EXLoc            

    # RangeLoc -> Extent | `(` Extent `)'
    # Loc -> PtLoc | RangeLoc
    RangeLoc = Extent ^ (Suppress('(') + Extent + Suppress(')'))
    Loc = Group(PtLoc("PtLoc"))("StartLoc") ^ \
          RangeLoc

    # FarLoc -> Ref Loc
    # Changed to:
    # FarLoc -> (RefSeqAcc | GeneSymbol) (`:' RefType? Extent)? 
    FarLoc = (RefSeqAcc ^ GeneSymbol) + Optional(Suppress(':') + \
             Optional(RefType) + Extent)

    # Subst -> PtLoc Nt `>' Nt
    Subst = Group(PtLoc("PtLoc"))("StartLoc") + Nt("Arg1") + \
        Literal('>').setParseAction(replaceWith("subst"))("MutationType") + \
        Nt("Arg2")

    # Del -> Loc `del' (Nt+ | Number)?
    Del = Loc + Literal("del")("MutationType") + \
          Optional(NtString ^ Number)("Arg1")

    # Dup -> Loc `dup' (Nt+ | Number | RangeLoc | FarLoc)?
    # Changed to:
    # Dup -> Loc `dup' (Nt+ | Number)? Nest?
    Dup = Loc + Literal("dup")("MutationType") + \
          Optional(NtString ^ Number) + Optional(Nest)

    # VarSSR -> PtLoc `(' Nt+ `)' Number `_' Number
    # Changed to:
    # AbrSSR -> PtLoc  Nt+ `(' Number `_' Number `)'
    # VarSSR -> (PtLoc  Nt+ `[' Number `]') | (RangeLoc `[' Number `]') 
    AbrSSR = PtLoc + NtString + Suppress('(') + Number + \
             Suppress('_') + Number + Suppress(')')
    VarSSR = (PtLoc + NtString + Suppress('[') + Number + \
             Suppress(']')) ^ (RangeLoc + Suppress('[') + Number + \
             Suppress(']')) ^ AbrSSR

    # Ins -> RangeLoc `ins' (Nt+ | Number | RangeLoc | FarLoc)
    # Changed to:
    # Ins -> RangeLoc `ins' (Nt+ | Number | RangeLoc | FarLoc) Nest?
    Ins = RangeLoc + Literal("ins")("MutationType") + \
          ((NtString ^ Number)("Arg1") ^ RangeLoc ^ \
          FarLoc("OptRef")) + Optional(Nest)

    # Indel -> RangeLoc `del' (Nt+ | Number)? 
    #          `ins' (Nt+ | Number | RangeLoc | FarLoc)
    # Changed to:
    # Indel -> (RangeLoc | PtLoc) `del' (Nt+ | Number)? 
    #          `ins' (Nt+ | Number | RangeLoc | FarLoc) Nest?
    Indel = (RangeLoc ^ Group(PtLoc("PtLoc"))("StartLoc")) + Literal("del") + \
        Optional(NtString ^ Number)("Arg1") + \
        Literal("ins").setParseAction(replaceWith("delins"))("MutationType") + \
        ((NtString ^ Number)("Arg2") ^ RangeLoc ^ FarLoc) + Optional(Nest)

    # Inv -> RangeLoc `inv' (Nt+ | Number)?
    # Changed to:
    # Inv -> RangeLoc `inv' (Nt+ | Number)? Nest?
    Inv = RangeLoc + Literal("inv")("MutationType") + \
          Optional(NtString ^ Number)("Arg1") + Optional(Nest)

    # Conv -> RangeLoc `con' FarLoc
    # Changed to:
    # Conv -> RangeLoc `con' FarLoc Nest?
    Conv = RangeLoc + Literal("con")("MutationType") + FarLoc + \
           Optional(Nest)
    
    # ChromBand -> (`p' | `q') Number `.' Number
    # ChromCoords -> `(' Chrom `;' Chrom `)' `(' ChromBand `;' ChromBand `)'
    # TransLoc -> `t' ChromCoords `(' FarLoc `)'
    ChromBand = Suppress(Word("pq", exact = 1)) + Number + Suppress('.') + \
                Number
    ChromCoords = \
        Suppress('(') + Chrom + Suppress(';') + Chrom + Suppress(')') + \
        Suppress('(') + ChromBand + Suppress(';') + ChromBand + Suppress(')')
    TransLoc = Suppress('t') + ChromCoords + Suppress('(') + FarLoc + \
               Suppress(')')

    # RawVar -> Subst | Del | Dup | VarSSR | Ins | Indel | Inv | Conv
    # Changed to:
    # CRawVar -> Subst | Del | Dup | VarSSR | Ins | Indel | Inv | Conv
    # RawVar -> (CRawVar | (`(' CRawVar `)')) `?'?
    CRawVar = Group(Subst ^ Del ^ Dup ^ VarSSR ^ \
                   Ins ^ Indel ^ Inv ^ Conv)("RawVar")
    RawVar = (CRawVar ^ (Suppress('(') + CRawVar + Suppress(')'))) + \
             Suppress(Optional('?'))

    # SingleVar -> Ref RawVar | TransLoc
    # ExtendedRawVar -> RawVar | `=' | `?'
    SingleVar = RefOne + RawVar ^ TransLoc
    ExtendedRawVar = RawVar ^ '=' ^ '?'

    # New:
    # CAlleleVarSet -> ExtendedRawVar (`;' ExtendedRawVar)*
    # UAlleleVarSet -> (CAlleleVarSet | (`(' CAlleleVarSet `)')) `?'?
    # SimpleAlleleVarSet -> (`[' UAlleleVarSet `]') |
    #                       ExtendedRawVar
    CAlleleVarSet = ExtendedRawVar + ZeroOrMore(Suppress(';') + ExtendedRawVar)
    UAlleleVarSet = (CAlleleVarSet ^ \
                    (Suppress('(') + CAlleleVarSet + Suppress(')'))) + \
                    Suppress(Optional('?'))
    SimpleAlleleVarSet << (Group(Suppress('[') + UAlleleVarSet + \
                          Suppress(']') ^ ExtendedRawVar)("SimpleAlleleVarSet"))

    # New:
    # MosaicSet -> (`[' SimpleAlleleVarSet (`/' SimpleAlleleVarSet)* `]') |
    #              SimpleAlleleVarSet
    # ChimeronSet -> (`[' MosaicSet (`//' MosaicSet)* `]') | MosaicSet
    MosaicSet = Group(Suppress('[') + SimpleAlleleVarSet + \
                ZeroOrMore(Suppress('/') + SimpleAlleleVarSet) + \
                Suppress(']'))("MosaicSet") ^ SimpleAlleleVarSet
    ChimeronSet = Group(Suppress('[') + MosaicSet + \
                ZeroOrMore(Suppress("//") + MosaicSet) + \
                Suppress(']'))("ChimeronSet") ^ MosaicSet

    # SingleAlleleVarSet -> `[` ExtendedRawVar (`;' ExtendedRawVar)+ `]'
    # UnkAlleleVars -> Ref `[` RawVar `(+)' RawVar `]'
    # Changed to:
    # SingleAlleleVarSet -> (`[` ChimeronSet ((`;' | `^') ChimeronSet)* 
    #                        (`(;)' ChimeronSet)* `]') | ChimeronSet
    SingleAlleleVarSet = Group(Suppress('[') + ChimeronSet + \
                         ZeroOrMore((Suppress(';') ^ Suppress('^')) + \
                         ChimeronSet) + ZeroOrMore(Suppress("(;)") + \
                         ChimeronSet) + Suppress(']'))("SingleAlleleVarSet") ^ \
                         ChimeronSet

    # SingleAlleleVars -> Ref SingleAlleleVarSet
    SingleAlleleVars = Ref + SingleAlleleVarSet

    # MultiAlleleVars -> Ref SingleAlleleVarSet `+' Ref? SingleAlleleVarSet
    # Changed to:
    # MultiAlleleVars -> Ref SingleAlleleVarSet (`;' Ref? SingleAlleleVarSet)+
    MultiAlleleVars = Ref + Group(SingleAlleleVarSet + \
                      OneOrMore(Suppress(';') + \
                      SingleAlleleVarSet))("MultiAlleleVars")

    # MultiVar -> SingleAlleleVars | MultiAlleleVars | UnkAlleleVars
    # Changed to:
    # MultiVar -> SingleAlleleVars | MultiAlleleVars
    MultiVar = SingleAlleleVars ^ MultiAlleleVars

    # MultiTranscriptVar -> Ref `[` ExtendedRawVar (`;' ExtendedRawVar)* 
    #                       (`,' ExtendedRawVar (`;' ExtendedRawVar)*)+ `]' 
    MultiTranscriptVar = Ref + Suppress('[') + ExtendedRawVar + \
        ZeroOrMore(Suppress(';') + ExtendedRawVar) + \
        OneOrMore(Suppress(',') + ExtendedRawVar + ZeroOrMore(Suppress(';') + \
        ExtendedRawVar)) + Suppress(']')

    # UnkEffectVar -> Ref `(=)' | Ref `?'
    # Changed to:
    # UnkEffectVar -> Ref (`(=)' | `?')
    UnkEffectVar = Ref + (Suppress("(=)") ^ Suppress('?'))

    # SplicingVar -> Ref `spl?' | Ref `(spl?)'
    # Changed to:
    # SplicingVar -> Ref (`spl?' | `(spl?)')
    SplicingVar = Ref + (Suppress("spl?") ^ Suppress("(spl?)"))

    # NoRNAVar -> Ref `0' `?'?
    NoRNAVar = Ref + Suppress('0') + Optional(Suppress('?'))

    # DNAVar -> SingleVar | MultiVar
    # RNAVar -> SingleVar | MultiVar | MultiTranscriptVar | 
    #           UnkEffectVar | NoRNAVar | SplicingVar
    # Changed to:
    # Var -> SingleVar | MultiVar | MultiTranscriptVar | 
    #        UnkEffectVar | NoRNAVar | SplicingVar
    Var = SingleVar ^ MultiVar ^ MultiTranscriptVar ^ \
          UnkEffectVar ^ NoRNAVar ^ SplicingVar

    def __init__(self, output) :
        """
            Initialise the class and enable packrat parsing.

            Arguments:
                output ; The output object.

            Private variables (altered):
                __output ; Set to the output object.
        """

        self.__output = output
        ParserElement.enablePackrat() # Speed up parsing considerably.
    #__init__

    def parse(self, variant) :
        """
            Parse the input string and return a parse tree if the parsing was
            successful. Otherwise print the parse error and the position in 
            the input where the error occurred.
            
            Arguments:
                variant ; The input string that needs to be parsed.

            Private variables:
                __output ; The output object.

            Public variables:
                Var ; The top-level rule of our parser.

            Returns:
                Object ; The parse tree containing the parse results.
        """

        try :
            return self.Var.parseString(variant, parseAll = True)
        except ParseException, err :
            self.__output.addMessage(__file__, 4, "EPARSE", str(err))

            # Log the input.
            #self.__output.addMessage(__file__, 4, "EPARSE", variant)

            # And log the position where the parsing error occurred.
            pos = int(str(err).split(':')[-1][:-1]) - 1
            #self.__output.addMessage(__file__, 4, "EPARSE", pos * ' ' + '^')

            #self.__output.addOutput("parseError", str(err))
            self.__output.addOutput("parseError", variant)
            self.__output.addOutput("parseError", pos * ' ' + '^')
            return None
        #except
    #parse
#Nomenclatureparser

#
# Unit test.
#
if __name__ == "__main__" :
    P = Nomenclatureparser()
    parsetree = P.parse("NM_002001.2:c.[12del]")
    parsetree = P.parse("NM_002001.2:c.[(12del)]")
    parsetree = P.parse("NM_002001.2:c.[(12del)?]")
    parsetree = P.parse("NM_002001.2:c.[(12del);(12del)]")
    parsetree = P.parse("NM_002001.2:c.[(12del;12del)]")
    parsetree = P.parse("NM_002001.2:c.[((12del)?;12del)?]")
    del P
#if
