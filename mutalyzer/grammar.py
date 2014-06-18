"""
Module for parsing a variant described using the HGVS nomenclature.

A context-free grammar is defined here, the nomenclature rules are specified
in Backus-Naur Form (BNF), which is used (with some minor modifications) as
source of this module.

The pyparsing [1] module is used to define the grammar and parser.
Documentation for this module seems a bit scattered over the web. For an
informal overview, see [2].

The grammar is described in [3].

@todo: Automatically generate a LaTeX BNF description from this.

[1] http://pyparsing.wikispaces.com/
[2] http://pyparsing.wikispaces.com/HowToUsePyparsing
[3] http://www.biomedcentral.com/1471-2105/12/S4/S5
"""


from pyparsing import *


class Grammar():
    """
    Defines the HGVS nomenclature grammar.

    Some notes about this definition of the grammar follow.

    To create recursion, temporarily define a token using Forward() and later
    on fill in its definition using << instead of =.

    The syntax
        Token('Name')
    is a shorthand for
        Token.setResultsName('Name', listAllMatches=False)
    and you cannot use the shorthand if you want listAllMatches to be True.
    """
    # Forward is used to define a recursive pattern later on (with <<).
    SimpleAlleleVarSet = Forward()

    # BNF: Nest -> `{' SimpleAlleleVarSet `}'
    Nest = Suppress('{') + Group(SimpleAlleleVarSet)('Nest') + Suppress('}')

    ##########################################################################
    # Basic lexemes
    ##########################################################################

    # BNF: Name -> ([a-z] | [a-Z] | [0-9])+
    Name = Word(alphanums, min=1)

    # BNF: Nt -> `a' | `c' | `g' | `u' | `A' | `C' | `G' | `T' | `U'
    #Nt = Word('acgtuACGTU', exact=1)

    # BNF: Nt -> `a' | `c' | `g' | `t' | `u' | `r' | `y' | `k' |
    #            `m' | `s' | `w' | `b' | `d' | `h' | `v' | `n' |
    #            `A' | `C' | `G' | `T' | `U' | `R' | `Y' | `K' |
    #            `M' | `S' | `W' | `B' | `D' | `H' | `V' | `N'
    #
    # For completeness, we should also include `i' here, but since you can
    # then generate `ins' (keyword for insertions) it is left out.
    Nt = Word('acgturykmswbdhvnACGTURYKMSWBDHVN', exact=1)

    # BNF: NtString -> Nt+
    NtString = Combine(OneOrMore(Nt))

    # BNF: Number -> [0-9]+
    Number = Word(nums)

    ##########################################################################
    # Reference sequences
    ##########################################################################

    # BNF: TransVar -> `_v' Number
    TransVar = Suppress('_v') + Number('TransVar')

    # BNF: ProtIso -> `_i' Number
    ProtIso = Suppress('_i') + Number('ProtIso')

    # BNF: GeneName -> ([a-Z] | [0-9] | `-')+
    GeneName = Word(alphanums + '-', min=1)

    # BNF: GeneSymbol -> `(' Name (TransVar | ProtIso)? `)'
    GeneSymbol = Suppress('(') + Group(GeneName('GeneSymbol') + \
                 Optional(TransVar ^ ProtIso))('Gene') + Suppress(')')

    # BNF: GI -> (`GI' | `GI:')? Number
    GI = Suppress(Optional('GI') ^ Optional('GI:') ^ Optional('gi') ^
         Optional('gi:')) + Number('RefSeqAcc')

    # BNF: Version -> `.' Number
    Version = Suppress('.') + Number('Version')

    # BNF: AccNo -> ([a-Z] Number `_')+ Version?
    AccNo = NotAny('LRG_') + \
            Combine(Word(alphas + '_') + Number)('RefSeqAcc') + \
            Optional(Version)

    # BNF: UD -> `UD_' [a-Z]+ (`_' Number)+
    UD = Combine('UD_' + Word(alphas) + OneOrMore('_' + Number))('RefSeqAcc')

    # BNF: LRGTranscriptID -> `t' [0-9]+
    LRGTranscriptID = Suppress('t') + Number('LRGTranscriptID')

    # BNF: LRGProteinID -> `p' [0-9]+
    LRGProteinID = Suppress('p') + Number('LRGProteinID')

    # BNF: LRG -> `LRG' [0-9]+ (`_' (LRGTranscriptID | LRGProteinID))?
    LRG = Combine('LRG_' + Number)('LrgAcc') + Optional(LRGTranscriptID ^ \
                                                        LRGProteinID)

    # RefSeqAcc  -> (GI | AccNo | UD | LRG) (`(' GeneSymbol `)')?
    GenBankRef = (GI ^ AccNo ^ UD) + Optional(GeneSymbol)
    RefSeqAcc = GenBankRef ^ LRG

    # BNF: Chrom -> Name
    Chrom = Name('Chrom')

    # BNF: RefType -> (`c' | `g' | `m' | `n' | `r') `.'
    RefType = Word('cgmnr', exact=1)('RefType') + Suppress('.')

    # BNF: Ref -> ((RefSeqAcc | GeneSymbol) `:')? RefType?
    DRef = Optional((RefSeqAcc ^ GeneSymbol) + Suppress(':')) + Optional(RefType)

    # BNF: RefOne -> RefSeqAcc `:' RefType?
    # Note that RefOne is not included in [3]. We use it in the SingleVar
    # definition below.
    RefOne = RefSeqAcc + Suppress(':') + Optional(RefType)

    ##########################################################################
    # Locations
    ##########################################################################

    # BNF: Offset -> (`+' | `-') (`u' | `d')? (Number | `?')
    Offset = Word('+-', exact=1)('OffSgn') + \
             Optional(Word('ud', exact=1))('OffOpt') + \
             (Number ^ '?')('Offset')

    # BNF: RealPtLoc -> ((`-' | `*')? Number Offset?) | `?'
    RealPtLoc = Group((Optional(Word('-*', exact=1))('MainSgn') + \
                Number('Main') + Optional(Offset)) ^ '?')

    # BNF: IVSLoc -> `IVS' Number (`+' | `-') Number
    IVSLoc = Group(Suppress('IVS') + Number('IVSNumber') + \
             Word('+-', exact=1)('OffSgn') + Number('Offset'))('IVSLoc')

    # BNF: PtLoc -> IVSLoc | RealPtLoc
    PtLoc = IVSLoc ^ RealPtLoc

    # BNF: RealExtent -> PtLoc `_' (`o'? (RefSeqAcc | GeneSymbol) `:')? RefType? PtLoc
    # Note that the optional RefType is not included in [3].
    RealExtent = Group(PtLoc('PtLoc'))('StartLoc') + \
                 Suppress('_') + Group(Optional(Group(Optional('o') + \
                 (RefSeqAcc ^ GeneSymbol) + Suppress(':') + \
                 Optional(RefType)))('OptRef') + \
                 PtLoc('PtLoc'))('EndLoc')

    # BNF: EXLoc -> `EX' Number (`-' Number)?
    EXLoc = Group(Suppress('EX') + Number('EXNumberStart') + \
            Optional(Suppress('-') + Number('EXNumberStop')))('EXLoc')

    # BNF: Extent -> RealExtent | EXLoc
    Extent = RealExtent ^ EXLoc

    # BNF: RangeLoc -> Extent | `(` Extent `)'
    RangeLoc = Extent ^ (Suppress('(') + Extent + Suppress(')'))

    # BNF: Loc -> PtLoc | RangeLoc
    Loc = Group(PtLoc('PtLoc'))('StartLoc') ^ RangeLoc

    # BNF: FarLoc -> (RefSeqAcc | GeneSymbol) (`:' RefType? Extent)?
    FarLoc = (RefSeqAcc ^ GeneSymbol) + Optional(Suppress(':') + \
             Optional(RefType) + Extent)

    # BNF: ChromBand -> (`p' | `q') Number `.' Number
    ChromBand = Suppress(Word('pq', exact=1)) + Number + Suppress('.') + \
                Number

    # BNF: ChromCoords -> `(' Chrom `;' Chrom `)' `(' ChromBand `;' ChromBand `)'
    ChromCoords = \
        Suppress('(') + Chrom + Suppress(';') + Chrom + Suppress(')') + \
        Suppress('(') + ChromBand + Suppress(';') + ChromBand + Suppress(')')

    ##########################################################################
    # Single variations
    ##########################################################################

    # BNF: Subst -> PtLoc Nt `>' Nt
    Subst = Group(PtLoc('PtLoc'))('StartLoc') + Nt('Arg1') + \
            Literal('>').setParseAction(replaceWith('subst'))('MutationType') + \
            Nt('Arg2')

    # BNF: Del -> Loc `del' (Nt+ | Number)?
    Del = Loc + Literal('del')('MutationType') + \
          Optional(NtString ^ Number)('Arg1')

    # BNF: Dup -> Loc `dup' (Nt+ | Number)? Nest?
    Dup = Loc + Literal('dup')('MutationType') + \
          Optional(NtString ^ Number) + Optional(Nest)

    # BNF: AbrSSR -> PtLoc Nt+ `(' Number `_' Number `)'
    AbrSSR = PtLoc + NtString + Suppress('(') + Number + \
             Suppress('_') + Number + Suppress(')')

    # BNF: VarSSR -> (PtLoc Nt+ `[' Number `]') | (RangeLoc `[' Number `]') | AbrSSR
    VarSSR = (PtLoc + NtString + Suppress('[') + Number + \
             Suppress(']')) ^ (RangeLoc + Suppress('[') + Number + \
             Suppress(']')) ^ AbrSSR

    # BNF: Seq -> (Nt+ | Number | RangeLoc `inv'? | FarLoc) Nest?
    Seq = Group(NtString('Sequence') ^ Number ^ \
                (RangeLoc('Range') + Optional(Literal('inv')('Inv'))) ^ \
                FarLoc('OptRef') + Optional(Nest))('Seq')

    # BNF: SeqList -> Seq (`;' Seq)*
    SeqList = delimitedList(Seq, delim=';')('SeqList')

    # BNF: SimpleSeqList -> (`[' SeqList `]') | Seq
    SimpleSeqList = Suppress('[') + SeqList + Suppress(']') ^ Seq

    # BNF: Ins -> RangeLoc `ins' SimpleSeqList
    Ins = RangeLoc + Literal('ins')('MutationType') + SimpleSeqList

    # BNF: Indel -> (RangeLoc | PtLoc) `del' (Nt+ | Number)?
    #          `ins' SimpleSeqList
    # Note that the alternative PtLoc is not included in [3].
    Indel = (RangeLoc ^ Group(PtLoc('PtLoc'))('StartLoc')) + Literal('del') + \
            Optional(NtString ^ Number)('Arg1') + \
            Literal('ins').setParseAction(replaceWith('delins'))('MutationType') + \
            SimpleSeqList

    # BNF: Inv -> RangeLoc `inv' (Nt+ | Number)? Nest?
    Inv = RangeLoc + Literal('inv')('MutationType') + \
          Optional(NtString ^ Number)('Arg1') + Optional(Nest)

    # BNF: Conv -> RangeLoc `con' FarLoc Nest?
    Conv = RangeLoc + Literal('con')('MutationType') + FarLoc + Optional(Nest)

    # BNF: TransLoc -> `t' ChromCoords `(' FarLoc `)'
    TransLoc = Suppress('t') + ChromCoords + Suppress('(') + FarLoc + \
               Suppress(')')

    # The RawVar rule has been changed from [3], were it is given as:
    #
    #     BNF: RawVar -> Subst | Del | Dup | VarSSR | Ins | Indel | Inv | Conv

    # BNF: CRawVar -> Subst | Del | Dup | VarSSR | Ins | Indel | Inv | Conv
    CRawVar = Group(Subst ^ Del ^ Dup ^ VarSSR ^ Ins ^ Indel ^ Inv ^ Conv)('RawVar')

    # We use originalTextFor to retain the original (unparsed) raw variant
    # descriptions. It can be retrieved as element[-1] and is useful for
    # error messages. See:
    # http://packages.python.org/pyparsing/pyparsing.pyparsing-module.html#originalTextFor

    # BNF: RawVar -> (CRawVar | (`(' CRawVar `)')) `?'?
    RawVar = originalTextFor((CRawVar ^ (Suppress('(') + CRawVar + \
                                         Suppress(')'))) + \
                             Suppress(Optional('?')), False)

    # BNF: SingleVar -> Ref RawVar | TransLoc
    # Note that we diverge from [3] in that we use RefOne here instead of Ref.
    SingleVar = RefOne + RawVar ^ TransLoc

    # BNF: ExtendedRawVar -> RawVar | `=' | `?'
    ExtendedRawVar = RawVar ^ '=' ^ '?'

    # BNF: UnkEffectVar -> Ref (`(=)' | `?')
    UnkEffectVar = DRef + (Suppress('(=)') ^ Suppress('?'))

    # BNF: SplicingVar -> Ref (`spl?' | `(spl?)')
    SplicingVar = DRef + (Suppress('spl?') ^ Suppress('(spl?)'))

    # BNF: NoRNAVar -> Ref `0' `?'?
    NoRNAVar = DRef + Suppress('0') + Optional(Suppress('?'))

    ##########################################################################
    # Multiple variations
    ##########################################################################

    # BNF: CAlleleVarSet -> ExtendedRawVar (`;' ExtendedRawVar)*
    CAlleleVarSet = ExtendedRawVar + ZeroOrMore(Suppress(';') + ExtendedRawVar)

    # BNF: UAlleleVarSet -> (CAlleleVarSet | (`(' CAlleleVarSet `)')) `?'?
    UAlleleVarSet = (CAlleleVarSet ^ \
                    (Suppress('(') + CAlleleVarSet + Suppress(')'))) + \
                    Suppress(Optional('?'))

    # BNF: SimpleAlleleVarSet -> (`[' UAlleleVarSet `]') | ExtendedRawVar
    SimpleAlleleVarSet << (Group(Suppress('[') + UAlleleVarSet + \
                          Suppress(']') ^ ExtendedRawVar)('SimpleAlleleVarSet'))

    # BNF: MosaicSet -> (`[' SimpleAlleleVarSet (`/' SimpleAlleleVarSet)* `]') | SimpleAlleleVarSet
    MosaicSet = Group(Suppress('[') + SimpleAlleleVarSet + \
                ZeroOrMore(Suppress('/') + SimpleAlleleVarSet) + \
                Suppress(']'))('MosaicSet') ^ SimpleAlleleVarSet

    # BNF: ChimeronSet -> (`[' MosaicSet (`//' MosaicSet)* `]') | MosaicSet
    ChimeronSet = Group(Suppress('[') + MosaicSet + \
                  ZeroOrMore(Suppress('//') + MosaicSet) + \
                  Suppress(']'))('ChimeronSet') ^ MosaicSet

    # BNF: SingleAlleleVarSet -> (`[` ChimeronSet ((`;' | `^') ChimeronSet)*
    #                            (`(;)' ChimeronSet)* `]') | ChimeronSet
    SingleAlleleVarSet = Group(Suppress('[') + ChimeronSet + \
                         ZeroOrMore((Suppress(';') ^ Suppress('^')) + \
                         ChimeronSet) + ZeroOrMore(Suppress('(;)') + \
                         ChimeronSet) + Suppress(']'))('SingleAlleleVarSet') ^ \
                         ChimeronSet

    # BNF: SingleAlleleVars -> Ref SingleAlleleVarSet
    SingleAlleleVars = DRef + SingleAlleleVarSet

    # BNF: MultiAlleleVars -> Ref SingleAlleleVarSet (`;' Ref? SingleAlleleVarSet)+
    MultiAlleleVars = DRef + Group(SingleAlleleVarSet + \
                      OneOrMore(Suppress(';') + \
                      SingleAlleleVarSet))('MultiAlleleVars')

    # BNF: MultiVar -> SingleAlleleVars | MultiAlleleVars
    MultiVar = SingleAlleleVars ^ MultiAlleleVars

    # BNF: MultiTranscriptVar -> Ref `[` ExtendedRawVar (`;' ExtendedRawVar)*
    #                            (`,' ExtendedRawVar (`;' ExtendedRawVar)*)+ `]'
    MultiTranscriptVar = DRef + Suppress('[') + ExtendedRawVar + \
                         ZeroOrMore(Suppress(';') + ExtendedRawVar) + \
                         OneOrMore(Suppress(',') + ExtendedRawVar + ZeroOrMore(Suppress(';') + \
                         ExtendedRawVar)) + Suppress(']')

    ##########################################################################
    # Protein level variants
    ##########################################################################

    # BNF: Ref -> (Name `:')? `p.'
    #ARef = Optional(Name + Suppress(':')) + Suppress('p.')

    # BNF: Ref -> ((RefSeqAcc | GeneSymbol) `:')? RefType?
    # Note that we diverge from [3] in that we accept RefSeqAcc and GeneSymbol
    # as reference, inherited from the DNA/RNA case.
    PRef = Optional((RefSeqAcc ^ GeneSymbol) + Suppress(':')) + Literal('p')('RefType') + Suppress('.')

    # BNF: AA3 -> `Ala' | `Arg' | `Asn' | `Asp' | `Cys' | `Gln' | `Glu' |
    #             `Gly' | `His' | `Ile' | `Leu' | `Lys' | `Met' | `Phe' |
    #             `Pro' | `Ser' | `Thr' | `Trp' | `Tyr' | `Val'
    AA3 = Literal('Ala') ^ Literal('Arg') ^ Literal('Asn') ^ Literal('Asp') ^ \
          Literal('Cys') ^ Literal('Gln') ^ Literal('Glu') ^ Literal('Gly') ^ \
          Literal('His') ^ Literal('Ile') ^ Literal('Leu') ^ Literal('Lys') ^ \
          Literal('Met') ^ Literal('Phe') ^ Literal('Pro') ^ Literal('Ser') ^ \
          Literal('Thr') ^ Literal('Trp') ^ Literal('Tyr') ^ Literal('Val')

    # BNF: AA1 -> `A' | `R' | `N' | `D' | `C' | `Q' | `E' | `G' | `H' | `I' |
    #             `L' | `K' | `M' | `F' | `P' | `S' | `T' | `W' | `Y' | `V'
    AA1 = Word('ARNDCQEGHILKMFPSTWYV', exact=1)

    # BNF: AA -> AA1 | AA3 | `X'
    # Todo: The '*' notation is preferred.
    AA = AA1 ^ AA3 ^ Word('X*', exact=1)('StopCodon')

    # BNF: PtLoc -> (`-' | `*')? Number | Number (`+' | `-') Number
    PPtLoc = (Optional(Word('-*', exact=1))('MainSgn') + Number('Main')) ^ \
             (Number('Main') + Word('+-', exact=1)('OffSgn') + Number('Offset'))

    # BNF: AAPtLoc -> AA PtLoc
    # Note: AA('Args') is a shorthand for AA.setResultsName('Args', listAllMatches=False)
    AAPtLoc = AA.setResultsName('Args', listAllMatches=True) + PPtLoc

    # BNF: Extent -> AAPtLoc `_' AAPtLoc
    PExtent = AAPtLoc('StartLoc') + Suppress('_') + AAPtLoc('EndLoc')

    # BNF: AARange -> Extent | `(' Extent `)'
    AARange = PExtent ^ (Suppress('(') + PExtent + Suppress(')'))

    # BNF: AALoc -> AAPtLoc | AARange
    AALoc = AAPtLoc('StartLoc') ^ AARange

    # BNF: Subst -> AAPtLoc AA (`extX' `*'? Number)? | (`Met1' | `M1') (`?' | `ext' `-' Number)
    # Todo: 'extX' -> 'ext*' (and loose the optional '*'?)
    # Todo: Optional AA before 'ext' and 'fMet' after 'ext'?
    PSubst = (AAPtLoc + AA.setResultsName('Args') + Optional(Literal('extX') + Optional('*') + Number)) ^ \
             ((Literal('Met1') ^ Literal('M1')) + (Literal('?') ^ (Literal('ext') + Literal('-') + Number)))

    # BNF: Del -> AALoc `del'
    PDel = AALoc + Literal('del')('MutationType')

    # BNF: Dup -> AALoc `dup'
    PDup = AALoc + Literal('dup')('MutationType')

    # BNF: VarSSR -> AALoc `(' Number `_' Number `)'
    PVarSSR = AALoc + Suppress('(') + Number + Suppress('_') + Number + Suppress(')')

    # BNF: Ins -> AARange `ins' (AA+ | Number)
    PIns = AARange + Literal('ins')('MutationType') + (OneOrMore(AA) ^ Number)('Args')

    # BNF: Indel -> AALoc `delins' (AA+ | Number)
    PIndel = AALoc + Literal('delins')('MutationType') + (OneOrMore(AA) ^ Number)('Args')

    # BNF: ShortFS -> AAPtLoc `fs'
    ShortFS = AAPtLoc + Suppress('fs')

    # BNF: LongFS -> AAPtLoc AA `fs' `X' Number
    # Todo: 'X' and '*'?
    LongFS = AAPtLoc + AA('Args') + Suppress('fs') + Suppress(Word('X*', exact=1)) + Number

    # BNF: FrameShift -> ShortFS | LongFS
    FrameShift = ShortFS ^ LongFS

    # The ARawVar rule has been changed from [3], were it is given as:
    #
    #     BNF: RawVar -> Subst | Del | Dup | VarSSR | Ins | Indel | FrameShift |
    #                    `=' | `?' | `0' | `0?'

    # BNF: CRawVar -> Subst | Del | Dup | VarSSR | Ins | Indel | FrameShift |
    #                 `=' | `?' | `0' | `0?'
    PCRawVar = Group(PSubst ^ PDel ^ PDup ^ PVarSSR ^ PIns ^ PIndel ^ FrameShift ^ \
                     Literal('=') ^ Literal('?') ^ Literal('0') ^ Literal('0?'))('RawVar')

    # BNF: RawVar -> CRawVar | `(' CRawVar `)'
    PRawVar = PCRawVar ^ (Suppress('(') + PCRawVar + Suppress(')'))

    # BNF: SingleVar -> Ref RawVar
    PSingleVar = PRef + PRawVar

    # BNF: UnkAlleleVars -> Ref `[' RawVar `(;)' RawVar `]'
    PUnkAlleleVars = PRef + Suppress('[') + PRawVar + Suppress('(;)') + PRawVar + Suppress(']')

    # BNF: SingleAlleleVarSet -> `[' RawVar (`;' RawVar)+  | (`,' RawVar)+ `]'
    PSingleAlleleVarSet = Suppress('[') + PRawVar + \
                          (OneOrMore(Suppress(';') + PRawVar) ^ OneOrMore(Suppress(',') + PRawVar)) + \
                          Suppress(']')

    # BNF: MultiAlleleVars -> Ref SingleAlleleVarSet `;' Ref? SingleAlleleVarSet
    PMultiAlleleVars = PRef + PSingleAlleleVarSet + Suppress(';') + Optional(PRef) + PSingleAlleleVarSet

    # BNF: SingleAlleleVars -> Ref SingleAlleleVarSet
    PSingleAlleleVars = PRef + PSingleAlleleVarSet

    # BNF: MultiVar -> SingleAlleleVars | MultiAlleleVars | UnkAlleleVars
    PMultiVar = PSingleAlleleVars ^ PMultiAlleleVars ^ PUnkAlleleVars

    # BNF: ProteinVar -> SingleVar | MultiVar
    ProteinVar = PSingleVar ^ PMultiVar

    ##########################################################################
    # Top-level rule
    ##########################################################################

    # BNF: Var -> SingleVar | MultiVar | MultiTranscriptVar |
    #             UnkEffectVar | NoRNAVar | SplicingVar
    Var = SingleVar ^ MultiVar ^ MultiTranscriptVar ^ \
          UnkEffectVar ^ NoRNAVar ^ SplicingVar ^ ProteinVar

    def __init__(self, output):
        """
        Initialise the class and enable packrat parsing. Packrat speeds up
        parsing considerably.

        @arg output: The output object.
        @type output: mutalyzer.output.Output
        """
        self._output = output
        ParserElement.enablePackrat()
    #__init__

    def parse(self, variant):
        """
        Parse the input string and return a parse tree if the parsing was
        successful. Otherwise print the parse error and the position in
        the input where the error occurred (and return None).

        @arg variant: The input string that needs to be parsed.
        @type variant: string

        @return: The parse tree containing the parse results, or None in
                 case of a parsing error.
        @rtype: pyparsing.ParseResults

        @todo: Use information in ParseException as described here:
            http://pyparsing.wikispaces.com/HowToUsePyparsing
        """
        try:
            return self.Var.parseString(variant, parseAll=True)
            # Todo: check .dump()
        except ParseException as err:
            print err.line
            print " "*(err.column-1) + "^"
            print err
            # Log parse error and the position where it occurred.
            self._output.addMessage(__file__, 4, 'EPARSE', str(err))
            pos = int(str(err).split(':')[-1][:-1]) - 1
            self._output.addOutput('parseError', variant)
            self._output.addOutput('parseError', pos * ' ' + '^')
            return None
    #parse
#Grammar
