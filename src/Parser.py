#!/usr/bin/python

from pyparsing import *

class Nomenclatureparser(object) :
    SimpleAlleleVarSet = Forward()
    Nest = Suppress('{') + Group(SimpleAlleleVarSet)("Nest") + \
           Suppress('}')

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
    Nt = Word("acgturykmswbdhvinACGTU", exact = 1)

    # New:
    NtString = Combine(OneOrMore(Nt));

    # Number -> [0-9]+
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
    # RefSeqAcc  -> (GI | AccNo) (`(' GeneSymbol `)')?
    TransVar = Suppress("_v") + Number("TransVar")
    ProtIso = Suppress("_i") + Number("ProtIso")
    GeneSymbol = Suppress('(') + Group(Name("GeneSymbol") + \
                 Optional(TransVar ^ ProtIso))("Gene") + \
                 Suppress(')')
    GI = Suppress(Optional("GI") ^ Optional("GI:")) + \
         Number("RefSeqAcc")
    Version = Suppress('.') + Number("Version")
    AccNo = Combine(Word(alphas + '_') + Number)("RefSeqAcc") + \
            Optional(Version)
    RefSeqAcc = (GI ^ AccNo) + Optional(GeneSymbol)

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
    # PtLoc -> ((`-' | `*')? Number Offset?) | `?'
    PtLoc = Group((Optional(Word("-*", exact = 1))("MainSgn") + \
            Number("Main") + \
            Optional(Offset)) ^ '?')

    # Ref -> ((RefSeqAcc | GeneSymbol) `:')? (`c.' | `g.' | `m.')
    # Changed to:
    # RefType -> (`c' | `g' | `m' | `n' | `r') `.'
    # Ref -> ((RefSeqAcc | GeneSymbol) `:')? RefType?
    RefType = Word("cgmnr", exact = 1)("RefType") + Suppress('.')
    Ref = Optional((RefSeqAcc ^ GeneSymbol) + Suppress(':')) + Optional(RefType)

    # Extent -> PtLoc `_' (`o'? (RefSeqAcc | GeneSymbol) `:')? PtLoc
    # Changed to:
    # Extent -> PtLoc `_' (`o'? (RefSeqAcc | GeneSymbol) `:')? RefType? PtLoc
    Extent = Group(PtLoc("PtLoc"))("StartLoc") + \
             Suppress('_') + Group(Optional(Group(Optional('o') + \
             (RefSeqAcc ^ GeneSymbol) + Suppress(':') + \
             Optional(RefType)))("OptRef") + \
             PtLoc("PtLoc"))("EndLoc")

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
    Subst = Group(PtLoc("PtLoc"))("StartLoc") + \
            Nt + Literal('>')("MutationType") + Nt

    # Del -> Loc `del' (Nt+ | Number)?
    Del = Loc + Literal("del")("MutationType") + \
          Optional(NtString ^ Number)

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
          (NtString ^ Number ^ RangeLoc ^ \
          FarLoc("OptRef")) + Optional(Nest)

    # Indel -> RangeLoc `del' (Nt+ | Number)? 
    #          `ins' (Nt+ | Number | RangeLoc | FarLoc)
    # Changed to:
    # Indel -> RangeLoc `del' (Nt+ | Number)? 
    #          `ins' (Nt+ | Number | RangeLoc | FarLoc) Nest?
    Indel = RangeLoc + Literal("del") + Optional(NtString ^ Number) + \
            Literal("ins")("MutationType") + \
            (NtString ^ Number ^ RangeLoc ^ FarLoc) + \
            Optional(Nest)

    # Inv -> RangeLoc `inv' (Nt+ | Number)?
    # Changed to:
    # Inv -> RangeLoc `inv' (Nt+ | Number)? Nest?
    Inv = RangeLoc + Literal("inv")("MutationType") + \
          Optional(NtString ^ Number) + Optional(Nest)

    # Conv -> RangeLoc `con' FarLoc
    # Changed to:
    # Conv -> RangeLoc `con' FarLoc Nest?
    Conv = RangeLoc + Literal("con")("MutationType") + FarLoc + \
           Optional(Nest)
    
    # ChromBand -> (`p' | `q') Number `.' Number
    # ChromCoords -> `(' Chrom `;' Chrom `)' `(' ChromBand `;' ChromBand `)'
    # TransLoc -> `t' ChromCoords `(' FarLoc `)'
    # RawVar -> Subst | Del | Dup | VarSSR | Ins | Indel | Inv | Conv
    # SingleVar -> Ref RawVar | TransLoc
    # ExtendedRawVar -> RawVar | `=' | `?'
    ChromBand = Suppress(Word("pq", exact = 1)) + Number + Suppress('.') + \
                Number
    ChromCoords = \
        Suppress('(') + Chrom + Suppress(';') + Chrom + Suppress(')') + \
        Suppress('(') + ChromBand + Suppress(';') + ChromBand + Suppress(')')
    TransLoc = Suppress('t') + ChromCoords + Suppress('(') + FarLoc + \
               Suppress(')')
    RawVar = Group(Subst ^ Del ^ Dup ^ VarSSR ^ \
                   Ins ^ Indel ^ Inv ^ Conv)("RawVar")
    SingleVar = Ref + Group(RawVar)("VarSet") ^ TransLoc
    ExtendedRawVar = RawVar ^ '=' ^ '?'

    # New:
    # SimpleAlleleVarSet -> (`[' ExtendedRawVar (`;' ExtendedRawVar)* `]') |
    #                       ExtendedRawVar
    SimpleAlleleVarSet << (Group(Suppress('[') + ExtendedRawVar + \
                          ZeroOrMore(Suppress(';') + ExtendedRawVar) + \
                          Suppress(']'))("SimpleAlleleVarSet") ^ ExtendedRawVar)

    # New:
    # MosaicSet -> (`[' SimpleAlleleVarSet (`/' SimpleAlleleVarSet)* `]') |
    #              SimpleAlleleVarSet
    # ChimeronSet -> (`[' MosaicSet (`//' MosaicSet)* `]') | MosaicSet
    MosaicSet = Group(Suppress('[') + SimpleAlleleVarSet + \
                ZeroOrMore(Suppress('/') + SimpleAlleleVarSet) + \
                Suppress(']'))("MosaicSet") ^ SimpleAlleleVarSet
    #MosaicSet = Forward()
    #MosaicSet << ((Suppress('[') + \
    #                  Group(MosaicSet)("MosaicSet") + \
    #                Suppress(']') + Suppress('/') + Suppress('[') + \
    #                  Group(MosaicSet)("MosaicSet") + \
    #              Suppress(']')) ^ \
    #              Group(Suppress('[') + SimpleAlleleVarSet + \
    #              ZeroOrMore(Suppress('/') + SimpleAlleleVarSet) + \
    #              Suppress(']'))("MosaicSet") ^ \
    #              Group(SimpleAlleleVarSet + ZeroOrMore(Suppress('/') + \
    #              SimpleAlleleVarSet))("MosaicSet") ^ \
    #              SimpleAlleleVarSet)
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
    #Var = StringStart() + (SingleVar ^ MultiVar ^ MultiTranscriptVar ^ \
    #      UnkEffectVar ^ NoRNAVar ^ SplicingVar) + StringEnd()

    def __init__(self) :
        ParserElement.enablePackrat() # Speed up parsing considerably.
        self._state = dict()

    def _parse_command(self, cmd) :
        try :
            return self.Var.parseString(cmd, parseAll = True)
            #return self.Var.parseString(cmd)
        except ParseException, err :
            print "Error: " + str(err)
            pos = int(str(err).split(':')[-1][:-1]) - 1
            print cmd
            print pos * ' ' + '^'
            exit()
        #except
    #_parse_command
#Nomenclatureparser

#def printp(string, depth) :
#    print (depth * "  ") + string
#
#def ppp(parts, depth) :
#    printp("+++ recurse +++", depth)
#    #printp(str(repr(parts)), depth)
#    printp(str(parts), depth)
#    printp("RefSeqAcc: " + str(parts.RefSeqAcc), depth)
#    printp("RefType: " + str(parts.RefType), depth)
#    printp("Version: " + str(parts.Version), depth)
#    if parts.Gene :
#        printp("Gene Symbol: " + str(parts.Gene.GeneSymbol), depth)
#        printp("Transcript variant: " + str(parts.Gene.TransVar), depth)
#        printp("Protein isoform: " + str(parts.Gene.ProtIso), depth)
#    if parts.ChimeronSet :
#        #printp(str(parts.ChimeronSet), depth)
#        for i in parts.ChimeronSet :
#            printp("ChimeronSet", depth)
#            ppp(i, depth + 1)
#    if parts.MosaicSet :
#        #printp(str(parts.MosaicSet), depth)
#        for i in parts.MosaicSet :
#            printp("MosaicSet", depth)
#            ppp(i, depth + 1)
#    if parts.SimpleAlleleVarSet :
#        #printp(str(parts.SimpleAlleleVarSet), depth)
#        for i in parts.SimpleAlleleVarSet :
#            printp("SimpleAlleleVarSet", depth)
#            printp(str(i), depth)
#            #ppp(i, depth + 1)
#    if parts.SingleAlleleVarSet :
#        #printp(str(parts.SingleAlleleVarSet), depth)
#        for i in parts.SingleAlleleVarSet :
#            printp("SingleAlleleVarSet", depth)
#            ppp(i, depth + 1)
#    if parts.MultiAlleleVars :
#        #printp(str(parts.MultiAlleleVars), depth)
#        for i in parts.MultiAlleleVars :
#            printp("MultiAlleleVars", depth)
#            ppp(i, depth + 1)
#            #for j in i :
#            #    printp("Test", depth)
#            #    ppp(j, depth + 1)
#    if parts.RawVar :
#        printp("RawVar", depth)
#        printp(str(parts.RawVar), depth + 1)
#        #for i in parts.RawVar :
#        #    printp(str(i), depth + 1)
#        #    ppp(i, depth + 1)
#    #if
#    """
#    for i in parts.VarSet :
#        printp(i, depth)
#        if i.Nest :
#            ppp(i.Nest, depth + 1)
#        #if
#        if i.StartLoc :
#            if i.StartLoc.OptRef :
#                printp(str(i.StartLoc.OptRef), depth)
#            printp("StartLoc: " + str(i.StartLoc.PtLoc.MainSgn) + \
#            str(i.StartLoc.PtLoc.Main) + \
#            str(i.StartLoc.PtLoc.OffSgn) + str(i.StartLoc.PtLoc.OffOpt) + \
#            str(i.StartLoc.PtLoc.Offset), depth)
#        #if
#        if i.EndLoc :
#            printp(i.EndLoc, depth)
#            if i.EndLoc.OptRef :
#                printp(i.EndLoc.OptRef, depth)
#            printp("EndLoc: " + str(i.EndLoc.PtLoc.MainSgn) + \
#                  str(i.EndLoc.PtLoc.Main) + \
#                  str(i.EndLoc.PtLoc.OffSgn) + str(i.EndLoc.PtLoc.OffOpt) + \
#                  str(i.EndLoc.PtLoc.Offset), depth)
#        #if
#        printp("Mutation type: " + str(i.MutationType), depth)
#    #for
#    """
#    printp("--- recurse ---", depth)
##ppp
#
#import sys
#
#parser = Nomenclatureparser()
#cmd = sys.argv[1]
#print cmd
#ParseObj = parser._parse_command(cmd)
#ppp(ParseObj, 0)
