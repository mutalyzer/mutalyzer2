#!/usr/bin/python

"""
    The nomenclature checker.
"""

import sys
import math
import types
import Bio

import Bio.Seq
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import seq3

from Modules import Retriever
from Modules import GenRecord
from Modules import Crossmap
from Modules import Parser
from Modules import Db
from Modules import Mutator
from Modules import Output
from Modules import Config

class newMut() :
    def __init__(self) :
        self.c = ""
        self.g = ""
    #__init__
#newMut

def __order(a, b) :
    if a < b :
        return a, b
    return b, a
#__order

def __roll(string, start, stop, orientation) :
    """
    """

    pattern = string[start:stop]
    if orientation == 1 :
        i = stop - 1

        while i < len(string) and string[i] == pattern[-1] :
            pattern = pattern[1:] + pattern[0]
            i += 1
        #while

        return i - len(pattern) + 1
    #if
    else :
        i = start

        while i > -1 and string[i] == pattern[0] :
            pattern = pattern[-1] + pattern[:-1]
            i -= 1
        #while

        return i + 2
    #else
#__roll

def roll2(ref, start, stop) :
    """
    """

    pattern = ref[start:stop]
    patternLength = len(pattern)

    g_min = start - 1
    j = patternLength - 1
    while g_min > -1 and ref[g_min] == pattern[j % patternLength] :
        j -= 1
        g_min -= 1
    #while 
    g_min += 1

    g_max = stop
    j = 0
    while g_max < len(ref) and ref[g_max] == pattern[j % patternLength] :
        j += 1
        g_max += 1
    #while 

    return g_min, g_max - patternLength
#roll2

def __palinsnoop(string) :
    """
    """

    revcomp = Bio.Seq.reverse_complement(string)

    for i in range(int(math.ceil(len(string) / 2.0))) :
        if string[i] != revcomp[i] :
            return i # The first i elements are palindromic.
    return -1        # Perfect palindrome.
#__palinsnoop

def __bprint(s) :
    """
    """

    if not s :
        return

    block = 10
    line = 6 * block

    m = int(math.floor(math.log(len(s), 10)) + 1)
    o = 1
    print "%s " % str(o).rjust(m),
    for i in range(0, len(s), block) :
        print s[i:i + block],
        if not (i + block) % line and i + block < len(s) :
            o += line
            print "\n%s " % str(o).rjust(m),
        #if
    #for
#__bprint

def __PtLoc2main(Loc) :
    """
    """

    main = int(Loc.Main)
    if Loc.MainSgn == '-' :
        main = -main

    return main
#__PtLoc2main

def __PtLoc2offset(Loc) :
    """
    """

    if Loc.Offset :
        offset = int(Loc.Offset)
        if Loc.OffSgn == '-' :
            offset = -offset
    #if
    else :
        offset = 0

    return offset
#__PtLoc2offset

def __splice(string, splice_sites) :
    """
        Construct the transcript or the coding sequence from a record and
        a list of splice sites.

        Arguments:
            record       ; A GenBank record (see the BioPython documentation).
            splice_sites ; A list of even length of integers.

        Returns:
            String ; The concatenation of slices from the sequence that is 
                     present in the GenBank record.
    """

    transcript = ""

    for i in range(0, len(splice_sites), 2) :
        transcript += string[splice_sites[i] - 1:splice_sites[i + 1]] 

    return transcript
#__splice

def __nsplice(string, splice_sites, CDS, orientation) :
    """
    """

    transcript = ""

    if orientation == 1 :
        for i in range(0, len(splice_sites), 2) :
            if CDS[0] >= splice_sites[i] and CDS[0] <= splice_sites[i + 1] :
                transcript += string[CDS[0] - 1:splice_sites[i + 1]]
            else :
                if splice_sites[i] > CDS[0] :
                    transcript += \
                        string[splice_sites[i] - 1:splice_sites[i + 1]] 
        #for
    #if        
    else :
        for i in range(0, len(splice_sites), 2) :
            if CDS[1] >= splice_sites[i] and CDS[1] <= splice_sites[i + 1] :
                transcript += string[splice_sites[i] - 1: CDS[1]]
            else :
                if splice_sites[i] < CDS[1] :
                    transcript += \
                     string[splice_sites[i] - 1:splice_sites[i + 1]] 
        #for
    #else        

    return transcript
#__nsplice

def __checkOptArg(ref, p1, p2, arg, M, O) :
    """
    """

    if arg :
        if arg.isdigit() :
            length = int(arg)
            interval = p2 - p1 + 1
            if length != interval :
                O.addMessage(__file__, 3, "EARGLEN", 
                    "The length (%i) differed from that of the range (%i)." % (
                    length, interval))
                return False
            #if
        #if
        else :
            ref_slice = str(ref[p1 - 1:p2])
            if ref_slice != str(arg) : # FIXME more informative.
                O.addMessage(__file__, 3, "EREF", 
                    "%s not found at position c.%s (g.%i), found %s instead." \
                    % (arg, M.g2c(p1), p1, ref_slice))
                return False
            #if
        #else
    #if
    return True
#__checkOptArg

def __lcp(str1, str2) :
    """
        Calculate the length of the longest common prefix of two strings.

        Arguments:
            str1 ; The first string.
            str2 ; The second string.

        Returns:
            integer ; The length of the longest common prefix of str1 and str2.
    """

    pos = 0
    s1l = len(str1) # Use the lengths to make sure we don't exceed the length
    s2l = len(str2) # of the strings.

    while pos < s1l and pos < s2l and str1[pos] == str2[pos] :
        pos += 1

    return pos
#__lcp

def __lcs(str1, str2) :
    """
        Calculate the length of the longest common suffix of two strings.

        Arguments:
            str1 ; The first string.
            str2 ; The second string.

        Returns:
            integer ; The length of the longest common suffix of str1 and str2.
    """

    t1 = str1[::-1] # Invert str1.
    t2 = str2[::-1] # Invert str2.

    # The lcp of the two inverted strings is the lcs of the original strings.
    return __lcp(t1, t2) 
#__lcs

def findInFrameDescription(str1, str2) :
    """
    """

    # Nothing happened.
    if str1 == str2 :
        return "p.(=)"

    lcp = __lcp(str1, str2)
    lcs = __lcs(str1[lcp:], str2[lcp:])
    str1_end = len(str1) - lcs
    str2_end = len(str2) - lcs

    # Insertion / Duplication.
    if not str1_end - lcp :
        inLen = str2_end - lcp

        if lcp - inLen >= 0 and str1[lcp - inLen:lcp] == str2[lcp:str2_end] :
            if inLen == 1 :
                return "p.(%s%idup)" % (seq3(str1[lcp - inLen]), 
                                        lcp - inLen + 1)
            return "p.(%s%i_%s%idup)" % (seq3(str1[lcp - inLen]), 
                                         lcp - inLen + 1,
                                         seq3(str1[lcp - 1], lcp))
        #if
        return "p.(%s%i_%s%iins%s)" % (seq3(str1[lcp - 1]), lcp, 
                                       seq3(str1[lcp]), lcp + 1,
                                       seq3(str2[lcp:str2_end]))
    #if

    # Deletion.
    if not str2_end - lcp :
        if lcp + 1 == str1_end :
            return "p.(%s%idel)" % (seq3(str1[lcp], lcp + 1))
        return "p.(%s%i_%s%idel)" % (seq3(str1[lcp - 1]), lcp + 1, 
                                     seq3(str1[str1_end - 1]), str1_end)
    #if

    # Substitution.
    if str1_end == str2_end and str1_end == lcp + 1 :
        if len(str1) > len(str2) :
            return "p.(*%i%sext*%i)" % (len(str1) + 1, seq3(str2[len(str1)]), 
                                        abs(len(str1) - len(str2)))
        if len(str1) > len(str2) :
            return "p.(%s%i*)" % (seq3(str1[len(str2)]), len(str2) + 1)
        return "p.(%s%i%s)" % (seq3(str1[lcp]), lcp + 1, seq3(str2[lcp]))
    #if

    # InDel.
    if lcp + 1 == str1_end :
        return "p.(%s%idelins%s)" % (seq3(str1[lcp]), lcp + 1, 
                                     seq3(str2[lcp:str2_end]))
    return "p.(%s%i_%s%idelins%s)" % (seq3(str1[lcp]), lcp + 1, 
                                      seq3(str1[str1_end - 1]),
                                      str1_end, seq3(str2[lcp:str2_end]))
#findInFrameDescription

def findFrameShift(str1, str2) :
    """
    """

    lcp = __lcp(str1, str2)

    return "p.(%s%i%sfs*%i)" % (seq3(str1[lcp]), lcp + 1, seq3(str2[lcp]),
                                len(str2) - lcp)
#findFrameShift

def __toProtDescr(CDSStop, orig, trans) :
    """
    """

    if CDSStop % 3 :
        return findFrameShift(str(orig), str(trans))
    return findInFrameDescription(str(orig), str(trans))
#__toProtDescr

def __trim(string, lcp, lcs) :
    """
    """

    return string[lcp:len(string) - lcs]
#__trim

def __rangeToC(M, g1, g2) :
    """
    """

    if M.orientation == -1 :
        return M.g2c(g2), M.g2c(g1)
    return M.g2c(g1), M.g2c(g2)
#__rangeToC

def __maybeInvert(M, string) :
    """
    """

    if M.orientation == -1 :
        return Bio.Seq.reverse_complement(string)
    return string
#__maybeInvert

def __searchFrameShift(orig, mutated) :
    pass
#__searchFrameShift

def __rv(MUU, record, RawVar, GenRecordInstance, RefType, O, NM) :
    """
    """

    for i in GenRecordInstance.record.geneList :
        for j in i.transcriptList :
            M = j.CM
            if j.CM :
                start_main = j.CM.main2int(RawVar.StartLoc.PtLoc.MainSgn + 
                                        RawVar.StartLoc.PtLoc.Main)
                start_offset = __PtLoc2offset(RawVar.StartLoc.PtLoc)
                end_main = start_main
                end_offset = start_offset
                if RawVar.EndLoc :
                    end_main = j.CM.main2int(RawVar.EndLoc.PtLoc.MainSgn + 
                                          RawVar.EndLoc.PtLoc.Main)
                    end_offset = __PtLoc2offset(RawVar.EndLoc.PtLoc)
                #if
                
                start_g = int(start_main)
                end_g = int(end_main)

                Arg1 = RawVar.Arg1
                Arg2 = RawVar.Arg2
                if RefType in ['c', 'n'] :
                    start_g = j.CM.x2g(start_main, start_offset)
                    end_g = j.CM.x2g(end_main, end_offset)
                    if j.CM.orientation == -1 :
                        Arg1 = Bio.Seq.reverse_complement(RawVar.Arg1)
                        Arg2 = Bio.Seq.reverse_complement(RawVar.Arg2)
                    #if
                #if
                
                start_g, end_g = __order(start_g, end_g)

                start_c, end_c = __rangeToC(M, start_g, end_g)
                if start_c.isdigit() and 1 <= int(start_c) <= 3 or \
                   end_c.isdigit() and 1 <= int(end_c) <= 3 :
                    O.addMessage(__file__, 2, "WSTART", 
                                 "Mutation in start codon.")

                # start_offset has to be calculated (not accepted from the 
                # parser)
                start_t_m, start_t_o = j.CM.g2x(start_g)
                t_s, t_e, c_s = j.CM.info()
                if start_t_o == -1 or start_t_o == -2 :
                    if start_t_m == t_s :
                        O.addMessage(__file__, 2, "WTXSTART", 
                            "Mutation hits transcription start.")
                    else :
                        O.addMessage(__file__, 2, "WSPLDON", 
                            "Mutation hits a splice donor site.")
                #if                            
                if start_t_o == 1 or start_t_o == 2 :
                    O.addMessage(__file__, 2, "WSPLACC", 
                        "Mutation hits a splice acceptor site.")
                
                if RawVar.MutationType in ["del", "dup", "subst", "delins"] :
                    __checkOptArg(record.seq, start_g, end_g, Arg1, M, O)
            #if                    
        #for
    #for

    # Substitution.
    if RawVar.MutationType == "subst" :
        if RawVar.Arg1 == RawVar.Arg2 :
            O.addMessage(__file__, 3, "ENOVAR", 
                "No mutation given (%c>%c) at position c.%s (g.%i)." % (
                RawVar.Arg1, RawVar.Arg1, M.g2c(start_g), start_g))

        MUU.subM(start_g, Arg2)

        GenRecordInstance.name(start_g, 0, "subst", record.seq[start_g - 1], 
                               Arg2)
        NM.g += str(start_g) + record.seq[start_g - 1] + '>' + Arg2
    #if
    
    # Deletion / Duplication.
    if RawVar.MutationType in ["del", "dup"] :
        rollposstart = __roll(record.seq, start_g - 1, end_g,
                              M.orientation)
        if rollposstart != start_g :
            rollposend = rollposstart + (end_g - start_g)
            O.addMessage(__file__, 2, "WROLL", 
                "Sequence %s at position c.%s (g.%i) was given, however, " \
                "the HGVS notation prescribes that it should be %s at " \
                "position c.%s (g.%i)." % (str(record.seq[start_g - 1:end_g]), 
                M.g2c(start_g), start_g, 
                str(record.seq[rollposstart - 1:rollposend]), 
                M.g2c(rollposstart), rollposstart))
            start_g = rollposstart
            end_g = rollposend
        #if
        if RawVar.MutationType == "del" :
            MUU.delM(start_g, end_g)
        else :
            MUU.dupM(start_g, end_g)

        if start_g != end_g :
            c1, c2 = __rangeToC(M, start_g, end_g)
            cpos = c1 + '_' + c2
            gpos = str(start_g) + '_' + str(end_g)
        #if
        else :
            cpos = str(M.g2c(start_g))
            gpos = str(start_g)
        #else
        GenRecordInstance.name(start_g, end_g, RawVar.MutationType, "", "")
        NM.g += gpos + RawVar.MutationType
    #if
    
    # Inversion.
    if RawVar.MutationType == "inv" :
        snoop = __palinsnoop(record.seq[start_g - 1:end_g])
        if snoop :
            if snoop == -1 :
                O.addMessage(__file__, 2, "WNOCHANGE", 
                    "Sequence %s at position c.%s (g.%i) is a palindrome " \
                    "(its own reverse complement)." % (
                    str(record.seq[start_g - 1:end_g]), 
                    M.g2c(start_g), start_g))
                return NM
            else :
                O.addMessage(__file__, 2, "WNOTMINIMAL", 
                    "Sequence %s at position c.%s (g.%i) is a partial " \
                    "palindrome (the first %i nucleotide(s) are the reverse " \
                    "complement of the last one(s)), the HGVS notation " \
                    "prescribes that it should be %s at position c.%s " \
                    "(g.%i)." % (str(record.seq[start_g - 1:end_g]), 
                    M.g2c(start_g), start_g, snoop, 
                    str(record.seq[start_g + snoop - 1: end_g - snoop]),
                    M.g2c(start_g + snoop), start_g + snoop))
                start_g += snoop
                end_g -= snoop
        #if
        MUU.invM(start_g, end_g)
            
        c1, c2 = __rangeToC(M, start_g, end_g)
        #NM.c += c1 + '_' + c2 + "inv"
        GenRecordInstance.name(start_g, end_g, "inv", "", "")
        NM.g += str(start_g) + '_' + str(end_g) + "inv"
    #if
    
    # Insertion.
    if RawVar.MutationType == "ins" :
        if start_g + 1 != end_g :
            O.addMessage(__file__, 3, "EINSRANGE", 
                "c.%s (g.%i) and c.%s (g.%i) are not consecutive positions." \
                % (M.g2c(start_g), start_g, M.g2c(end_g), end_g))
    
        MUU.insM(start_g, Arg1)
    
        way = M.orientation
        l = len(Arg1)
        rs1 = MUU.shiftpos(start_g)
        re1 = MUU.shiftpos(start_g) + l

        rs2 = __roll(MUU.mutated, rs1, re1, way) - 1

        shiftlen = (((rs2 - rs1) * way) - l + 1) * way

        c1 = rs2
        c2 = rs2 + ((l - 1) * way)
        c1, c2 = __order(c1, c2)
        corr = 0
        if rs1 != rs2 or \
           str(MUU.mutated[c1:c2 + 1]) == str(MUU.mutated[c1-l:(c2-l)+1]) :
            O.addMessage(__file__, 2, "WINSDUP", 
                "Insertion of %s at position c.%s (g.%i) was given, " \
                "however, the HGVS notation prescribes that it should be a " \
                "duplication of %s at position c.%s (g.%i)." % (RawVar.Arg1, 
                M.g2c(start_g), start_g, str(MUU.mutated[c1:c2 + 1]),
                M.g2c(start_g + shiftlen), start_g + shiftlen))
            start_g += shiftlen
            end_g += shiftlen
            corr = 1
        #if

        c1, c2 = __rangeToC(M, start_g, end_g)
        if corr :
            #NM.c += c1 + '_' + c2 + "dup"
            GenRecordInstance.name(start_g, end_g, "dup", "", "")
            NM.g += str(start_g) + '_' + str(end_g) + "dup"
        else :
            #NM.c += c1 + '_' + c2 + "ins" + RawVar.Arg1
            GenRecordInstance.name(start_g, end_g, "ins", "", "")
            NM.g += str(start_g) + '_' + str(end_g) + "ins" + Arg1
    #if

    # DelIns.
    if RawVar.MutationType == "delins" :
        Arg1 = RawVar.Arg1
        if not RawVar.Arg1 :
            Arg1 = MUU.orig[start_g - 1:end_g]

        lcp =  __lcp(Arg1, RawVar.Arg2)
        lcs =  __lcs(Arg1, RawVar.Arg2)

        if str(Arg1) == str(RawVar.Arg2) :
            O.addMessage(__file__, 2, "WNOCHANGE", 
                "Sequence %s at position c.%s (g.%i) is identical to the " \
                "variant." % (
                str(record.seq[start_g - 1:end_g]), 
                M.g2c(start_g), start_g))
            return NM
        ins_part = RawVar.Arg2
        if lcp or lcs :
            del_part = __trim(Arg1, lcp, lcs)
            ins_part = __trim(RawVar.Arg2, lcp, lcs)
            #O.addMessage(__file__, 2, "WNOTMINIMAL", 
            #    "")

        start_g += lcp
        end_g -= lcs
        MUU.delinsM(start_g, end_g, ins_part)

        if start_g != end_g :
            c1, c2 = __rangeToC(M, start_g, end_g)
            cpos = c1 + '_' + c2
            gpos = str(start_g) + '_' + str(end_g)
        #if
        else :
            cpos = str(M.g2c(start_g))
            gpos = str(start_g)
        #else
        #NM.c += cpos + "delins" + __maybeInvert(M, ins_part)
        GenRecordInstance.name(start_g, end_g, "delins", ins_part, "")
        NM.g += gpos + "delins" + ins_part
    #if
    
    return NM
#__rv

def __ppp(MUU, record, parts, GenRecordInstance, refseq, depth, O) :
    if parts.RefSeqAcc :
        refseq = parts.RefSeqAcc

    if parts.RawVar or parts.SingleAlleleVarSet :
        NM = newMut()

        print GenRecordInstance.record.mol_type
        #if parts.RefType == 'n' :
        #    #NM.c += "n."
        #else :
        #    #NM.c += "c."
        if GenRecordInstance.record.organelle and \
           GenRecordInstance.record.organelle == "mitochondrion" :
            NM.g += "m."
        else :
            if GenRecordInstance.record.geneList :
                NM.g += "g."
            #else : # EST
            #    NM.g += ""
        #else   

        GS = GenRecordInstance.record.geneList[0].name
        if parts.SingleAlleleVarSet :
            #NM.c += '['
            NM.g += '['
            for i in parts.SingleAlleleVarSet :
                __rv(MUU, record, i.RawVar, GenRecordInstance, 
                     parts.RefType, O, NM)
                #NM.c += ';'
                NM.g += ';'
            #for
            #NM.c = NM.c[0:-1] + ']'
            NM.g = NM.g[0:-1] + ']'
        #if
        else :
            NM = __rv(MUU, record, parts.RawVar, GenRecordInstance, 
                      parts.RefType, O, NM)

        #O.addOutput("variantdescription", NM.c)
        O.addOutput("variantdescription", NM.g)
        del NM

        W = GenRecordInstance.record.geneList[0].transcriptList[0]
        if not W.CDS : # Noncoding.
            return
        #if noTrans :
        #    return
        if not GenRecordInstance.record.geneList : # EST
            return

        cds = Seq(str(__splice(MUU.orig, W.CDS.positionList)), 
                  IUPAC.unambiguous_dna)
        cdsm = Seq(str(__nsplice(MUU.mutated, 
                                 MUU.newSplice(W.mRNA.positionList), 
                                 MUU.newSplice(W.CDS.location), 
                                 W.CM.orientation)),
                   IUPAC.unambiguous_dna)
        if W.CM.orientation == -1 :
            cds = Bio.Seq.reverse_complement(cds)
            cdsm = Bio.Seq.reverse_complement(cdsm)

        if '*' in cds.translate()[:-1] :
            O.addMessage(__file__, 3, "ESTOP", "In frame stop codon found.")
            return
        #if
        orig = cds.translate(table = W.txTable, cds = True, to_stop = True)
        O.addOutput("oldprotein", orig + '*')
        trans = cdsm.translate(table = W.txTable, to_stop = True)

        if not trans or trans[0] != 'M' :
            if str(cdsm[0:3]) in \
                Bio.Data.CodonTable.unambiguous_dna_by_id[
                    W.txTable].start_codons :
                O.addOutput("newprotein", '?')
                O.addOutput("altstart", str(cdsm[0:3]))
                O.addOutput("altprotein", 'M' + trans[1:] + '*')
            else :
                __bprint('?')
                O.addOutput("newprotein", '?')
        else :
            O.addOutput("newprotein", trans + '*')

        if not parts.SingleAlleleVarSet :
            O.addOutput("proteindescription", __toProtDescr(
                W.CM.g2x(MUU.newSplice(W.CDS.location)[1])[0], orig, trans))
        else :
            O.addOutput("proteindescription", "p.?")
            
        del W.CM
    #if                
#__ppp

def process(cmd, C, O) :
    parser = Parser.Nomenclatureparser(O)
    O.addOutput("inputvariant", cmd)
    ParseObj = parser.parse(cmd)
    del parser

    if ParseObj :
        if ParseObj.Version :
            RetrieveRecord = ParseObj.RefSeqAcc + '.' + ParseObj.Version
        else :
            RetrieveRecord = ParseObj.RefSeqAcc
        O.addOutput("reference", RetrieveRecord)
        
        D = Db.Cache(C.Db)
        retriever = Retriever.Retriever(C.Retriever, O, D)
        record = retriever.loadrecord(RetrieveRecord)
        if not record :
            return
        del retriever
        del D
        
        GenRecordInstance = GenRecord.GenRecord(C.GenRecord, O)
        GenRecordInstance.parseRecord(record)

        MUU = Mutator.Mutator(record.seq, C.Mutator, O)
        __ppp(MUU, record, ParseObj, GenRecordInstance, "",  0, O)
        del MUU
        return GenRecordInstance
    #if
#process

def main(cmd) :
    C = Config.Config()
    O = Output.Output(__file__, C.Output)

    O.addMessage(__file__, -1, "INFO", "Received variant " + cmd)

    RD = process(cmd, C, O)

    O.addMessage(__file__, -1, "INFO", "Finished processing variant " + cmd)

    ### OUTPUT BLOCK ###
    gn = O.getOutput("genename")
    if gn :
        print "Gene Name: " + gn[0]
    tv = O.getOutput("transcriptvariant")
    if tv :
        print "Transcript variant: " + tv[0]
        print
    #if
    
    O.getMessages()
    errors, warnings, summary = O.Summary()
    print summary
    print

    if not errors :
        visualisation = O.getOutput("visualisation")
        if visualisation :
            for i in range(len(visualisation)) :
                if i and not i % 3 :
                    print
                print visualisation[i]
            #for
            print
        #if

        reference = O.getOutput("reference")[0]

        for i in RD.record.geneList :
            for j in i.transcriptList :
                if ';' in j.description :
                    print "%s(%s_%s):%c.[%s]" % (reference, i.name, j.name, 
                                                 j.molType, j.description)
                else :
                    print "%s(%s_%s):%c.%s" % (reference, i.name, j.name, 
                                               j.molType, j.description)

        vd = O.getOutput("variantdescription")
        if vd :
            for i in vd :
                print "%s:%s" % (reference, i)

        pd = O.getOutput("proteindescription")
        if pd :
            if O.getOutput("altprotein") :
                print "%s:p.(0)" % reference
            else :
                print "%s:%s" % (reference, pd[0])
        #if

        op = O.getOutput("oldprotein")
        if op :
            print "\n<b>Old protein:</b>"
            __bprint(op[0])
            print
        #if
        np = O.getOutput("newprotein")
        if np :
            print "\n<b>New protein:</b>"
            __bprint(np[0])
            print
        #if
        ap = O.getOutput("altprotein")
        if ap :
            print "\n<b>Alternative protein using start codon %s:</b>" % \
                O.getOutput("altstart")[0]
            __bprint(ap[0])
            print
        #if
    #if
    ### OUTPUT BLOCK ###
    del O
#main

if __name__ == "__main__" :
    main(sys.argv[1])
#if
