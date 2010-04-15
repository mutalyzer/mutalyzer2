#!/usr/bin/python

from Modules import Retriever
from Modules import GenRecord
from Modules import Crossmap
from Modules import Parser

import types
from Modules import Output

import Bio.Seq

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

def __palinsnoop(string) :
    import math

    revcomp = Bio.Seq.reverse_complement(string)

    for i in range(int(math.ceil(len(string) / 2.0))) :
        if string[i] != revcomp[i] :
            return i # The first i elements are palindromic.
    return -1        # Perfect palindrome.
#__palinsnoop

def __bprint(s) :
    import math

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
    main = int(Loc.Main)
    if Loc.MainSgn == '-' :
        main = -main

    return main
#__PtLoc2main

def __PtLoc2offset(Loc) :
    if Loc.Offset :
        offset = int(Loc.Offset)
        if Loc.OffSgn == '-' :
            offset = -offset
    #if
    else :
        offset = 0

    return offset
#__PtLoc2offset

"""
def IsInt(string) :
    try :
        num = int(string)
        return 1
    #try
    except ValueError :
        return 0
#IsInt
"""

"""
def printp(string, depth) :
    print (depth * "  ") + str(string)
"""

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
    transcript = ""

    if orientation == 1 :
        for i in range(0, len(splice_sites), 2) :
            if CDS[0] >= splice_sites[i] and CDS[0] <= splice_sites[i + 1] :
                transcript += string[CDS[0] - 1:splice_sites[i + 1]]
            else :
                if splice_sites[i] > CDS[0] :
                    transcript += string[splice_sites[i] - 1:splice_sites[i + 1]] 
    else :
        for i in range(0, len(splice_sites), 2) :
            if CDS[1] >= splice_sites[i] and CDS[1] <= splice_sites[i + 1] :
                transcript += string[splice_sites[i] - 1: CDS[1]]
            else :
                if splice_sites[i] < CDS[1] :
                    transcript += string[splice_sites[i] - 1:splice_sites[i + 1]] 


    return transcript
#__nsplice


def __toProtDescr(orig, trans) :
    from Bio.SeqUtils import seq3

    if str(trans) == str(orig) :
        print "p.="
        return

    if len(trans) > len(orig) :
        ext = abs(len(orig) - len(trans))
        print "p.*%i%sext*%i" % (len(orig) + 1, seq3(trans[len(orig)]), ext)
        return

    if len(orig) > len(trans) :
        print "p.%s%i*" % (seq3(orig[len(trans)]), len(trans) + 1)
        #print "p.%s%i*" % (seq3(orig[len(trans) - 1]), i + 1)
        return

    i = 0
    while i < len(orig) - 1 and orig[i] == trans[i] :
        i += 1

    print "p.%s%i%s" % (seq3(orig[i]), i + 1, seq3(trans[i]))
#__toProtDescr

def __constructCDS(mRNA, CDSpos) :
    #print mRNA
    #print CDSpos
    i = 1
    ret = [CDSpos[0]]

    while CDSpos[0] > mRNA[i] :
        i += 2

    j = i
    while CDSpos[1] > mRNA[j] :
        j += 2

    ret.extend(mRNA[i:j])
    ret.append(CDSpos[1])

    #print ret
    return ret
#__constructCDS

"""
def __isStringThere(ref, p1, p2, string) :
    if ref[p1 - 1:p2] == string :
        return True
    return False
#__isStringThere

def __checkStringLength(p1, p2, length) :
    if p2 - p1 + 1 == int(length) :
        return True
    return False
#__checkStringLength
"""

def __checkOptArg(ref, p1, p2, arg, M, O) :
    if arg :
        if arg.isdigit() :
            length = int(arg)
            interval = p2 - p1 + 1
            if length != interval :
                O.ErrorMsg(__file__, "The length (%i) differed from that of " \
                    "the range (%i)." % (length, interval))
                return False
            #if
        #if
        else :
            #revcomp = reverse_complement(string)

            ref_slice = str(ref[p1 - 1:p2])
            #if M.orientation == -1 :
            #    ref_slice = Bio.Seq.reverse_complement(ref_slice)
            if ref_slice != str(arg) : # FIXME more informative.
                O.ErrorMsg(__file__, "%s not found at position c.%s (g.%i), " \
                    "found %s instead." % (arg, M.g2c(p1), p1, ref_slice))
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

def __trim(string, lcp, lcs) :
    if lcp and lcs :
        return string[lcp:-lcs]
    if lcp :
        return string[lcp:]
    if lcs :
        return string[:-lcs]
    return string
#__trim

def __rangeToC(M, g1, g2) :
    if M.orientation == -1 :
        return M.g2c(g2), M.g2c(g1)
    return M.g2c(g1), M.g2c(g2)
#__rangeToC

def __maybeInvert(M, string) :
    if M.orientation == -1 :
        return Bio.Seq.reverse_complement(string)
    return string
#__maybeInvert

def __rv(MUU, record, GeneSymbol, RawVar, M, RefType, O, NM) :
    #start_main = __PtLoc2main(RawVar.StartLoc.PtLoc)
    start_main = M.main2int(RawVar.StartLoc.PtLoc.MainSgn + 
                            RawVar.StartLoc.PtLoc.Main)
    start_offset = __PtLoc2offset(RawVar.StartLoc.PtLoc)
    end_main = start_main
    end_offset = start_offset
    if RawVar.EndLoc :
        #end_main = __PtLoc2main(RawVar.EndLoc.PtLoc)
        end_main = M.main2int(RawVar.EndLoc.PtLoc.MainSgn + 
                              RawVar.EndLoc.PtLoc.Main)
        end_offset = __PtLoc2offset(RawVar.EndLoc.PtLoc)
    #if
    
    start_g = int(start_main)
    end_g = int(end_main)

    Arg1 = RawVar.Arg1
    Arg2 = RawVar.Arg2
    if RefType in ['c', 'n'] :
        start_g = M.x2g(start_main, start_offset)
        end_g = M.x2g(end_main, end_offset)
        if M.orientation == -1 :
            Arg1 = Bio.Seq.reverse_complement(RawVar.Arg1)
            Arg2 = Bio.Seq.reverse_complement(RawVar.Arg2)
        #if
    #if
    
    start_g, end_g = __order(start_g, end_g)

    start_c, end_c = __rangeToC(M, start_g, end_g)
    if start_c.isdigit() and 1 <= int(start_c) <= 3 or \
       end_c.isdigit() and 1 <= int(end_c) <= 3 :
        O.WarningMsg(__file__, "Mutation in start codon.")

    # start_offset has to be calculated (not accepted from the parser)
    start_t_m, start_t_o = M.g2x(start_g)
    t_s, t_e, c_s = M.info()
    if start_t_o == -1 or start_t_o == -2 :
        if start_t_m == t_s :
            O.WarningMsg(__file__, "Mutation hits transcription start.")
        else :
            O.WarningMsg(__file__, "Mutation hits a splice donor site.")
    if start_t_o == 1 or start_t_o == 2 :
        O.WarningMsg(__file__, "Mutation hits a splice acceptor site.")
    
    #print str(record.seq[start_g - 20:start_g + 20])
    
    if RawVar.MutationType in ["del", "dup", "subst", "delins"] :
        __checkOptArg(record.seq, start_g, end_g, Arg1, M, O)

    global protDescr
    protDescr = False
    # Substitution.
    if RawVar.MutationType == "subst" :
        if RawVar.Arg1 == RawVar.Arg2 :
            O.ErrorMsg(__file__, "No mutation given (%c>%c) at position " \
                "c.%s (g.%i)." % (RawVar.Arg1, RawVar.Arg1, M.g2c(start_g), 
                start_g))
        MUU.subM(start_g, Arg2)

        NM.c += str(M.g2c(start_g)) + \
                    __maybeInvert(M, record.seq[start_g - 1]) + \
                    '>' + __maybeInvert(M, Arg2)
        NM.g += str(start_g) + \
                   record.seq[start_g - 1] + \
                   '>' + Arg2
        protDescr = True                   
    #if
    
    # Deletion / Duplication.
    if RawVar.MutationType in ["del", "dup"] :
        rollposstart = __roll(record.seq, start_g - 1, end_g,
                              M.orientation)
        if rollposstart != start_g :
            rollposend = rollposstart + (end_g - start_g)
            O.WarningMsg(__file__, "Sequence %s at position c.%s (g.%i) was " \
                "given, however, the HGVS notation prescribes that it should " \
                "be %s at position c.%s (g.%i)." % (
                str(record.seq[start_g - 1:end_g]), M.g2c(start_g), start_g,
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
        NM.c += cpos + RawVar.MutationType
        NM.g += gpos + RawVar.MutationType
    #if
    
    # Inversion.
    if RawVar.MutationType == "inv" :
        snoop = __palinsnoop(record.seq[start_g - 1:end_g])
        if snoop :
            if snoop == -1 :
                O.WarningMsg(__file__, "Sequence %s at position c.%s (g.%i) " \
                    "is a palindrome (its own reverse complement)." % (
                    str(record.seq[start_g - 1:end_g]), 
                    M.g2c(start_g), start_g))
                # Do nothing...
            else :
                O.WarningMsg(__file__, "Sequence %s at position c.%s (g.%i) " \
                    "is a partial palindrome (the first %i " \
                    "nucleotide(s) are the reverse complement of " \
                    "the last one(s)), the HGVS notation " \
                    "prescribes that it should be %s at position " \
                    "c.%s (g.%i)." % (
                    str(record.seq[start_g - 1:end_g]), M.g2c(start_g),
                    start_g, snoop, 
                    str(record.seq[start_g + snoop - 1: end_g - snoop]),
                    M.g2c(start_g + snoop), start_g + snoop))
                start_g += snoop
                end_g -= snoop
        #if
        MUU.invM(start_g, end_g)
            
        c1, c2 = __rangeToC(M, start_g, end_g)
        NM.c += c1 + '_' + c2 + "inv"
        NM.g += str(start_g) + '_' + str(end_g) + "inv"
    #if
    
    # Insertion.
    if RawVar.MutationType == "ins" :
        if start_g + 1 != end_g :
            O.ErrorMsg(__file__, "c.%s (g.%i) and c.%s (g.%i) are not " \
                "consecutive positions." % (
                M.g2c(start_g), start_g, M.g2c(end_g), end_g))
    
        #inserted = RawVar.Arg1
        #if M.orientation == -1 :
        #    inserted = Bio.Seq.reverse_complement(RawVar.Arg1)

        MUU.insM(start_g, Arg1)
    
        way = M.orientation
        l = len(Arg1)
        rs1 = MUU.shiftpos(start_g)
        re1 = MUU.shiftpos(start_g) + l

        #rs1, re1 = __order(rs1, re1)
        #print "+++", MUU.mutated[rs1:re1], rs1, re1

        rs2 = __roll(MUU.mutated, rs1, re1, way) - 1

        #shiftlen = ((rs2 - rs1 - l + 1) * way)
        shiftlen = (((rs2 - rs1) * way) - l + 1) * way
        #print "+++", rs2, shiftlen

        c1 = rs2
        c2 = rs2 + ((l - 1) * way)
        c1, c2 = __order(c1, c2)
        corr = 0
        if rs1 != rs2 or \
           str(MUU.mutated[c1:c2 + 1]) == str(MUU.mutated[c1-l:(c2-l)+1]) :
            O.WarningMsg(__file__, "Insertion of %s at position c.%s (g.%i) " \
                "was given, however, the HGVS notation prescribes that it " \
                "should be a duplication of %s at position c.%s (g.%i)." % (
                RawVar.Arg1, M.g2c(start_g), start_g, 
                str(MUU.mutated[c1:c2 + 1]), M.g2c(start_g + shiftlen), 
                start_g + shiftlen))
            start_g += shiftlen
            end_g += shiftlen
            corr = 1
        #if

        c1, c2 = __rangeToC(M, start_g, end_g)
        if corr :
            NM.c += c1 + '_' + c2 + "dup"
            NM.g += str(start_g) + '_' + str(end_g) + "dup"
        else :
            NM.c += c1 + '_' + c2 + "ins" + RawVar.Arg1
            NM.g += str(start_g) + '_' + str(end_g) + "ins" + Arg1
    #if

    if RawVar.MutationType == "delins" :
        lcp =  __lcp(RawVar.Arg1, RawVar.Arg2)
        lcs =  -__lcs(RawVar.Arg1, RawVar.Arg2)

        if lcp or lcs :
            del_part = __trim(RawVar.Arg1, lcp, lcs)
            ins_part = __trim(RawVar.Arg2, lcp, lcs)
            print start_g + lcp, end_g - lcs
            print M.g2c(start_g + lcp), M.g2c(end_g - lcs)
            print "del%sins%s" % (del_part, ins_part)

        MUU.delinsM(start_g, end_g, Arg2)

        #NM.c += str(M.g2c(start_g)) + '_' + str(M.g2c(end_g)) + "delins" + RawVar.Arg2
        #NM.g += str(start_g) + '_' + str(end_g) + "delins" + RawVar.Arg2

        if start_g != end_g :
            c1, c2 = __rangeToC(M, start_g, end_g)
            cpos = c1 + '_' + c2
            gpos = str(start_g) + '_' + str(end_g)
        #if
        else :
            cpos = str(M.g2c(start_g))
            gpos = str(start_g)
        #else
        #NM.c += cpos + "delins" + RawVar.Arg2 
        #NM.g += gpos + "delins" + RawVar.Arg2
        NM.c += cpos + "delins" + __maybeInvert(M, Arg2)
        NM.g += gpos + "delins" + Arg2
    #if
    
    #print MUU.mutated[start_g - 20:start_g + 20]
    return NM
#__rv

def __ppp(MUU, record, parts, recordDict, refseq, depth, O) :
    #printp("+++ recurse +++", depth)
    #printp(repr(parts), depth)
    #printp(parts, depth)
    """
    printp(parts[4], depth)
    printp(repr(parts[4]), depth)
    printp(parts[4][0][0][0].Main, 2)
    printp(parts, depth)
    """
    if parts.RefSeqAcc :
        refseq = parts.RefSeqAcc
    #printp("RefSeqAcc: " + str(refseq), depth)
    #printp("RefType: " + str(parts.RefType), depth)


    #printp("Version: " + str(parts.Version), depth)
    if parts.Gene :
        print "Gene Symbol: " + str(parts.Gene.GeneSymbol)
        #printp("Transcript variant: " + str(parts.Gene.TransVar), depth)
        #printp("Protein isoform: " + str(parts.Gene.ProtIso), depth)
    #if
    """
    if parts.ChimeronSet :
        #printp(str(parts.ChimeronSet), depth)
        for i in parts.ChimeronSet :
            printp("ChimeronSet", depth)
            __ppp(MUU, record, i, recordDict, refseq, depth + 1, O)
        #for
    #if
    if parts.MosaicSet :
        #printp(str(parts.MosaicSet), depth)
        for i in parts.MosaicSet :
            printp("MosaicSet", depth)
            __ppp(MUU, record, i, recordDict, refseq, depth + 1, O)
        #for
    #if
    if parts.SimpleAlleleVarSet :
        #printp(str(parts.SimpleAlleleVarSet), depth)
        for i in parts.SimpleAlleleVarSet :
            printp("SimpleAlleleVarSet", depth)
            __ppp(MUU, record, i, recordDict, refseq, depth + 1, O)
        #for
            #__ppp(MUU, record, i, recordDict, refseq, depth + 1)
    #if
    if parts.MultiAlleleVars :
        printp(str(parts.MultiAlleleVars), depth)
        for i in parts.MultiAlleleVars :
            printp("MultiAlleleVars", depth)
            print repr(i)
            __ppp(MUU, record, i, recordDict, refseq, depth + 1, O)
        #for
    #if
    """
    if parts.RawVar or parts.SingleAlleleVarSet :
        GS = ""
        if recordDict.genelist :
            GS = recordDict.genelist.keys()[0]

        if parts.Gene and parts.Gene.GeneSymbol :
            GS = parts.Gene.GeneSymbol
        print "Gene Name: " + GS

        transcriptvariant = "001" 
        if parts.Gene and parts.Gene.TransVar :
            transcriptvariant = parts.Gene.TransVar
        print "Transcript variant: " + transcriptvariant
        print
            
        if recordDict.genelist :
            if recordDict.genelist.has_key(GS) :
                currentGene = recordDict.genelist[GS]
            else :
                print "No such gene %s in record." % GS
                return
            #else
        else :
            currentGene = recordDict.source
        W = currentGene.list[transcriptvariant]

        noTrans = False
        if not W.mRNA :
            if not W.exon:
                O.WarningMsg(__file__, "No mRNA field found for gene %s, " \
                    "transcript variant %s in GenBank record %s, " \
                    "constructing it from CDS." % (GS, transcriptvariant, 
                    record.id))
                if W.CDS :
                    if not W.CDS.list :
                        print "Extra warning"
                        W.mRNA = W.CDS
                        W.mRNA.list = W.CDS.location
                        noTrans = True
                    else :
                        W.mRNA = W.CDS
                #if
                else :
                    print currentGene.location
                    W.CDS = GenRecord.Locus()
                    W.CDS.location = W.location
                    W.mRNA = W.CDS
                    W.mRNA.list = currentGene.location
                    noTrans = True
            #if
            else :
                O.WarningMsg(__file__, "No mRNA field found for gene %s, " \
                    "transcript variant %s in GenBank record %s, " \
                    "constructing it from gathered exon information." % (GS, 
                    transcriptvariant, record.id))
                W.mRNA = W.exon
        #if
        #print W.mRNA.list
        if not W.mRNA.list :
            W.mRNA.list = W.mRNA.location
        if W.CDS :
            if not W.CDS.list :
                O.WarningMsg(__file__, "No CDS list found for gene %s, " \
                    "transcript variant %s in GenBank record %s, " \
                    "constructing it from mRNA list and CDS location." % (GS, 
                    transcriptvariant, record.id))
                if W.mRNA.list :
                    W.CDS.list = __constructCDS(W.mRNA.list, W.CDS.location)
                    #print W.mRNA.list, W.CDS.location
                    #print W.CDS.list
                else :
                    W.CDS.list = __constructCDS(W.mRNA.location, W.CDS.location)
        #else :
        #    pass # Noncoding RNA?

        if parts.RefType == 'n' :
            M = Crossmap.Crossmap(
                W.mRNA.list,
                [],
                currentGene.orientation)
        else :                
            if not W.CDS :
                O.ErrorMsg(__file__, "No CDS information found for gene %s, " \
                                     "transcript variant %s in GenBank " \
                                     "record %s." % (GS, transcriptvariant, 
                                                     record.id))
                return
            #if
            M = Crossmap.Crossmap(
                W.mRNA.list,
                W.CDS.location,
                currentGene.orientation)
        #else
        #print W.mRNA

        #print recordDict["organelle"], recordDict["mol_type"]
        NM = newMut()
        NM.c = parts.RefSeqAcc 
        NM.g = parts.RefSeqAcc 
        if parts.Version :
            NM.c += '.' + parts.Version
            NM.g += '.' + parts.Version
        #if
        if parts.RefType == 'n' :
            NM.c += ":n."
        else :
            NM.c += ":c."
        if recordDict.organelle and recordDict.organelle == "mitochondrion" :
            NM.g += ":m."
        else :
            if recordDict.genelist :
                NM.g += ":g."
            else : # EST
                NM.g += ":"
        #else   

        if parts.SingleAlleleVarSet :
            NM.c += '['
            NM.g += '['
            for i in parts.SingleAlleleVarSet :
                __rv(MUU, record, GS, i.RawVar, M, parts.RefType, O, NM)
                NM.c += ';'
                NM.g += ';'
            #for
            NM.c = NM.c[0:-1] + ']'
            NM.g = NM.g[0:-1] + ']'
        #if
        else :
            NM = __rv(MUU, record, GS, parts.RawVar, M, parts.RefType, O, NM)

        print
        print "+++", NM.c
        print "+++", NM.g
        del NM

        if not W.CDS : # Noncoding.
            return
        if noTrans :
            return
        if not recordDict.genelist : # EST
            return

        import Bio
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC

        cds = Seq(str(__splice(MUU.orig, W.CDS.list)), IUPAC.unambiguous_dna)
        cdsm = Seq(str(__nsplice(MUU.mutated, MUU.newSplice(W.mRNA.list), 
                                 MUU.newSplice(W.CDS.location), M.orientation)),
                   IUPAC.unambiguous_dna)
        if M.orientation == -1 :
            cds = Bio.Seq.reverse_complement(cds)
            cdsm = Bio.Seq.reverse_complement(cdsm)
        del M

        print "\n<b>Old protein:</b>"
        if '*' in cds.translate()[:-1] :
            print "In frame stop codon found."
            return
        #if
        orig = cds.translate(table = W.txTable, cds = True, to_stop = True)
        __bprint(orig + '*')
        print "\n\n<b>New protein:</b>"
        trans = cdsm.translate(table = W.txTable, to_stop = True)

        if not trans or trans[0] != 'M' :
            if str(cdsm[0:3]) in \
                Bio.Data.CodonTable.unambiguous_dna_by_id[
                    W.txTable].start_codons :
                __bprint('?')
                print "\n\n<b>Alternative protein using start codon %s:</b>" % \
                    str(cdsm[0:3])
                __bprint('M' + trans[1:] + '*')
            else :
                __bprint('?')
        else :
            __bprint(trans + '*')

        if protDescr :
            print
            print
            __toProtDescr(orig, trans)
        #if
    #if                
#__ppp

def main(cmd) :
    O = Output.Output(__file__)

    O.LogMsg(__file__, "Received variant " + cmd)

    parser = Parser.Nomenclatureparser(O)
    print cmd
    print
    ParseObj = parser.parse(cmd)
    del parser

    if ParseObj :
        if ParseObj.Version :
            RetrieveRecord = ParseObj.RefSeqAcc + '.' + ParseObj.Version
        else :
            RetrieveRecord = ParseObj.RefSeqAcc
        
        #print "Retrieving..."

        retriever = Retriever.Retriever()
        record = retriever.loadrecord(RetrieveRecord)
        if not record :
            return
        del retriever
        
        #print "Dicting..."
        D = GenRecord.GenRecord()
        d = D.record2dict(record)
        del D
        #print "Printing..."
        #D.printRecordDict(d, record)
        
        from Modules import Mutator
        
        MUU = Mutator.Mutator(record.seq)
        __ppp(MUU, record, ParseObj, d, "",  0, O)
        del MUU
    #if
    print "\n\n"
    O.Summary()
    O.LogMsg(__file__, "Finished processing variant " + cmd)
    del O
#main

if __name__ == "__main__" :
    import sys

    main(sys.argv[1])
#if
