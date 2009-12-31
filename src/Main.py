#!/usr/bin/python

import Retriever
import GenRecord
import Crossmap
import Parser

import types
import Config
import Output

def roll(string, start, stop, orientation) :
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
#roll

def palinsnoop(string) :
    from Bio.Seq import reverse_complement
    import math

    revcomp = reverse_complement(string)

    for i in range(int(math.ceil(len(string) / 2.0))) :
        if string[i] != revcomp[i] :
            return i # The first i elements are palindromic.
    return -1        # Perfect palindrome.
#palinsnoop

def bprint(s) :
    import math

    block = 10
    line = 6 * block

    m = int(math.floor(math.log(len(s), 10)) + 1)
    o = 1
    print "%s " % str(o).rjust(m),
    for i in range(0, len(s), block) :
        print s[i:i + block],
        if not (i + block) % line :
            o += line
            print "\n%s " % str(o).rjust(m),
        #if
    #for
#bprint

def PtLoc2main(Loc) :
    main = int(Loc.Main)
    if Loc.MainSgn == '-' :
        main = -main

    return main
#PtLoc2int

def PtLoc2offset(Loc) :
    if Loc.Offset :
        offset = int(Loc.Offset)
        if Loc.OffSgn == '-' :
            offset = -offset
    #if
    else :
        offset = 0

    return offset
#PtLoc2int

"""
def ErrorMsg(message) :
    print "%s: Error: %s" % (__file__.split('/')[-1].split('.')[0], message)
#ErrorMsg

def WarningMsg(message) :
    print "%s: Warning: %s" % (__file__.split('/')[-1].split('.')[0], message)
#WarningMsg
"""

def IsInt(string) :
    try :
        num = int(string)
        return 1
    #try
    except ValueError :
        return 0
#IsInt

def printp(string, depth) :
    print (depth * "  ") + str(string)

def constructCDS(mRNA, CDSpos) :
    print mRNA
    print CDSpos
    i = 1
    ret = [CDSpos[0]]

    while CDSpos[0] > mRNA[i] :
        i += 2

    j = i
    while CDSpos[1] > mRNA[j] :
        j += 2

    ret.extend(mRNA[i:j])
    ret.append(CDSpos[1])

    print ret
    return ret
#constructCDS

def splice(str, splice_sites) :
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

    transcript = ''

    for i in range(0, len(splice_sites), 2) :
        transcript += str[splice_sites[i] - 1:splice_sites[i + 1]] 

    return transcript
#splice

def rv(MUU, record, recordDict, GeneSymbol, RawVar, M, RefType, O) :
    start_main = PtLoc2main(RawVar.StartLoc.PtLoc)
    start_offset = PtLoc2offset(RawVar.StartLoc.PtLoc)
    end_main = start_main
    end_offset = start_offset
    if RawVar.EndLoc :
        end_main = PtLoc2main(RawVar.EndLoc.PtLoc)
        end_offset = PtLoc2offset(RawVar.EndLoc.PtLoc)
    #if
    
    start_g = int(start_main)
    end_g = int(end_main)
    if RefType == 'c' :
        start_g = M.x2g(start_main, start_offset)
        end_g = M.x2g(end_main, end_offset)
    #if
    
    if recordDict[GeneSymbol].orientation == -1 :
        temp = start_g
        start_g = end_g
        end_g = temp
    #if

    start_c = M.g2c(start_g)
    end_c = M.g2c(end_g)
    if IsInt(start_c) and 1 <= int(start_c) <= 3 or \
       IsInt(end_c) and 1 <= int(end_c) <= 3 :
        O.WarningMsg(__file__, "Mutation in start codon.")

    if start_offset == -1 or start_offset == -2 :
        O.WarningMsg(__file__, "Mutation hits a splice donor site.")
    if start_offset == 1 or start_offset == 2 :
        O.WarningMsg(__file__, "Mutation hits a splice acceptor site.")
    
    #print str(record.seq[start_g - 20:start_g + 20])
    
    # Substitution.
    if RawVar.MutationType == "subst" :
        if record.seq[start_g - 1] != RawVar.Arg1 :
            O.ErrorMsg(__file__, "No nucleotide %c at position c.%s (g.%i), found a " \
                "%c instead." % (RawVar.Arg1, M.g2c(start_g), 
                start_g, record.seq[start_g - 1]))
        if RawVar.Arg1 == RawVar.Arg2 :
            O.ErrorMsg(__file__, "No mutation given (%c>%c) at position c.%s (g.%i)." % 
                (RawVar.Arg1, RawVar.Arg1, M.g2c(start_g), 
                start_g))
        MUU.subM(start_g, RawVar.Arg2)
    #if
    
    # Deletion / Duplication.
    if RawVar.MutationType == "del" or \
       RawVar.MutationType == "dup" :
        if RawVar.Arg1 and \
            RawVar.Arg1 != str(record.seq[start_g - 1:end_g]) :
            if not IsInt(RawVar.Arg1) :
                O.ErrorMsg(__file__, "String %s not found at position c.%s (g.%i), " \
                    "found %s instead." % (RawVar.Arg1, 
                    M.g2c(start_g), start_g, 
                    str(record.seq[start_g - 1:end_g])))
            else :
                if end_g - start_g + 1 != int(RawVar.Arg1) :
                    O.ErrorMsg(__file__, "The length of the deletion (%i) at " \
                        "position c.%s (g.%i) differed from that of the " \
                        "range (%i)." % (int(RawVar.Arg1), 
                        M.g2c(start_g), start_g, end_g - start_g + 1))
        #if
        rollposstart = roll(record.seq, start_g - 1, end_g,
                       recordDict[GeneSymbol].orientation)
        if rollposstart != start_g :
            rollposend = rollposstart + (end_g - start_g)
            O.WarningMsg(__file__, "String %s at position c.%s (g.%i) was given, " \
                "however, the HGVS notation prescribes that it should be " \
                "%s at position c.%s (g.%i)." % (
                str(record.seq[start_g - 1:end_g]), M.g2c(start_g), start_g,
                str(record.seq[rollposstart - 1:rollposend]), 
                M.g2c(rollposstart), rollposstart))
        #if
        if RawVar.MutationType == "del" :
            MUU.delM(start_g, end_g)
        else :
            MUU.dupM(start_g, end_g)
            print start_g, " ", end_g
            print "oioioi"
    #if
    
    # Inversion.
    if RawVar.MutationType == "inv" :
        snoop = palinsnoop(record.seq[start_g - 1:end_g])
        if snoop :
            if snoop == -1 :
                O.WarningMsg(__file__, "String %s at position c.%s (g.%i) is a " \
                    "palindrome (its own reverse complement)." % (
                    str(record.seq[start_g - 1:end_g]), 
                    M.g2c(start_g), start_g))
            else :
                O.WarningMsg(__file__, "String %s at position c.%s (g.%i) is a " \
                    "partial palindrome (the first %i " \
                    "nucleotide(s) are the reverse complement of " \
                    "the last one(s)), the HGVS notation " \
                    "prescribes that it should be %s at position " \
                    "c.%s (g.%i)." % (
                    str(record.seq[start_g - 1:end_g]), M.g2c(start_g),
                    start_g, snoop, 
                    str(record.seq[start_g + snoop - 1: end_g - snoop]),
                    M.g2c(start_g + snoop), start_g + snoop))
        #if
        MUU.invM(start_g, end_g)
    #if
    
    # Insertion.
    if RawVar.MutationType == "ins" :
        if start_g + 1 != end_g :
            O.ErrorMsg(__file__, "c.%s (g.%i) and c.%s (g.%i) are not consecutive " \
                "positions." % (M.g2c(start_g), start_g, M.g2c(end_g),
                end_g))
    
        MUU.insM(start_g, RawVar.Arg1)
    
        # Niet record.seq, maar de gemuteerde seq!
        rollposstart = roll(record.seq, start_g - 1, 
                       start_g + len(RawVar.Arg1),
                       recordDict[GeneSymbol].orientation)
        if rollposstart != start_g :
            rollposend = rollposstart + len(RawVar.Arg1)
            O.WarningMsg(__file__, "Insertion of %s at position c.%s (g.%i) was " \
                "given, however, the HGVS notation prescribes that it " \
                "should be a duplication of %s at position c.%s (g.%i)." % (
                RawVar.Arg1, M.g2c(start_g), start_g, 
                str(record.seq[rollposstart - 1:rollposend]), 
                M.g2c(rollposstart), rollposstart))
        #if
    #if
    
    #print MUU.mutated[start_g - 20:start_g + 20]
#rv

def ppp(MUU, record, parts, recordDict, refseq, depth, O) :
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
        printp("Gene Symbol: " + str(parts.Gene.GeneSymbol), depth)
        #printp("Transcript variant: " + str(parts.Gene.TransVar), depth)
        #printp("Protein isoform: " + str(parts.Gene.ProtIso), depth)
    #if
    if parts.ChimeronSet :
        #printp(str(parts.ChimeronSet), depth)
        for i in parts.ChimeronSet :
            printp("ChimeronSet", depth)
            ppp(MUU, record, i, recordDict, refseq, depth + 1, O)
        #for
    #if
    if parts.MosaicSet :
        #printp(str(parts.MosaicSet), depth)
        for i in parts.MosaicSet :
            printp("MosaicSet", depth)
            ppp(MUU, record, i, recordDict, refseq, depth + 1, O)
        #for
    #if
    if parts.SimpleAlleleVarSet :
        #printp(str(parts.SimpleAlleleVarSet), depth)
        for i in parts.SimpleAlleleVarSet :
            printp("SimpleAlleleVarSet", depth)
            ppp(MUU, record, i, recordDict, refseq, depth + 1, O)
        #for
            #ppp(MUU, record, i, recordDict, refseq, depth + 1)
    #if
    if parts.MultiAlleleVars :
        printp(str(parts.MultiAlleleVars), depth)
        for i in parts.MultiAlleleVars :
            printp("MultiAlleleVars", depth)
            print repr(i)
            ppp(MUU, record, i, recordDict, refseq, depth + 1, O)
        #for
    #if
    if parts.RawVar or parts.SingleAlleleVarSet :
        GS = recordDict.keys()[0]
        if parts.Gene and parts.Gene.GeneSymbol :
            GS = parts.Gene.GeneSymbol
        printp("Gene Name: " + GS, depth)

        transcriptvariant = "001" 
        if parts.Gene and parts.Gene.TransVar :
            transcriptvariant = parts.Gene.TransVar
        printp("Transcript variant: " + transcriptvariant, depth)
        print
            
        W = recordDict[GS].list[transcriptvariant]

        if not W.mRNA :
            if not W.exon:
                O.WarningMsg(__file__, "No mRNA field found for gene %s, " \
                    "transcript variant %s in GenBank record %s, " \
                    "constructing it from CDS." % (GS, transcriptvariant, 
                    record.id))
                W.mRNA = W.CDS
            #if
            else :
                W.mRNA = W.exon
        #if
        #print W.mRNA.list
        #if W.CDS :
        #    if not W.CDS.list :
        #        W.CDS.list = constructCDS(W.mRNA.list, W.CDS.location)
        #else :
        #    pass # Noncoding RNA?

        M = Crossmap.Crossmap(
          W.mRNA.list,
          W.CDS.location,
          recordDict[GS].orientation)
        #print W.mRNA

        if parts.SingleAlleleVarSet :
            for i in parts.SingleAlleleVarSet :
                rv(MUU, record, recordDict, GS, i.RawVar, M, parts.RefType, O)
        else :
            rv(MUU, record, recordDict, GS, parts.RawVar, M, parts.RefType, O)

        import Bio
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC

        #print W.CDS.list
        #print splice(MUU.orig, W.CDS.list)
        cds = Seq(str(splice(MUU.orig, W.CDS.list)), IUPAC.unambiguous_dna)
        print "\n<b>Old protein:</b>"
        bprint(cds.translate(cds = True, to_stop = True) + '*')
        
        #print MUU.newSplice(W.CDS.list)
        #print splice(MUU.mutated, MUU.newSplice(W.CDS.list))
        cdsm = Seq(str(splice(MUU.mutated, MUU.newSplice(W.CDS.list))), IUPAC.unambiguous_dna)
        print "\n\n<b>New protein:</b>"
        trans = cdsm.translate(to_stop = True)
        if trans[0] != 'M' :
            if str(cdsm[0:3]) in \
               Bio.Data.CodonTable.unambiguous_dna_by_id[1].start_codons :
                bprint('?')
                print "\n\n<b>Alternative protein using start codon %s:</b>" % \
                    str(cdsm[0:3])
                bprint('M' + trans[1:] + '*')
            else :
                bprint('?')
        else :
            bprint(trans + '*')
    #if                
#ppp

def main(cmd) :
    C = Config.Config()

    O = Output.Output(C, __file__)

    O.LogMsg(__file__, "Received variant " + cmd)

    parser = Parser.Nomenclatureparser(O)
    print cmd
    print
    ParseObj = parser.parse(cmd)

    if ParseObj :
        if ParseObj.Version :
            RetrieveRecord = ParseObj.RefSeqAcc + '.' + ParseObj.Version
        else :
            RetrieveRecord = ParseObj.RefSeqAcc
        
        #print "Retrieving..."

        retriever = Retriever.Retriever(C, O)
        record = retriever.loadrecord(RetrieveRecord)
        
        #print "Dicting..."
        D = GenRecord.GenRecord()
        d = D.record2dict(record)
        #print "Printing..."
        #D.printRecordDict(d, record)
        
        import Mutator
        
        MUU = Mutator.Mutator(record.seq)
        ppp(MUU, record, ParseObj, d, "",  0, O)
    #if
    print "\n\n"
    O.Summary()
    O.LogMsg(__file__, "Finished processing variant " + cmd)
#main

if __name__ == "__main__" :
    import sys

    main(sys.argv[1])
