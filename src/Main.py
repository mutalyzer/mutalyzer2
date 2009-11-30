#!/usr/bin/python

import sys

import Retriever
import GenRecord
import Crossmap
import Parser

import types
import Config


def roll(string, start, stop, orientation) :
    pattern = string[start:stop]
    if orientation == 1 :
        i = stop - 1

        while i < len(string) and string[i] == pattern[-1] :
            pattern = pattern[1:] + pattern[0]
            i += 1
        #while

        return i - len(pattern) 
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

def ErrorMsg(message) :
    print "An error occurred:"
    print "    ", message
#ErrorMsg

def WarningMsg(message) :
    print "Warning:"
    print "    ", message
#WarningMsg

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

def ppp(parts, recordDict, refseq, depth) :
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
    printp("RefSeqAcc: " + str(refseq), depth)
    printp("RefType: " + str(parts.RefType), depth)


    printp("Version: " + str(parts.Version), depth)
    if parts.Gene :
        printp("Gene Symbol: " + str(parts.Gene.GeneSymbol), depth)
        #printp("Transcript variant: " + str(parts.Gene.TransVar), depth)
        #printp("Protein isoform: " + str(parts.Gene.ProtIso), depth)
    #if
    if parts.ChimeronSet :
        #printp(str(parts.ChimeronSet), depth)
        for i in parts.ChimeronSet :
            printp("ChimeronSet", depth)
            ppp(i, recordDict, refseq, depth + 1)
        #for
    #if
    if parts.MosaicSet :
        #printp(str(parts.MosaicSet), depth)
        for i in parts.MosaicSet :
            printp("MosaicSet", depth)
            ppp(i, recordDict, refseq, depth + 1)
        #for
    #if
    if parts.SimpleAlleleVarSet :
        #printp(str(parts.SimpleAlleleVarSet), depth)
        for i in parts.SimpleAlleleVarSet :
            printp("SimpleAlleleVarSet", depth)
            ppp(i, recordDict, refseq, depth + 1)
        #for
            #ppp(i, recordDict, refseq, depth + 1)
    #if
    if parts.SingleAlleleVarSet :
        #printp(str(parts.SingleAlleleVarSet), depth)
        print parts.SingleAlleleVarSet.RawVar
        for i in parts.SingleAlleleVarSet :
            printp("SingleAlleleVarSet", depth)
            ppp(i, recordDict, refseq, depth + 1)
        #for
    #if
    if parts.MultiAlleleVars :
        #printp(str(parts.MultiAlleleVars), depth)
        for i in parts.MultiAlleleVars :
            printp("MultiAlleleVars", depth)
            ppp(i, recordDict, refseq, depth + 1)
        #for
    #if
    if parts.RawVar :
        start_main = PtLoc2main(parts.RawVar.StartLoc.PtLoc)
        start_offset = PtLoc2offset(parts.RawVar.StartLoc.PtLoc)
        end_main = start_main
        end_offset = start_offset
        if parts.RawVar.EndLoc :
            end_main = PtLoc2main(parts.RawVar.EndLoc.PtLoc)
            end_offset = PtLoc2offset(parts.RawVar.EndLoc.PtLoc)
        #if

        transcriptvariant = "001" 
        if parts.Gene and parts.Gene.TransVar :
            transcriptvariant = parts.Gene.TransVar

        M = Crossmap.Crossmap(
          recordDict[parts.Gene.GeneSymbol].list[transcriptvariant].mRNA.list,
          recordDict[parts.Gene.GeneSymbol].list[transcriptvariant].CDS.list,
          recordDict[parts.Gene.GeneSymbol].orientation)

        start_g = int(start_main)
        end_g = int(end_main)
        if parts.RefType == 'c' :
            start_g = M.x2g(start_main, start_offset)
            end_g = M.x2g(end_main, end_offset)
        #if

        if recordDict[parts.Gene.GeneSymbol].orientation == -1 :
            temp = start_g
            start_g = end_g
            end_g = temp
        #if

        print str(record.seq[start_g - 20:start_g + 20])

        # Substitution.
        if parts.RawVar.MutationType == "subst" :
            if record.seq[start_g - 1] != parts.RawVar.Arg1 :
                ErrorMsg("No nucleotide %c at position c.%s (g.%i), found a " \
                    "%c instead." % (parts.RawVar.Arg1, M.g2c(start_g), 
                    start_g, record.seq[start_g - 1]))
            if parts.RawVar.Arg1 == parts.RawVar.Arg2 :
                ErrorMsg("No mutation given (%c>%c) at position c.%s (g.%i)." % 
                    (parts.RawVar.Arg1, parts.RawVar.Arg1, M.g2c(start_g), 
                    start_g))
            MUU.subM(start_g, parts.RawVar.Arg2)
        #if

        # Deletion / Duplication.
        if parts.RawVar.MutationType == "del" or \
           parts.RawVar.MutationType == "dup" :
            if parts.RawVar.Arg1 and \
                parts.RawVar.Arg1 != str(record.seq[start_g - 1:end_g]) :
                if not IsInt(parts.RawVar.Arg1) :
                    ErrorMsg("String %s not found at position c.%s (g.%i), " \
                        "found %s instead." % (parts.RawVar.Arg1, 
                        M.g2c(start_g), start_g, 
                        str(record.seq[start_g - 1:end_g])))
                else :
                    if end_g - start_g + 1 != int(parts.RawVar.Arg1) :
                        ErrorMsg("The length of the deletion (%i) at " \
                            "position c.%s (g.%i) differed from that of the " \
                            "range (%i)." % (int(parts.RawVar.Arg1), 
                            M.g2c(start_g), start_g, end_g - start_g + 1))
            #if
            rollposstart = roll(record.seq, start_g - 1, end_g,
                           recordDict[parts.Gene.GeneSymbol].orientation)
            if rollposstart != start_g :
                rollposend = rollposstart + (end_g - start_g)
                WarningMsg("String %s at position c.%s (g.%i) was given, " \
                    "however, the HGVS notation prescribes that it should be " \
                    "%s at position c.%s (g.%i)." % (
                    str(record.seq[start_g - 1:end_g]), M.g2c(start_g), start_g,
                    str(record.seq[rollposstart - 1:rollposend]), 
                    M.g2c(rollposstart), rollposstart))
            #if
            if parts.RawVar.MutationType == "del" :
                MUU.delM(start_g, end_g)
            else :
                MUU.dupM(start_g, end_g)
                print start_g, " ", end_g
                print "oioioi"
        #if

        # Inversion.
        if parts.RawVar.MutationType == "inv" :
            snoop = palinsnoop(record.seq[start_g - 1:end_g])
            if snoop :
                if snoop == -1 :
                    WarningMsg("String %s at position c.%s (g.%i) is a " \
                        "palindrome (its own reverse complement)." % (
                        str(record.seq[start_g - 1:end_g]), 
                        M.g2c(start_g), start_g))
                else :
                    WarningMsg("String %s at position c.%s (g.%i) is a " \
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
        if parts.RawVar.MutationType == "ins" :
            if start_g + 1 != end_g :
                ErrorMsg("c.%s (g.%i) and c.%s (g.%i) are not consecutive " \
                    "positions." % (M.g2c(start_g), start_g, M.g2c(end_g),
                    end_g))

            MUU.insM(start_g, parts.RawVar.Arg1)

            # Niet record.seq, maar de gemuteerde seq!
            rollposstart = roll(record.seq, start_g - 1, 
                           start_g + len(parts.RawVar.Arg1),
                           recordDict[parts.Gene.GeneSymbol].orientation)
            if rollposstart != start_g :
                rollposend = rollposstart + len(parts.RawVar.Arg1)
                WarningMsg("Insertion of %s at position c.%s (g.%i) was " \
                    "given, however, the HGVS notation prescribes that it " \
                    "should be a duplication of %s at position c.%s (g.%i)." % (
                    parts.RawVar.Arg1, M.g2c(start_g), start_g, 
                    str(record.seq[rollposstart - 1:rollposend]), 
                    M.g2c(rollposstart), rollposstart))
            #if
        #if

        print MUU.mutated[start_g - 20:start_g + 20]
    #if
#ppp

parser = Parser.Nomenclatureparser()
cmd = sys.argv[1]
print cmd
ParseObj = parser.parse(cmd)

RetrieveRecord = ParseObj.RefSeqAcc + '.' + ParseObj.Version

print "Retrieving..."
C = Config.Config()
retriever = Retriever.Retriever(C)
record = retriever.loadrecord(RetrieveRecord)

print "Dicting..."
D = GenRecord.GenRecord()
d = D.record2dict(record)
print "Printing..."
#D.printRecordDict(d, record)

import Mutator

MUU = Mutator.Mutator(record.seq)
ppp(ParseObj, d, "",  0)
