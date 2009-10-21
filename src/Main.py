#!/usr/bin/python

import sys

import Retriever
import GenRecord
import Crossmap
import Parser

import types
import Config


def roll(string, start, stop) :
    pattern = string[start:stop]
    i = stop - 1

    while string[i] == pattern[-1] :
        pattern = pattern[1:] + pattern[0]
        i += 1
    #while

    return i - len(pattern)
#roll

def printp(string, depth) :
    print (depth * "  ") + str(string)

def ppp(parts, recordDict, depth) :
    printp("+++ recurse +++", depth)
    printp(repr(parts), depth)
    printp(parts, depth)
    """
    printp(parts[4], depth)
    printp(repr(parts[4]), depth)
    printp(parts[4][0][0][0].Main, 2)
    printp(parts, depth)
    """
    printp("RefSeqAcc: " + str(parts.RefSeqAcc), depth)
    printp("RefType: " + str(parts.RefType), depth)


    printp("Version: " + str(parts.Version), depth)
    if parts.Gene :
        printp("Gene Symbol: " + str(parts.Gene.GeneSymbol), depth)
        printp("Transcript variant: " + str(parts.Gene.TransVar), depth)
        printp("Protein isoform: " + str(parts.Gene.ProtIso), depth)
    if parts.ChimeronSet :
        #printp(str(parts.ChimeronSet), depth)
        for i in parts.ChimeronSet :
            printp("ChimeronSet", depth)
            ppp(i, depth + 1)
        #for
    if parts.MosaicSet :
        #printp(str(parts.MosaicSet), depth)
        for i in parts.MosaicSet :
            printp("MosaicSet", depth)
            ppp(i, depth + 1)
        #for
    if parts.SimpleAlleleVarSet :
        #printp(str(parts.SimpleAlleleVarSet), depth)
        for i in parts.SimpleAlleleVarSet :
            printp("SimpleAlleleVarSet", depth)
            ppp(i, depth + 1)
        #for
            #ppp(i, depth + 1)
    if parts.SingleAlleleVarSet :
        #printp(str(parts.SingleAlleleVarSet), depth)
        print parts.SingleAlleleVarSet.RawVar
        for i in parts.SingleAlleleVarSet :
            printp("SingleAlleleVarSet", depth)
            ppp(i, depth + 1)
        #for
    if parts.MultiAlleleVars :
        #printp(str(parts.MultiAlleleVars), depth)
        for i in parts.MultiAlleleVars :
            printp("MultiAlleleVars", depth)
            ppp(i, depth + 1)
    if parts.RawVar :
        printp("RawVar", depth)
        printp(parts.RawVar, depth + 1)
        printp(parts.RawVar.StartLoc, depth + 1)
        printp(parts.RawVar.StartLoc.PtLoc.Main, depth + 1)
        printp(parts.RawVar.StartLoc.PtLoc.Offset, depth + 1)

        main = int(parts.RawVar.StartLoc.PtLoc.Main)
        if parts.RawVar.StartLoc.PtLoc.MainSgn == '-' :
            main = -main
        if parts.RawVar.StartLoc.PtLoc.Offset :
            offset = int(parts.RawVar.StartLoc.PtLoc.Offset)
        else :
            offset = 0

        printp("mRNA: " + str(recordDict[parts.Gene.GeneSymbol].list["RP11-149I2.2-001"].mRNA.list), depth);
        printp("CDS: " + str(recordDict[parts.Gene.GeneSymbol].list["RP11-149I2.2-001"].CDS.list), depth);
        M = Crossmap.Crossmap(recordDict[parts.Gene.GeneSymbol].list["RP11-149I2.2-001"].mRNA.list,
                              recordDict[parts.Gene.GeneSymbol].list["RP11-149I2.2-001"].CDS.list,
                              recordDict[parts.Gene.GeneSymbol].orientation)
        if parts.RefType == 'c' :
            g_notation = M.x2g(main, offset)
            print "g: " + str(g_notation)
            print "c: " + str(M.g2x(g_notation))
        else :
            g_notation = int(main)
            print "g: " + str(g_notation)
            print "c: " + str(M.g2x(g_notation))

        printp(parts.RawVar.MutationType, depth + 1)
    printp("--- recurse ---", depth)
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
D.printRecordDict(d, record)

ppp(ParseObj, d, 0)
