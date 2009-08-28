#!/usr/bin/python

import sys

import Retriever
import Crossmap
import Parser

import types


def dict2class(d):
    """
        Return a class that has same attributes/values and
        dictionaries key/value
    """

    #see if it is indeed a dictionary
    if type(d) != types.DictType:
       return d

    #define a dummy class
    class Dummy:
        pass

    c = Dummy
    for elem in d.keys():
        c.__dict__[elem] = dict2class(d[elem])

    return c
#dict2class


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
    print (depth * "  ") + string

def ppp(parts, depth) :
    printp("+++ recurse +++", depth)
    #printp(str(repr(parts)), depth)
    printp(str(parts), depth)
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
            printp(str(i), depth)
        #for
            #ppp(i, depth + 1)
    if parts.SingleAlleleVarSet :
        #printp(str(parts.SingleAlleleVarSet), depth)
        for i in parts.SingleAlleleVarSet :
            printp("SingleAlleleVarSet", depth)
            ppp(i, depth + 1)
        #for
    if parts.MultiAlleleVars :
        #printp(str(parts.MultiAlleleVars), depth)
        for i in parts.MultiAlleleVars :
            printp("MultiAlleleVars", depth)
            ppp(i, depth + 1)
        #for
            #for j in i :
            #    printp("Test", depth)
            #    ppp(j, depth + 1)
    if parts.RawVar :
        printp("RawVar", depth)
        printp(str(parts.RawVar), depth + 1)
        #for i in parts.RawVar :
        #    printp(str(i), depth + 1)
        #    ppp(i, depth + 1)
    #if
    """
    for i in parts.VarSet :
        printp(i, depth)
        if i.Nest :
            ppp(i.Nest, depth + 1)
        #if
        if i.StartLoc :
            if i.StartLoc.OptRef :
                printp(str(i.StartLoc.OptRef), depth)
            printp("StartLoc: " + str(i.StartLoc.PtLoc.MainSgn) + \
            str(i.StartLoc.PtLoc.Main) + \
            str(i.StartLoc.PtLoc.OffSgn) + str(i.StartLoc.PtLoc.OffOpt) + \
            str(i.StartLoc.PtLoc.Offset), depth)
        #if
        if i.EndLoc :
            printp(i.EndLoc, depth)
            if i.EndLoc.OptRef :
                printp(i.EndLoc.OptRef, depth)
            printp("EndLoc: " + str(i.EndLoc.PtLoc.MainSgn) + \
                  str(i.EndLoc.PtLoc.Main) + \
                  str(i.EndLoc.PtLoc.OffSgn) + str(i.EndLoc.PtLoc.OffOpt) + \
                  str(i.EndLoc.PtLoc.Offset), depth)
        #if
        printp("Mutation type: " + str(i.MutationType), depth)
    #for
    """
    printp("--- recurse ---", depth)
#ppp

def record2dict(record) :
    recordDict = {}
    for i in  record.features :
        if i.qualifiers and i.qualifiers.has_key("gene") :
            if i.qualifiers.has_key("locus_tag") :
                locus_tag = i.qualifiers["locus_tag"][0]
            else :
                locus_tag = 0
            gene = i.qualifiers["gene"][0]
            if not recordDict.has_key(gene) :
                recordDict[gene] = {}
            if not recordDict[gene].has_key(locus_tag) :
                recordDict[gene][locus_tag] = {}
            if i.type == "mRNA" :
                recordDict[gene][locus_tag]["mRNA"] = {}
                recordDict[gene][locus_tag]["mRNA"]["location"] = i.location
                recordDict[gene][locus_tag]["mRNA"]["list"] = []
                for j in i.sub_features :
                    recordDict[gene][locus_tag]["mRNA"]["list"].append(\
                        j.location)
            #if
            if i.type == "CDS" :
                recordDict[gene][locus_tag]["CDS"] = {}
                recordDict[gene][locus_tag]["CDS"]["location"] = i.location
                recordDict[gene][locus_tag]["CDS"]["list"] = []
                for j in i.sub_features :
                    recordDict[gene][locus_tag]["CDS"]["list"].append(\
                        j.location)
            #if
            if i.type == "gene" :
                recordDict[gene][locus_tag]["orientation"] = i.strand
        #if
    #for
    return recordDict
#record2dict

def printRecordDict(recordDict) :
    for i in recordDict :
        print i
        for j in recordDict[i] :
            print "  " + str(j)
            if recordDict[i][j].has_key("mRNA") :
                print "    mRNA: " + str(recordDict[i][j]["mRNA"]["location"])
                print "    mRNA: " + str(recordDict[i][j]["mRNA"]["list"])
            if recordDict[i][j].has_key("CDS") :
                print "      CDS: " + str(recordDict[i][j]["CDS"]["location"])
                print "      CDS: " + str(recordDict[i][j]["CDS"]["list"])
            #if
        #for
    #for
#printRecordDict

parser = Parser.Nomenclatureparser()
cmd = sys.argv[1]
print cmd
ParseObj = parser._parse_command(cmd)
ppp(ParseObj, 0)

RetrieveRecord = ParseObj.RefSeqAcc + '.' + ParseObj.Version

print "Retrieving..."
retriever = Retriever.Retriever()
record = retriever.loadrecord(RetrieveRecord)
print "Dicting..."
d = record2dict(record)
print "Printing..."
#printRecordDict(d)

def splice(record, splice_sites) :
    transcript = ''

    for i in range(0, len(splice_sites), 2) :
        transcript += record.seq[splice_sites[i]:splice_sites[i + 1]] 

    return transcript
#splice

if ParseObj.Gene :
    for i in d[ParseObj.Gene.GeneSymbol] :
        print i
        print d[ParseObj.Gene.GeneSymbol][i]["orientation"]
        print "trans"
        mRNA = []
        if d[ParseObj.Gene.GeneSymbol][i].has_key("mRNA") :
            # Can we find a nice mRNA list (/join() tag)?
            # If so, make it into a normal list.
            if d[ParseObj.Gene.GeneSymbol][i]["mRNA"]["list"] :
                #print d[ParseObj.Gene.GeneSymbol][i]["mRNA"]["list"]
                for j in d[ParseObj.Gene.GeneSymbol][i]["mRNA"]["list"] :
                    mRNA.append(j.start.position)
                    mRNA.append(j.end.position)
                #for
            # Otherwise use the begin and end position (/location tag).
            else :
                mRNA.append(d[ParseObj.Gene.GeneSymbol]\
                            [i]["mRNA"]["location"].start.position)
                mRNA.append(d[ParseObj.Gene.GeneSymbol]\
                            [i]["mRNA"]["location"].end.position)
            #else
        #if
        print "mRNA: " + str(mRNA)
        print splice(record, mRNA)
        CDS = []
        if d[ParseObj.Gene.GeneSymbol][i].has_key("CDS") :
            # Can we find a nice CDS list (/join() tag)?
            # If so, make it into a normal list.
            if d[ParseObj.Gene.GeneSymbol][i]["CDS"]["list"] :
                for j in d[ParseObj.Gene.GeneSymbol][i]["CDS"]["list"] :
                    CDS.append(j.start.position)
                    CDS.append(j.end.position)
                #for
            # Otherwise use the begin and end position (/location tag).
            else :
                CDS.append(d[ParseObj.Gene.GeneSymbol]\
                            [i]["CDS"]["location"].start.position)
                CDS.append(d[ParseObj.Gene.GeneSymbol]\
                            [i]["CDS"]["location"].end.position)
            #else
        #if
        print "CDS: " + str(CDS)
        print splice(record, CDS)
    #for


#d = {}
#d["bla"] = {}
#d["bla"]["xxx"] = "test"
#c = dict2class(d)

#print "---"
#
#RNAf = [5002, 5125, 27745, 27939, 58661, 58762, 74680, 74767, 103409, 103528, 119465, 119537, 144687, 144810, 148418, 149215]
#CDSf = [27925, 27939, 58661, 58762, 74680, 74736]
#
#map = Crossmap.Crossmap(RNAf, CDSf, 1)
#print map.test("*31+d100")
#print map.g2x(map.test("*31+d100"))
#print map.test("*31-u100")
#print map.g2x(map.test("*31-u100"))
#print map.test("30+100")
#print map.g2x(map.test("31+100"))
#print map.test("-31-100")
#print map.g2x(map.test("-31-100"))

roll("AATAATAATAATAATCCCCCC", 3, 6)
