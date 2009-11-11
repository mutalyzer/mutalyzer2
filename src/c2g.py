#!/usr/bin/python

import sys
import Config
import Db
import Crossmap
import Parser

Conf = Config.Config()
Database = Db.Db(Conf)

LOVD_ver = sys.argv[1]
build = sys.argv[2]
accno = sys.argv[3].split('.')[0]
version = sys.argv[3].split('.')[1]
db_version = Database.get_NM_version(accno)
if int(db_version) != int(version) :
    print 0
    print 0
    print 0
    print 0
    exit()
#if
result = Database.get_NM_info(accno)

exonStarts = result[0].split(',')
exonEnds = result[1].split(',')
cdsStart = result[2]
cdsEnd = result[3]
strand = result[4]

#if cdsStart > cdsEnd :
#    temp = cdsStart
#    cdsStart = cdsEnd
#    cdsEnd = temp
##if

mRNA = []
CDS = [int(cdsStart) + 1]
for i in range(len(exonStarts) - 1) :
    mRNA.append(int(exonStarts[i]) + 1)
    mRNA.append(int(exonEnds[i]))
    if int(exonStarts[i]) >= int(cdsStart) and \
        int(exonStarts[i]) <= int(cdsEnd) :
        CDS.append(int(exonStarts[i]) + 1)
    if int(exonEnds[i]) >= int(cdsStart) and \
        int(exonEnds[i]) <= int(cdsEnd) :
        CDS.append(int(exonEnds[i]))
#for
CDS.append(int(cdsEnd))

orientation = 1
if strand == '-' :
    orientation = -1

Cross = Crossmap.Crossmap(mRNA, CDS, orientation)

P = Parser.Nomenclatureparser()
parsetree = P.parse("NM_0000:" + sys.argv[4])
print Cross.c2str(parsetree.RawVar.StartLoc.PtLoc.MainSgn,
                  parsetree.RawVar.StartLoc.PtLoc.Main,
                  parsetree.RawVar.StartLoc.PtLoc.OffSgn,
                  parsetree.RawVar.StartLoc.PtLoc.Offset)
if parsetree.RawVar.EndLoc :
    print Cross.c2str(parsetree.RawVar.EndLoc.PtLoc.MainSgn,
                      parsetree.RawVar.EndLoc.PtLoc.Main,
                      parsetree.RawVar.EndLoc.PtLoc.OffSgn,
                      parsetree.RawVar.EndLoc.PtLoc.Offset)
print Cross.c2g(parsetree.RawVar.StartLoc.PtLoc.MainSgn,
                parsetree.RawVar.StartLoc.PtLoc.Main,
                parsetree.RawVar.StartLoc.PtLoc.OffSgn,
                parsetree.RawVar.StartLoc.PtLoc.Offset)
if parsetree.RawVar.EndLoc :
    print Cross.c2g(parsetree.RawVar.EndLoc.PtLoc.MainSgn,
                    parsetree.RawVar.EndLoc.PtLoc.Main,
                    parsetree.RawVar.EndLoc.PtLoc.OffSgn,
                    parsetree.RawVar.EndLoc.PtLoc.Offset)
else :
    print "0"
    print "0"
#else
