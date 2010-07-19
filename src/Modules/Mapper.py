#!/usr/bin/python

"""
    Search for an NM number in the MySQL database, if the version number
    matches, get the start and end positions in a variant and translate these
    positions to g. notation if the variant is in c. notation and vice versa.

    - If no end position is present, the start position is assumed to be the
      end position.
    - If the version number is not found in the database, an error message is
      generated and a suggestion for an other version is given.
    - If the reference sequence is not found at all, an error is returned.
    - If no variant is present, the transcription start and end and CDS end
      in c. notation is returned.
    - If the variant is not accepted by the nomenclature parser, a parse error
      will be printed.

"""
import sys                     # argv
from Modules import Config     # Config()
from Modules import Db         # Db(), get_NM_version(), get_NM_info()
from Modules import Crossmap   # Crossmap(), g2x(), x2g(), main2int(),
                               # offset2int(), info()
from Modules import Parser     # Nomenclatureparser(), parse()
from Modules import Output     # Output(), LogMsg()
from ZSI.fault import Fault    # Fault()
from ZSI import TC             # Struct()

from Bio.Seq import reverse_complement
from collections import defaultdict

from soaplib.serializers.primitive import String, Integer, Array
from soaplib.serializers.clazz import ClassSerializer

#def debug(f):
#    """Debug Wrapper for functions called from within the daemon"""
#    def _tempf(*args):
#        with open("/tmp/mapper.out", "a+") as of:
#            try:
#                of.write("\nFunction %s\n\targs: %s\n\t" % (`f`, `args`))
#                ret = f(*args)
#                of.write("Returns: %s" % `ret`)
#                return ret
#            except Exception, e:
#                import traceback
#                of.write("\nEXCEPTION:\n")
#                traceback.print_exc(file=of)
#    return _tempf

class SoapMessage(ClassSerializer):
    """Send info message over the soapline"""

    class types():
        errorcode = String
        message = String

    def __init__(self):
        self.typecode = TC.Struct(SoapMessage, [
            TC.String("errorcode"),
            TC.String("message")], "SoapMessage")
#SoapMessage

class Mapping(ClassSerializer) :
    """
        Extended ClassSerializer object with mixed types of attributes

        Attributes:
            startmain ; Define the type of startmain.
            startoffset ; Define the type of startoffset.
            endmain ; Define the type of endmain value.
            endoffset ; Define the type of endoffset value.
            start_g ; Define the type of start_g value.
            end_g ; Define the type of end_g value.
            mutationType ; Define the type of mutation type
    """

    class types() :
        """
            Types are defined here for the soaplib module.
        """

        startmain = Integer
        startoffset = Integer
        endmain = Integer
        endoffset = Integer
        start_g = Integer
        end_g = Integer
        mutationType = String
        errorcode = Integer
        messages = Array(SoapMessage)
    #types

    def __init__(self) :
        """
            Types are defined here for the TC module.
        """

        self.typecode = TC.Struct(Mapping, [
            TC.Integer('startmain'),
            TC.Integer('startoffset'),
            TC.Integer('endmain'),
            TC.Integer('endoffset'),
            TC.Integer('start_g'),
            TC.Integer('end_g'),
            TC.String('mutationType'),
            TC.Integer("errorcode"),
            TC.Array("SoapMessage", TC.Struct(SoapMessage, [
                TC.String("errorcode"),
                TC.String("message")], "SoapMessage"), "messages")
            ], 'Mapping')
    #__init__
#Mapping

class Transcript(ClassSerializer) :
    """
        Extended ClassSerializer object with mixed types of attributes

        Attributes:
            trans_start ; Define the type of trans_start
            trans_stop  ; Define the type of trans_stop
            CDS_stop    ; Define the type of CDS_stop
    """

    class types() :
        """
        """

        trans_start = Integer
        trans_stop = Integer
        CDS_stop = Integer
    #types

    def __init__(self) :
        """
        """

        self.typecode = TC.Struct(Transcript, [
            TC.Integer('trans_start'),
            TC.Integer('trans_stop'),
            TC.Integer('CDS_stop')
            ], 'Transcript')
    #__init__
#Transcript

def _sl2il(l) :
    """
        Convert a list of strings to a list of integers.

        Arguments: l ; A list of strings.

        Returns: list ; A list of integers.
    """

    for i in range(len(l)) :
        l[i] = int(l[i])
    return l
#__sl2il

def _getcoords(C, Loc, Type) :
    """
        Return main, offset and g positions given either a position in
        c. or in g. notation.

        Arguments:
            C    ; A crossmapper.
            Loc  ; Either a location in g. or c. notation.
            Type ; The reference type.
        Returns:
            triple:
                0 ; Main coordinate in c. notation.
                1 ; Offset coordinate in c. notation.
                2 ; Position in g. notation.
    """

    if Type == 'c' :
        main = C.main2int(Loc.MainSgn +  Loc.Main)
        offset = C.offset2int(Loc.OffSgn +  Loc.Offset)
        g = C.x2g(main, offset)
        main, offset = C.g2x(g)
    #if
    else :
        g = int(Loc.Main)
        main, offset = C.g2x(g)
    #else
    return (main, offset, g)
#__getcoords

class Converter(object):
    """
        Converter object docstring
    """
    def __init__(self, build, C, O):
        self.build = None
        self.__output = O
        self.__config = C
        self.__database = Db.Mapping(build, C.Db)

        # Populated arguments
        self.parseTree = None
        self.crossmap = None
        self.dbFields = {}

    def _changeBuild(self, build):
        if self.build != build:
            self.crossmap = None
            self.dbFields = {}
            self.build = build
            self.__database = Db.Mappin(build, C.Db)
    
    def _reset(self):
        self.crossmap = None
        self.dbFields = {}

    def _parseInput(self, variant):
        P = Parser.Nomenclatureparser(self.__output)
        parseTree = P.parse(variant)
        if not parseTree:
            self.__output.addMessage(__file__, 4, "EPARSE",
                    "Could not parse the given variant")
            return None
        if parseTree.SingleAlleleVarSet:
            #Only simple mutations
            self.__output.addMessage(__file__, 4, "EPARSE",
                    "Can not process multiple mutation variant")
            return None
        if not parseTree.RefSeqAcc: #In case of LRG for example
            self.__output.addMessage(__file__, 4, "EONLYGB",
                "Currently we only support GenBank Records")
            return None
        self.parseTree = parseTree
        return parseTree
    #_parseInput

    def _populateFields(self, Fields):
        #TODO: ADD Error Messages, unlikely that CDS info is missing


        Fields["exonStarts"] =\
            _sl2il(Fields["exonStarts"].split(',')[:-1])
        Fields["exonEnds"] =\
            _sl2il(Fields["exonEnds"].split(',')[:-1])
        Fields["cdsStart"] =\
             int(Fields["cdsStart"])
        Fields["cdsEnd"] =\
            int(Fields["cdsEnd"])

        # Create Mutalyzer compatible exon list
        Fields["exons"] = []
        exons = zip(Fields["exonStarts"], Fields["exonEnds"])
        for exon in exons:
            Fields["exons"].extend(exon)

        self.dbFields = Fields
        return Fields

    def _FieldsFromValues(self, values):
        Fields = dict(zip(
            ("acc", "txStart", "txEnd", "cdsStart", "cdsEnd",
             "exonStarts", "exonEnds", "geneName",
             "chrom", "strand", "protAcc", "version"),
            values))
        return self._populateFields(Fields)

    def _FieldsFromDb(self, acc, version):
        """Get data from database and populate dbFields dict"""
        Values = self.__database.getAllFields(acc, version)
        if Values[-1] is None:
            versions = self.__database.get_NM_version(acc)
            if not versions:
                self.__output.addMessage(__file__, 4, "ENOTFOUND",
                        "The accession number: %s could not be "
                        "found in our database." % acc)
            else:
                if not version:
                    self.__output.addMessage(__file__, 4, "ENOVERSION",
                        "Version number missing for %s" % acc)
                else:
                    self.__output.addMessage(__file__, 4, "ENOTFOUND",
                        "The accession number: %s version %s "
                        "could not be found in our database." %
                        (acc, version))
                self.__output.addMessage(__file__, 2, "WDIFFFOUND",
                    "We found these versions: %s" %
                    (", ".join("%s.%s" % (acc, i) for i in sorted(versions))))
            return None     #Excplicit return of None in case of error
        return self._FieldsFromValues(Values)

    def makeCrossmap(self) :
        ''' Build the crossmapper

            Arguments:
                build   ; The human genome build
                acc     ; The NM accession number, including version.

            Returns:
                Cross ; A Crossmap object.
        '''
        #TODO: ADD Error Messages
        if not self.dbFields: return None

        CDS = []
        if self.dbFields["cdsStart"] != self.dbFields["cdsEnd"]:
            #The counting from 0 conversion.
            CDS = [self.dbFields["cdsStart"] + 1, self.dbFields["cdsEnd"]]

        mRNA = self.dbFields["exons"]

        # Convert the strand information to orientation.
        orientation = 1
        if self.dbFields["strand"] == '-' :
            orientation = -1

        # Build the crossmapper.
        self.crossmap = Crossmap.Crossmap(mRNA, CDS, orientation)
        return self.crossmap
    #makeCrossmap

    def _coreMapping(self):
        Cross = self.makeCrossmap()
        if not Cross: return None

        mutation = self.parseTree.RawVar

        startmain, startoffset, start_g = \
                _getcoords(Cross, mutation.StartLoc.PtLoc,
                            self.parseTree.RefType)

        # Assume there is no end position given.
        end_g = start_g
        endmain = startmain
        endoffset = startoffset

        # If there is an end position, calculate the coordinates.
        if mutation.EndLoc :
            endmain, endoffset, end_g = \
                _getcoords(Cross, mutation.EndLoc.PtLoc,
                            self.parseTree.RefType)

        # Assign these values to the Mapping V types.
        V = Mapping()
        V.startmain = startmain
        V.startoffset = startoffset
        V.endmain = endmain
        V.endoffset = endoffset
        V.start_g = start_g
        V.end_g = end_g
        V.mutationType = mutation.MutationType

        return V

    def giveInfo(self, accNo) :
        self._FieldsFromDb(*accNo.split('.'))
        self.makeCrossmap()
        return self.crossmap.info()

    def mainMapping(self, accNo, mutation) :
        """
            One of the entry points (called by the HTML publisher).
        """
        variant = "%s:%s" % (accNo, mutation)
        if self._parseInput(variant):
            acc = self.parseTree.RefSeqAcc
            version = self.parseTree.Version
            self._FieldsFromDb(acc, version)

        mapping = self._coreMapping()
        soaperrors = self.__output.getSoapMessages()

        if mapping is None:         # Something went wrong
            mapping = Mapping()
            mapping.errorcode = len(soaperrors)
        else:
            mapping.errorcode = 0

        mapping.messages = soaperrors

        return mapping

    #main_Mapping

    def c2chrom(self, variant) :
        """
           Converts a complete HGVS c. notation into a chromosomal notation

           Arguments:
               variant  ; The variant in HGVS c.notation.

           Returns:
               var_in_g ; The variant in HGVS g. notation (string).

        """
        if self._parseInput(variant):
            acc = self.parseTree.RefSeqAcc
            version = self.parseTree.Version
            self._FieldsFromDb(acc, version)
        M = self._coreMapping()
        if M is None: return None

        # construct the variant description
        chromAcc = self.__database.chromAcc(self.dbFields["chrom"])
        f_change = self._constructChange(False)
        r_change = self._constructChange(True)
        if self.dbFields["strand"] == "+":
            change = f_change
        else:
            change = r_change

        if M.start_g != M.end_g:
            if self.dbFields["strand"] == '+' :
                var_in_g = "g.%s_%s%s" % (M.start_g, M.end_g, change)
            else :
                var_in_g = "g.%s_%s%s" % (M.end_g, M.start_g, change)
        else :
            var_in_g = "g.%s%s" % (M.start_g, change)

        return "%s:%s" % (chromAcc, var_in_g)
    #cTog

    def correctChrVariant(self, variant):
        if variant.startswith("chr"):
            if ':' not in variant:
                self.__output.addMessage(__file__, 4, "EPARSE",
                    "The variant needs a colon")
                return None
            preco, postco = variant.split(":")
            chrom = self.__database.chromAcc(preco)
            if chrom is None:
                self.__output.addMessage(__file__, 4, "ENOTINDB",
                    "Chromosome %s could not be found in our database" %
                    preco)
                return None
            else:
                variant = "%s:%s" % (chrom, postco)
        return variant

    def chrom2c(self, variant):
        if not self._parseInput(variant):
             return None

        acc = self.parseTree.RefSeqAcc
        version = self.parseTree.Version
        chrom = self.__database.chromName("%s.%s" % (acc, version))
        if not chrom:
            self.__output.addMessage(__file__, 4, "ENOTINDB",
                    "Accession number: %s could not be found in our database" %
                acc)
            return None
        f_change = self._constructChange(False)
        r_change = self._constructChange(True)

        loc = int(self.parseTree.RawVar.StartLoc.PtLoc[0])
        if self.parseTree.RawVar.EndLoc:
            loc2 = int(self.parseTree.RawVar.EndLoc.PtLoc[0])
        else:
            loc2 = loc
        transcripts = self.__database.get_Transcripts(\
                chrom, loc-5000, loc2+5000, 1)

        HGVS_notatations = defaultdict(list)
        for transcript in transcripts:
            self._reset()
            self._FieldsFromValues(transcript)
            M = self._coreMapping()
            if M is None:
                #balen
                continue
            # construct the variant description
            accNo = "%s.%s" % (self.dbFields["acc"],self.dbFields["version"])
            geneName = self.dbFields["geneName"] or ""
            strand = self.dbFields["strand"] == '+'
            startp = self.crossmap.tuple2string((M.startmain, M.startoffset))
            endp = self.crossmap.tuple2string((M.endmain, M.endoffset))

            if strand:
                change = f_change
            else:
                change = r_change
                startp, endp = endp, startp

            #Check if n or c type
            info = self.crossmap.info()
            if info[0] == '1' and info[1] == info[2]:
                mtype = 'n'
            else:
                mtype = 'c'

            if M.start_g != M.end_g:
                loca = "%s_%s" % (startp, endp)
            else:
                loca = "%s" % startp

            variant = "%s:%c.%s%s" % (accNo, mtype, loca, change)
            HGVS_notatations[geneName].append(variant)
        return HGVS_notatations

    #chrom2c

    def _constructChange(self, revc=False):
        p = self.parseTree
        if not p or p.SingleAlleleVarSet:
            return None
        var = p.RawVar

        if revc:
            arg1 = reverse_complement(var.Arg1 or "") #imported from Bio.Seq
            arg2 = reverse_complement(var.Arg2 or "")
        else:
            arg1 = var.Arg1
            arg2 = var.Arg2

        if var.MutationType == "subst":
            change = "%s>%s" % (arg1, arg2)
        else:
            change = "%s%s" % (var.MutationType, arg1 or arg2 or "")
        return change
    #_constructChange


def mainTranscript(build, acc, C, Cross) :
    """
        One of the entry points (called by the HTML publisher).

        Arguments:
            build   ; The human genome build.
            acc     ; The full NM accession number (including version).
            C       ; Conf object

        Returns:
            T       ; ClassSerializer object with the types trans_start,
                      trans_stop and CDS_stop.

    """

    # Initiate ClassSerializer object
    T = Transcript()
    # Return transcription start, transcription end and
    # CDS stop in c. notation.
    info = Cross.info()
    # Assign these values to the Transcript T types.
    T.trans_start = info[0]
    T.trans_stop  = info[1]
    T.CDS_stop    = info[2]

    return T
#mainTranscript
