#!/usr/bin/python

"""
Search for an NM number in the MySQL database, if the version number
matches, get the start and end positions in a variant. Translate these
positions to I{g.} notation if the variant is in I{c.} notation or vice versa.

    - If no end position is present, the start position is assumed to be the end
    position.
    - If the version number is not found in the database, an error message is
    generated and a suggestion for an other version is given.
    - If the reference sequence is not found at all, an error is returned.
    - If no variant is present, the transcription start and end and CDS end in
    I{c.} notation is returned.
    - If the variant is not accepted by the nomenclature parser, a parse error
    will be printed.

@requires: sys
@requires: Modules.Config 
@requires: Modules.Db 
@requires: Modules.Crossmap
@requires: Modules.Serializers.SoapMessage
@requires: Modules.Serializers.Mapping
@requires: Modules.Serializers.Transcript
@requires: Bio.Seq.reverse_complement
@requires: collections.defaultdict
"""

import sys                     # argv
from mutalyzer.grammar import Grammar
from mutalyzer import Config     # Config()
from mutalyzer import Db         # Db(), get_NM_version(), get_NM_info()
from mutalyzer import Crossmap   # Crossmap(), g2x(), x2g(), main2int(),
                               # offset2int(), info()
from mutalyzer.Serializers import SoapMessage, Mapping, Transcript

from Bio.Seq import reverse_complement
from collections import defaultdict

def _sl2il(l) :
    """
    Convert a list of strings to a list of integers.

    @arg l: A list of strings
    @type l: list

    @returns: A list of integers
    @rtype: list
    """

    return [int(s) for s in l]
#__sl2il

def _getcoords(C, Loc, Type) :
    """
    Return main, offset and g positions given either a position in
    I{c.} or in I{g.} notation.

    @arg C: A crossmapper
    @type C: object
    @arg Loc: A location in either I{g.} or I{c.} notation
    @type Loc: object
    @arg Type: The reference type
    @type Type: string
    @returns: triple:
        0. Main coordinate in I{c.} notation
        1. Offset coordinate in I{c.} notation
        2. Position in I{g.} notation
    @rtype: triple (integer, integer, integer)
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

class Converter(object) :
    """
    Converter object docstring
    """

    def __init__(self, build, C, O) :
        """
        Initialise the class.
        
        @arg build: the genome build version of the organism (e.g. hg19 for
        human genome build version 19)
        @type build: string
        @arg C: crossmapper object
        @type C: object
        @arg O: output object
        @type O: object        
        """

        self.build = None
        self.__output = O
        self.__config = C
        self.__database = Db.Mapping(build, C.Db)

        # Populated arguments
        self.parseTree = None
        self.crossmap = None
        self.dbFields = {}
    #__init__

    def _changeBuild(self, build) :
        """
        @todo document me (figure out what is does)
        Change the build if it different from the one previously set?????.
        
        @arg build: the genome build version of the organism (e.g. hg19 for
        human genome build version 19)
        @type build: string
        """

        if self.build != build :
            self.crossmap = None
            self.dbFields = {}
            self.build = build
            self.__database = Db.Mapping(build, C.Db)
        #if
    #_changeBuild

    def _reset(self) :
        self.crossmap = None
        self.dbFields = {}
    #_reset

    def _parseInput(self, variant) :
        """
        Parse a variant.
        
        @arg variant: variant description
        @type variant: string
        
        @return: parsetree object
        @rtype: object
        """
        grammar = Grammar(self.__output)
        parseTree = grammar.parse(variant)
        if not parseTree :
            self.__output.addMessage(__file__, 4, "EPARSE",
                    "Could not parse the given variant")
            return None
        #if
        if parseTree.SingleAlleleVarSet :
            #Only simple mutations
            self.__output.addMessage(__file__, 4, "EPARSE",
                    "Can not process multiple mutation variant")
            return None
        #if
        if not parseTree.RefSeqAcc: #In case of LRG for example
            self.__output.addMessage(__file__, 4, "EONLYGB",
                "Currently we only support GenBank Records")
            return None
        #if
        self.parseTree = parseTree
        return parseTree
    #_parseInput

    def _populateFields(self, Fields) :
        #TODO: ADD Error Messages, unlikely that CDS info is missing
        """
        Create a Mutalyzer compatible exon list.
        
        @todo: ADD Error Messages, unlikely that CDS info is missing.
        
        @arg Fields: dictionary with exon start and end positions taken from the
        MySQL database
        @type Fields: dictionary
        
        @return: Exon list
        @rtype: list
        """

        Fields["exonStarts"] =\
            _sl2il(Fields["exonStarts"].split(',')[:-1])
        Fields["exonEnds"] =\
            _sl2il(Fields["exonEnds"].split(',')[:-1])
        assert(len(Fields["exonStarts"]) == len(Fields["exonEnds"]))

        Fields["cdsStart"] = int(Fields["cdsStart"])
        Fields["cdsEnd"]   = int(Fields["cdsEnd"])

        for i in range(len(Fields["exonStarts"])) :
            Fields["exonStarts"][i] += 1

        # Create Mutalyzer compatible exon list
        Fields["exons"] = []
        for exon in zip(Fields["exonStarts"], Fields["exonEnds"]) :
            Fields["exons"].extend(exon)

        self.dbFields = Fields
        return Fields
    #_populateFields

    def _FieldsFromValues(self, values) :
        """
        Combines labels with the given values to a dictionary.
        (zip returns a list of tuples, where the i-th tuple contains the i-th
        element from each of the argument sequences or iterables.
        dict(arg) creates a new data dictionary, with items taken from arg.)
        
        @arg values: list of values take from the MySQL database
        @type values: list
        
        @return: dictionary with values taken from the MySQL database
        @rtype: dictionary
        """

        Fields = dict(zip(
            ("acc", "txStart", "txEnd", "cdsStart", "cdsEnd",
             "exonStarts", "exonEnds", "geneName",
             "chrom", "strand", "protAcc", "version"),
            values))
        return self._populateFields(Fields)
    #_FieldsFromValues

    def _FieldsFromDb(self, acc, version) :
        """
        Get data from database and populate dbFields dict.
        
        @arg acc: NM_ accession number (without version)
        @type acc: string
        @arg version: version number
        @type version: integer
        """

        if not version :
            version = 0
        version = int(version)
        versions = self.__database.get_NM_version(acc)
        if not versions :
            self.__output.addMessage(__file__, 4, "EACCNOTINDB",
                    "The accession number: %s could not be "
                    "found in our database." % acc)
            self.__output.addOutput("LOVDERR",
                    "Reference sequence not found.")
            return None     #Explicit return of None in case of error
        #if
        else :
            if version in versions :
                Values = self.__database.getAllFields(acc, version)
                return self._FieldsFromValues(Values)
            #if
            if not version :
                self.__output.addMessage(__file__, 4, "ENOVERSION",
                    "Version number missing for %s" % acc)
            else :
                self.__output.addMessage(__file__, 4, "EACCNOTINDB",
                    "The accession number: %s version %s "
                    "could not be found in our database." %
                    (acc, version))
            self.__output.addMessage(__file__, 2, "WDIFFFOUND",
                "We found these versions: %s" %
                (", ".join("%s.%s" % (acc, i) for i in sorted(versions))))

            #LOVD list of versions available
            self.__output.addOutput("LOVDERR",
                    "Reference sequence version not found. "
                    "Available: %s" %
                (", ".join("%s.%s" % (acc, i) for i in sorted(versions))))

            #LOVD Only newest
            #self.__output.addOutput("LOVDERR",
            #        "Reference sequence version not found. "
            #        "Available: %s.%s" % (acc, sorted(versions)[-1]))
            return None
        #else
    #_FieldsFromDb

    def makeCrossmap(self) :
        """ 
        Build the crossmapper.
        
        @todo: ADD Error Messages

        @return: Cross ; A Crossmap object
        @rtype: object
        """

        #TODO: ADD Error Messages
        if not self.dbFields: return None

        CDS = []
        if self.dbFields["cdsStart"] != self.dbFields["cdsEnd"] :
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

    def _coreMapping(self) :
        """
        Build the Mapping ClassSerializer.

        @return: Mapping ; A ClassSerializer object
        @rtype: object
        """

        Cross = self.makeCrossmap()
        if not Cross : 
            return None

        mutation = self.parseTree.RawVar

        if not mutation.StartLoc :
            return None

        # Get the coordinates of the start position
        startmain, startoffset, start_g = \
                _getcoords(Cross, mutation.StartLoc.PtLoc,
                            self.parseTree.RefType)

        # If there is an end position, calculate the coordinates.
        if mutation.EndLoc :
            endmain, endoffset, end_g = \
                _getcoords(Cross, mutation.EndLoc.PtLoc,
                            self.parseTree.RefType)
        else :
            end_g, endmain, endoffset = start_g, startmain, startoffset

        # Assign these values to the Mapping ClassSerializer
        V = Mapping()
        V.startmain     = startmain
        V.startoffset   = startoffset
        V.endmain       = endmain
        V.endoffset     = endoffset
        V.start_g       = start_g
        V.end_g         = end_g
        V.mutationType  = mutation.MutationType

        return V
    #_coreMapping

    def giveInfo(self, accNo) :
        """
        Returns transcription start, transcription end and CDS stop, if
        available.
        
        @arg accNo: transcript (NM_) accession number (with or without version)
        @type accNo: string
        
        @return: transcription start, transcription end and CDS stop
        @rtype: triple
        """

        if '.' not in accNo :
            acc, ver = accNo, None
        else :
            acc, ver = accNo.split('.')
        self._FieldsFromDb(acc, ver)
        CM = self.makeCrossmap()
        if CM : 
            return CM.info()
    #giveInfo

    def mainTranscript(self, accNo) :
        """
        One of the entry points (called by the HTML publisher).

        @arg accNo: The full NM accession number (including version)
        @type accNo: string

        @return: T ; ClassSerializer object with the types trans_start,
        trans_stop and CDS_stop
        @rtype: object

        """

        # Initiate ClassSerializer object
        info = self.giveInfo(accNo)
        T = Transcript()
        if info :
            T.trans_start = info[0]
            T.trans_stop  = info[1]
            T.CDS_stop    = info[2]
        return T
    #mainTranscript

    def mainMapping(self, accNo, mutation) :
        """
        One of the entry points (called by the HTML publisher).
        
        @arg accNo: transcript (NM_) accession number (with version?)
        @type accNo: string
        @arg mutation: the 'mutation' (e.g. c.123C>T)
        @type mutation: string
        
        @return: ClassSerializer object
        @rtype: object
        """

        variant = "%s:%s" % (accNo, mutation)
        if self._parseInput(variant) :
            acc = self.parseTree.RefSeqAcc
            version = self.parseTree.Version
            self._FieldsFromDb(acc, version)

        mapping = self._coreMapping()
        soaperrors = self.__output.getSoapMessages()

        if mapping is None :         # Something went wrong
            mapping = Mapping()
            mapping.errorcode = len(soaperrors)
        else :
            mapping.errorcode = 0

        mapping.messages = soaperrors

        return mapping
    #main_Mapping

    def c2chrom(self, variant) :
        """
        Converts a complete HGVS I{c.} notation into a chromosomal notation.

        @arg variant: The variant in HGVS I{c.} notation
        @type variant: string

        @return: var_in_g ; The variant in HGVS I{g.} notation
        @rtype: string
        """

        if self._parseInput(variant) :
            acc = self.parseTree.RefSeqAcc
            version = self.parseTree.Version
            self._FieldsFromDb(acc, version)
        #if
        M = self._coreMapping()
        if M is None : 
            return None

        # construct the variant description
        chromAcc = self.__database.chromAcc(self.dbFields["chrom"])
        f_change = self._constructChange(False)
        r_change = self._constructChange(True)
        if self.dbFields["strand"] == "+" :
            change = f_change
        else :
            change = r_change

        if M.start_g != M.end_g :
            if self.dbFields["strand"] == '+' :
                var_in_g = "g.%s_%s%s" % (M.start_g, M.end_g, change)
            else :
                var_in_g = "g.%s_%s%s" % (M.end_g, M.start_g, change)
        #if
        else :
            var_in_g = "g.%s%s" % (M.start_g, change)

        return "%s:%s" % (chromAcc, var_in_g)
    #c2chrom

    def correctChrVariant(self, variant) :
        """
        @arg variant:
        @type variant: string
        
        @return: variant ; 
        @rtype: string
        """

        #Pre split check
        if ':' not in variant :
            self.__output.addMessage(__file__, 4, "EPARSE",
                "The variant needs a colon")
            return None

        #Remove whitespace
        variant = variant.replace(" ","")

        if variant.startswith("chr") :
            preco, postco = variant.split(":")
            chrom = self.__database.chromAcc(preco)
            if chrom is None :
                self.__output.addMessage(__file__, 4, "ENOTINDB",
                    "Chromosome %s could not be found in our database" %
                    preco)
                return None
            #if
            else :
                variant = "%s:%s" % (chrom, postco)
        #if
        return variant
    #correctChrVariant

    def chrom2c(self, variant, rt) :
        """
        @arg variant: a variant description
        @type variant: string
        @arg rt: the return type
        @type rt: string
        
        @return: HGVS_notatations ; 
        @rtype: dictionary or list
        """

        if not self._parseInput(variant) :
             return None

        acc = self.parseTree.RefSeqAcc
        version = self.parseTree.Version
        chrom = self.__database.chromName("%s.%s" % (acc, version))
        if not chrom :
            self.__output.addMessage(__file__, 4, "ENOTINDB",
                    "Accession number: %s could not be found in our database" %
                acc)
            return None
        #if
        f_change = self._constructChange(False)
        r_change = self._constructChange(True)

        #FIXME This should be a proper conversion.
        loc = int(self.parseTree.RawVar.StartLoc.PtLoc.Main)
        if self.parseTree.RawVar.EndLoc :
            loc2 = int(self.parseTree.RawVar.EndLoc.PtLoc.Main)
        else :
            loc2 = loc
        transcripts = self.__database.get_Transcripts(\
                chrom, loc-5000, loc2+5000, 1)

        HGVS_notatations = defaultdict(list)
        NM_list = []
        for transcript in transcripts :
            self._reset()
            self._FieldsFromValues(transcript)
            M = self._coreMapping()
            if M is None :
                #balen
                continue
            # construct the variant description
            accNo = "%s.%s" % (self.dbFields["acc"],self.dbFields["version"])
            geneName = self.dbFields["geneName"] or ""
            strand = self.dbFields["strand"] == '+'
            startp = self.crossmap.tuple2string((M.startmain, M.startoffset))
            endp = self.crossmap.tuple2string((M.endmain, M.endoffset))

            if strand :
                change = f_change
            else :
                change = r_change
                startp, endp = endp, startp

            #Check if n or c type
            info = self.crossmap.info()
            if info[0] == '1' and info[1] == info[2] :
                mtype = 'n'
            else :
                mtype = 'c'

            if M.start_g != M.end_g :
                loca = "%s_%s" % (startp, endp)
            else :
                loca = "%s" % startp

            variant = "%s:%c.%s%s" % (accNo, mtype, loca, change)
            HGVS_notatations[geneName].append(variant)
            NM_list.append(variant)
        #for
        if rt == "list" :
            return NM_list
        return HGVS_notatations
    #chrom2c

    def _constructChange(self, revc = False) :
        """
        @todo document me

        @arg revc:
        @type revc:
        
        @return:
        @rtype: string
        """

        p = self.parseTree
        if not p or p.SingleAlleleVarSet :
            return None
        var = p.RawVar

        if revc :
            # todo: if var.Arg1 is unicode, this crashes
            arg1 = reverse_complement(var.Arg1 or "") #imported from Bio.Seq
            arg2 = reverse_complement(var.Arg2 or "")
        #if
        else :
            arg1 = var.Arg1
            arg2 = var.Arg2
        #else

        if var.MutationType == "subst" :
            change = "%s>%s" % (arg1, arg2)
        else :
            change = "%s%s" % (var.MutationType, arg1 or arg2 or "")
        return change
    #_constructChange
#Converter
