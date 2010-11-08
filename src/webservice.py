#!/usr/bin/python

"""
    Mutalyzer webservices.

    Public classes:
        MutalyzerService ; Mutalyzer webservices.
"""

from soaplib.wsgi_soap import SimpleWSGISoapApp
from soaplib.service import soapmethod
from soaplib.serializers.primitive import String, Integer, Array
from ZSI.fault import Fault

import Mutalyzer                                
from Modules import Web
from Modules import Db
from Modules import Output
from Modules import Config
from Modules import Parser
from Modules import Mapper
from Modules.Serializers import SoapMessage, Mapping, Transcript, \
                                MutalyzerOutput

class MutalyzerService(SimpleWSGISoapApp) :
    """
        Mutalyzer webservices.

        These methods are made public via a SOAP interface.

        Private methods:
            __checkBuild(L, D, build) ; Check if the build is supported.
            __checkChrom(L, D, chrom) ; Check if the chromosome is in our
                                        database.
            __checkPos(L, pos)        ; Check if the position is valid.

        Public methods:
            getTranscripts(build, chrom, ; Get all transcripts that overlap
                           pos)            with a chromosomal position.
            getTranscriptsRange(build,   ; Get all transcripts that overlap
                                chrom,     with a range on a chromosome.
                                pos1,
                                pos2,
                                method)
            getGeneName(build, accno)    ; Find the gene name associated with a
                                           transcript.
            mappingInfo(LOVD_ver, build, ; Convert a transcript coordinate to a
                        accNo, variant)    chromosomal one, or vice versa.
            transcriptInfo(LOVD_ver,     ; Find transcription start and end,
                           build,          and CDS end (in c. notation) for a
                           accNo)          given transcript.
            cTogConversion(self, build,  ; Convert c. to g.
                           variant)
            gTocConversion(self, build,  ; Convert g. to c.
                           variant)
    """

    def __checkBuild(self, build, config) :
        """
            Check if the build is supported (hg18 or hg19).

            Arguments:
                L     ; An output object for logging.
                D     ; A handle to the database.
                build ; The build name that needs to be checked.

            Returns:
                Nothing (but raises an EARG exception).
        """

        if not build in config.dbNames :
            L.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault(Fault.Client, "EARG",
                detail = "The build argument (%s) was not a valid " \
                         "build name." % build)
        #if
    #__checkBuild

    def __checkChrom(self, L, D, chrom) :
        """
            Check if the chromosome is in our database.

            Arguments:
                L     ; An output object for logging.
                D     ; A handle to the database.
                chrom ; The name of the chromosome.

            Returns:
                Nothing (but raises an EARG exception).
        """

        if not D.isChrom(chrom) :
            L.addMessage(__file__, 4, "EARG", "EARG %s" % chrom)
            raise Fault(Fault.Client, "EARG",
                detail = "The chrom argument (%s) was not a valid " \
                         "chromosome name." % chrom)
        #if
    #__checkChrom

    def __checkPos(self, L, pos) :
        """
            Check if the position is valid.

            Arguments:
                L   ; An output object for logging.
                pos ; The position.

            Returns:
                Nothing (but raises an ERANGE exception).
        """

        if pos < 1 :
            L.addMessage(__file__, 4, "ERANGE", "ERANGE %i" % pos)
            raise Fault(Fault.Client, "ERANGE",
                detail = "The pos argument (%i) is out of range." % pos)
        #if
    #__checkPos

    def __checkVariant(self, L, variant) :
        """
            Check if a variant is provided.

            Arguments:
                L       ; An output object for logging.
                variant ; The variant.

            Returns:
                Nothing (but raises an EARG exception).
        """

        if not variant :
            L.addMessage(__file__, 4, "EARG", "EARG no variant")
            raise Fault(Fault.Client, "EARG",
                detail = "The variant argument is not provided.")
        #if
    #__checkVariant

    @soapmethod(String, String, Integer, _returns = Array(String))
    def getTranscripts(self, build, chrom, pos) :
        """
            Get all the transcripts that overlap with a chromosomal position.

            Arguments:
                string build ; The build name encoded as "hg18" or "hg19".
                string chrom ; A chromosome encoded as "chr1", ..., "chrY".
                int    pos   ; A postion on the chromosome.

            Returns:
                string ; A list of transcripts.

            On error an exception is raised:
                detail       ; Human readable description of the error.
                faultstring: ; A code to indicate the type of error.
                    EARG   ; The argument was not valid.
                    ERANGE ; An invalid range was given.
        """

        C = Config.Config()
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request getTranscripts(%s %s %s)" % (build,
                     chrom, pos))

        self.__checkBuild(build, C.Db)
        D = Db.Mapping(build, C.Db)

        self.__checkChrom(L, D, chrom)
        self.__checkPos(L, pos)

        ret = D.get_Transcripts(chrom, pos, pos, True)

        #filter out the accNo
        ret = [r[0] for r in ret]


        L.addMessage(__file__, -1, "INFO",
                     "Finished processing getTranscripts(%s %s %s)" % (build,
                     chrom, pos))

        del D, L, C
        return [ret]
    #getTranscripts

    @soapmethod(String, String, Integer, Integer, Integer,
                _returns = Array(String))
    def getTranscriptsRange(self, build, chrom, pos1, pos2, method) :
        """
            Get all the transcripts that overlap with a range on a chromosome.

            Arguments:
                string build ; The build name encoded as "hg18" or "hg19".
                string chrom  ; A chromosome encoded as "chr1", ..., "chrY".
                int    pos1   ; The first postion of the range.
                int    pos2   ; The last postion of the range.
                int    method ; The method of determining overlap:
                                0 ; Return only the transcripts that completely
                                    fall in the range [pos1, pos2].
                                1 ; Return all hit transcripts.

            Returns:
                string ; A list of transcripts.
        """

        C = Config.Config()
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
            "Received request getTranscriptsRange(%s %s %s %s %s)" % (build,
            chrom, pos1, pos2, method))

        D = Db.Mapping(build, C.Db)
        self.__checkBuild(build, C.Db)

        ret = D.get_Transcripts(chrom, pos1, pos2, method)

        #filter out the accNo
        ret = [r[0] for r in ret]

        L.addMessage(__file__, -1, "INFO",
            "Finished processing getTranscriptsRange(%s %s %s %s %s)" % (
            build, chrom, pos1, pos2, method))

        del D, L, C
        return [ret]
    #getTranscriptsRange

    @soapmethod(String, String, _returns = String)
    def getGeneName(self, build, accno) :
        """
            Find the gene name associated with a transcript.

            Arguments:
                string build ; The build name encoded as "hg18" or "hg19".
                string accno ; The identifier of a transcript.

            Returns:
                string ; The name of the associated gene.
        """

        C = Config.Config()
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request getGeneName(%s %s)" % (build, accno))

        D = Db.Mapping(build, C.Db)
        self.__checkBuild(build, C.Db)

        ret = D.get_GeneName(accno.split('.')[0])

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing getGeneName(%s %s)" % (build, accno))

        del D, L, C
        return ret
    #getGeneName


    @soapmethod(String, String, String, String, _returns = Mapping)
    def mappingInfo(self, LOVD_ver, build, accNo, variant) :
        """
            Search for an NM number in the MySQL database, if the version
            number matches, get the start and end positions in a variant and
            translate these positions to g. notation if the variant is in c.
            notation and vice versa.

            - If no end position is present, the start position is assumed to
              be the end position.
            - If the version number is not found in the database, an error
              message is generated and a suggestion for an other version is
              given.
            - If the reference sequence is not found at all, an error is
              returned.
            - If no variant is present, an error is returned.
            - If the variant is not accepted by the nomenclature parser, a
              parse error will be printed.


            Arguments (all strings):
                LOVD_ver ; The LOVD version.
                build ; The human genome build (hg19 or hg18).
                accNo ; The NM accession number and version.
                variant ; The variant.

            Returns:
                complex object:
                    start_main   ; The main coordinate of the start position
                                   in c. (non-star) notation.
                    start_offset ; The offset coordinate of the start position
                                   in c. notation (intronic position).
                    end_main     ; The main coordinate of the end position in
                                   c. (non-star) notation.
                    end_offset   ; The offset coordinate of the end position in
                                   c. notation (intronic position).
                    start_g      ; The g. notation of the start position.
                    end_g        ; The g. notation of the end position.
                    type         ; The mutation type.

        """

        C = Config.Config()
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Reveived request mappingInfo(%s %s %s %s)" % (
                        LOVD_ver, build, accNo, variant))

        conv = Mapper.Converter(build, C, L)
        result = conv.mainMapping(accNo, variant)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing mappingInfo(%s %s %s %s)" % (
                        LOVD_ver, build, accNo, variant))

        del L, C
        return result
    #mappingInfo

    @soapmethod(String, String, String, _returns = Transcript)
    def transcriptInfo(self, LOVD_ver, build, accNo) :
        """
            Search for an NM number in the MySQL database, if the version
            number matches, the transcription start and end and CDS end
            in c. notation is returned.


            Arguments (all strings:
                LOVD_ver ; The LOVD version.
                build ; The human genome build (hg19 or hg18).
                accNo ; The NM accession number and version.

            Returns:
                complex object:
                    trans_start  ; Transcription start in c. notation.
                    trans_stop   ; Transcription stop in c. notation.
                    CDS_stop     ; CDS stop in c. notation.
        """

        C = Config.Config()
        O = Output.Output(__file__, C.Output)

        O.addMessage(__file__, -1, "INFO",
                     "Received request transcriptInfo(%s %s %s)" % (LOVD_ver,
                     build, accNo))

        converter = Mapper.Converter(build, C, O)
        T = converter.mainTranscript(accNo)

        O.addMessage(__file__, -1, "INFO",
                     "Finished processing transcriptInfo(%s %s %s)" % (
                     LOVD_ver, build, accNo))
        return T
    #transcriptInfo

    @soapmethod(String, String, _returns = String)
    def chromAccession(self, build, name) :
        """
            Get the accession number of a chromosome, given a name.

            Arguments:
                build   ; The human genome build.
                name    ; The name of a chromosome.

            Returns:
                string ; The accession number of a chromosome.
        """
        C = Config.Config() # Read the configuration file.
        D = Db.Mapping(build, C.Db)
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request chromAccession(%s %s)" % (build, name))

        self.__checkBuild(build, C.Db)
        self.__checkChrom(L, D, name)

        result = D.chromAcc(name)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing chromAccession(%s %s)" % (build,
                     name))

        del D,L,C
        return result
    #chromAccession

    @soapmethod(String, String, _returns = String)
    def chromosomeName(self, build, accNo) :
        """
            Get the name of a chromosome, given a chromosome accession number.

            Arguments:
                build   ; The human genome build.
                accNo    ; The accession number of a chromosome (NC_...).

            Returns:
                string ; The name of a chromosome.
        """
        C = Config.Config() # Read the configuration file.
        D = Db.Mapping(build, C.Db)
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request chromName(%s %s)" % (build, accNo))

        self.__checkBuild(build, C.Db)
#        self.__checkChrom(L, D, name)

        result = D.chromName(accNo)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing chromName(%s %s)" % (build,
                     accNo))

        del D,L,C
        return result
    #chromosomeName

    @soapmethod(String, String, _returns = String)
    def getchromName(self, build, acc) :
        """
            Get the chromosome name, given a transcript identifier (NM number).

            Arguments:
                build ; The human genome build.
                acc   ; The NM accession number (version NOT included)

            Returns:
                string ; The name of a chromosome.
        """
        C = Config.Config() # Read the configuration file.
        D = Db.Mapping(build, C.Db)
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request getchromName(%s %s)" % (build, acc))

        self.__checkBuild(build, C.Db)
#        self.__checkChrom(L, D, name)

        result = D.get_chromName(acc)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing getchromName(%s %s)" % (build,
                     acc))

        del D,L,C
        return result
    #chromosomeName

    @soapmethod(String, String, _returns = Array(String))
    def numberConversion(self, build, variant) :
        """
            Converts c. to g. notation or vice versa


            Arguments (all strings:
                build   ; The human genome build (hg19 or hg18).
                variant ; The variant in either c. or g. notation, full HGVS
                          notation, including NM_ or NC_ accession number.

            Returns:
                string; The variant in either g. or c. notation.

        """

        C = Config.Config() # Read the configuration file.
        D = Db.Mapping(build, C.Db)
        O = Output.Output(__file__, C.Output)
        O.addMessage(__file__, -1, "INFO",
                     "Received request cTogConversion(%s %s)" % (
                     build, variant))
        converter = Mapper.Converter(build, C, O)
        variant = converter.correctChrVariant(variant)

        if "c." in variant :
            result = [converter.c2chrom(variant)]
        elif "g." in variant :
            result = converter.chrom2c(variant)
        else:
            result = [""]

        O.addMessage(__file__, -1, "INFO",
                     "Finished processing cTogConversion(%s %s)" % (
                     build, variant))
        return result
    #numberConversion

    @soapmethod(String, _returns = String)
    def checkSyntax(self, variant) :
        """
        """
        C = Config.Config() # Read the configuration file.
        L = Output.Output(__file__, C.Output)
        L.addMessage(__file__, -1, "INFO",
                     "Received request checkSyntax(%s)" % (variant))

        self.__checkVariant(L, variant)

        P = Parser.Nomenclatureparser(L)
        parsetree = P.parse(variant)
        del C, P
        if not parsetree :
            result = "This variant does not have the right syntax. "\
                     "Please try again."
        #if
        else :
            result = "The syntax of this variant is OK!"

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing checkSyntax(%s)" % (variant))

        del L
        return result
    #checkSyntax
    
    @soapmethod(String, _returns = MutalyzerOutput)
    def runMutalyzer(self, variant) :
        C = Config.Config() # Read the configuration file.
        L = Output.Output(__file__, C.Output)
        L.addMessage(__file__, -1, "INFO",
                     "Received request runMutalyzer(%s)" % (variant))
        Mutalyzer.process(variant, C, L)                     
        M = MutalyzerOutput()
        M.original = L.getOutput("original")[0]
        M.mutated = L.getOutput("mutated")[0]
        L.addMessage(__file__, -1, "INFO",
                     "Finished processing runMutalyzer(%s)" % (variant))
        return M
    #runMutalyzer
#MutalyzerService
