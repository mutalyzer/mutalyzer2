#!/usr/bin/env python

"""
Mutalyzer webservices.

The SOAP webservice is exposed through a WSGI interface.

Example Apache/mod_wsgi configuration:

   WSGIScriptAlias /services /var/www/mutalyzer/src/webservice.py

Be sure to have this line first if you also define a / alias, like this:

   WSGIScriptAlias /services /var/www/mutalyzer/src/webservice.py
   WSGIScriptAlias / /var/www/mutalyzer/src/wsgi.py

@todo: Do we really use namespaces correctly?
@todo: For some reason, the server exposes its location including ?wsdl.
@todo: More thourough input checking. The @soap decorator does not do any
       kind of strictness checks on the input. For example, in
       transcriptInfo, the build argument must really be present. (Hint:
       use __checkBuild.)
"""

import logging; logging.basicConfig()

# We now use very current soaplib:
#   $ git clone https://github.com/soaplib/soaplib.git
#   $ cd soaplib
#   $ sudo python setup.py install

from soaplib.core import Application
from soaplib.core.service import soap
from soaplib.core.service import DefinitionBase
from soaplib.core.model.primitive import String, Integer
from soaplib.core.model.clazz import Array
from soaplib.core.model.exception import Fault
from soaplib.core.server import wsgi
import os
import site

# Add /src to Python path
site.addsitedir(os.path.dirname(__file__))

# Todo: fix Mutalyzer to not depend on working directory
os.chdir(os.path.split(os.path.dirname(__file__))[0])

import Mutalyzer
from Modules import Db
from Modules import Output
from Modules import Config
from Modules import Parser
from Modules import Mapper
from Modules import Retriever
from Modules.Serializers import Mapping, Transcript, MutalyzerOutput, Mandatory, TranscriptNameInfo, CheckSyntaxOutput, SoapMessage


class MutalyzerService(DefinitionBase) :
    """
    Mutalyzer webservices.

    These methods are made public via a SOAP interface.

    Note: Don't use leading newline in the docstrings of SOAP methods. These
    are visible in the generated documentation.
    """

    def __checkBuild(self, L, build, config) :
        """
        Check if the build is supported (hg18 or hg19).

        Returns:
            - Nothing (but raises an EARG exception).

        @arg L: An output object for logging.
        @type L: object
        @arg build: The human genome build name that needs to be checked.
        @type build: string
        @arg config: Configuration object of the Db module.
        @type config: object
        """

        if not build in config.dbNames :
            L.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                        "build name." % build)
        #if
    #__checkBuild

    def __checkChrom(self, L, D, chrom) :
        """
        Check if the chromosome is in our database.

        Returns:
            - Nothing (but raises an EARG exception).

        @arg L: An output object for logging.
        @type L: object
        @arg D: A handle to the database.
        @type D: object
        @arg chrom: The name of the chromosome.
        @type chrom: string
        """

        if not D.isChrom(chrom) :
            L.addMessage(__file__, 4, "EARG", "EARG %s" % chrom)
            raise Fault("EARG",
                        "The chrom argument (%s) was not a valid " \
                        "chromosome name." % chrom)
        #if
    #__checkChrom

    def __checkPos(self, L, pos) :
        """
        Check if the position is valid.

        Returns:
            - Nothing (but raises an ERANGE exception).

        @arg L: An output object for logging.
        @type L: object
        @arg pos: The position.
        @type pos: integer
        """

        if pos < 1 :
            L.addMessage(__file__, 4, "ERANGE", "ERANGE %i" % pos)
            raise Fault("ERANGE",
                        "The pos argument (%i) is out of range." % pos)
        #if
    #__checkPos

    def __checkVariant(self, L, variant) :
        """
        Check if a variant is provided.

        Returns:
            - Nothing (but raises an EARG exception).

        @arg L: An output object for logging.
        @type L: object
        @arg variant: The variant.
        @type variant: string
        """

        if not variant :
            L.addMessage(__file__, 4, "EARG", "EARG no variant")
            raise Fault("EARG", "The variant argument is not provided.")
        #if
    #__checkVariant

    @soap(Mandatory.String, Mandatory.String, Mandatory.Integer,
        _returns = Array(Mandatory.String))
    def getTranscripts(self, build, chrom, pos) :
        """Get all the transcripts that overlap with a chromosomal position.

        On error an exception is raised:
          - detail       ; Human readable description of the error.
          - faultstring: ; A code to indicate the type of error.
              - EARG   ; The argument was not valid.
              - ERANGE ; An invalid range was given.

        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg chrom: A chromosome encoded as "chr1", ..., "chrY".
        @type chrom: string
        @arg pos: A position on the chromosome.
        @type pos: integer

        @return: A list of transcripts.
        @rtype: list
        """

        C = Config.Config()
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request getTranscripts(%s %s %s)" % (build,
                     chrom, pos))

        self.__checkBuild(L, build, C.Db)
        D = Db.Mapping(build, C.Db)

        self.__checkChrom(L, D, chrom)
        self.__checkPos(L, pos)

        ret = D.get_Transcripts(chrom, pos, pos, True)

        #filter out the accNo
        ret = [r[0] for r in ret]


        L.addMessage(__file__, -1, "INFO",
                     "Finished processing getTranscripts(%s %s %s)" % (build,
                     chrom, pos))

        L.addMessage(__file__, -1, "INFO",
                     "We return %s" % ret)

        del D, L, C
        return ret
    #getTranscripts

    @soap(Mandatory.String, Mandatory.String, _returns = Array(Mandatory.String))
    def getTranscriptsByGeneName(self, build, name) :
        """Todo: documentation.
        """

        C = Config.Config()
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request getTranscriptsByGene(%s %s)" % (build,
                     name))

        self.__checkBuild(L, build, C.Db)
        D = Db.Mapping(build, C.Db)

        ret = D.get_TranscriptsByGeneName(name)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing getTranscriptsByGene(%s %s)" % (
                     build, name))
        return ret
    #getTranscriptsByGene

    @soap(Mandatory.String, Mandatory.String, Mandatory.Integer,
        Mandatory.Integer, Mandatory.Integer, _returns = Array(Mandatory.String))
    def getTranscriptsRange(self, build, chrom, pos1, pos2, method) :
        """Get all the transcripts that overlap with a range on a chromosome.

        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg chrom: A chromosome encoded as "chr1", ..., "chrY".
        @type chrom: string
        @arg pos1: The first postion of the range.
        @type pos1: integer
        @arg pos2: The last postion of the range.
        @type pos2: integer
        @arg method: The method of determining overlap:
            - 0 ; Return only the transcripts that completely fall in the range
                  [pos1, pos2].
            - 1 ; Return all hit transcripts.

        @return: A list of transcripts.
        @rtype: list
        """

        C = Config.Config()
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
            "Received request getTranscriptsRange(%s %s %s %s %s)" % (build,
            chrom, pos1, pos2, method))

        D = Db.Mapping(build, C.Db)
        self.__checkBuild(L, build, C.Db)

        ret = D.get_Transcripts(chrom, pos1, pos2, method)

        #filter out the accNo
        ret = [r[0] for r in ret]

        L.addMessage(__file__, -1, "INFO",
            "Finished processing getTranscriptsRange(%s %s %s %s %s)" % (
            build, chrom, pos1, pos2, method))

        del D, L, C
        return ret
    #getTranscriptsRange

    @soap(Mandatory.String, Mandatory.String, _returns = Mandatory.String)
    def getGeneName(self, build, accno) :
        """Find the gene name associated with a transcript.

        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg accno: The identifier of a transcript.
        @type accno: string

        @return: The name of the associated gene.
        @rtype: string
        """

        C = Config.Config()
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request getGeneName(%s %s)" % (build, accno))

        D = Db.Mapping(build, C.Db)
        self.__checkBuild(L, build, C.Db)

        ret = D.get_GeneName(accno.split('.')[0])

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing getGeneName(%s %s)" % (build, accno))

        del D, L, C
        return ret
    #getGeneName


    @soap(Mandatory.String, Mandatory.String, Mandatory.String,
        Mandatory.String, _returns = Mapping)
    def mappingInfo(self, LOVD_ver, build, accNo, variant) :
        """Search for an NM number in the MySQL database, if the version
        number matches, get the start and end positions in a variant and
        translate these positions to I{g.} notation if the variant is in I{c.}
        notation and vice versa.

          - If no end position is present, the start position is assumed to be
            the end position.
          - If the version number is not found in the database, an error message
            is generated and a suggestion for an other version is given.
          - If the reference sequence is not found at all, an error is returned.
          - If no variant is present, an error is returned.
          - If the variant is not accepted by the nomenclature parser, a parse
            error will be printed.

        @arg LOVD_ver: The LOVD version.
        @type LOVD_ver: string
        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg accNo: The NM accession number and version.
        @type accNo: string
        @arg variant: The variant.
        @type variant: string

        @return: Complex object:
          - start_main   ; The main coordinate of the start position
                           in I{c.} (non-star) notation.
          - start_offset ; The offset coordinate of the start position
                           in I{c.} notation (intronic position).
          - end_main     ; The main coordinate of the end position in
                           I{c.} (non-star) notation.
          - end_offset   ; The offset coordinate of the end position in
                           I{c.} notation (intronic position).
          - start_g      ; The I{g.} notation of the start position.
          - end_g        ; The I{g.} notation of the end position.
          - type         ; The mutation type.
        @rtype: object
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

    @soap(Mandatory.String, Mandatory.String, Mandatory.String,
        _returns = Transcript)
    def transcriptInfo(self, LOVD_ver, build, accNo) :
        """Search for an NM number in the MySQL database, if the version
        number matches, the transcription start and end and CDS end
        in I{c.} notation is returned.

        @arg LOVD_ver: The LOVD version.
        @type LOVD_ver: string
        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg accNo: The NM accession number and version.
        @type accNo: string

        @return: Complex object:
          - trans_start  ; Transcription start in I{c.} notation.
          - trans_stop   ; Transcription stop in I{c.} notation.
          - CDS_stop     ; CDS stop in I{c.} notation.
        @rtype: object
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

    @soap(Mandatory.String, Mandatory.String, _returns = Mandatory.String)
    def chromAccession(self, build, name) :
        """Get the accession number of a chromosome, given a name.

        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg name: The name of a chromosome (e.g. chr1).
        @type name: string

        @return: The accession number of a chromosome.
        @rtype: string
        """
        C = Config.Config() # Read the configuration file.
        D = Db.Mapping(build, C.Db)
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request chromAccession(%s %s)" % (build, name))

        self.__checkBuild(L, build, C.Db)
        self.__checkChrom(L, D, name)

        result = D.chromAcc(name)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing chromAccession(%s %s)" % (build,
                     name))

        del D,L,C
        return result
    #chromAccession

    @soap(Mandatory.String, Mandatory.String, _returns = Mandatory.String)
    def chromosomeName(self, build, accNo) :
        """Get the name of a chromosome, given a chromosome accession number.

        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg accNo: The accession number of a chromosome (NC_...).
        @type accNo: string

        @return: The name of a chromosome.
        @rtype: string
        """
        C = Config.Config() # Read the configuration file.
        D = Db.Mapping(build, C.Db)
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request chromName(%s %s)" % (build, accNo))

        self.__checkBuild(L, build, C.Db)
#        self.__checkChrom(L, D, name)

        result = D.chromName(accNo)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing chromName(%s %s)" % (build,
                     accNo))

        del D,L,C
        return result
    #chromosomeName

    @soap(Mandatory.String, Mandatory.String, _returns = Mandatory.String)
    def getchromName(self, build, acc) :
        """Get the chromosome name, given a transcript identifier (NM number).

        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg acc: The NM accession number (version NOT included).
        @type acc: string

        @return: The name of a chromosome.
        @rtype: string
        """
        C = Config.Config() # Read the configuration file.
        D = Db.Mapping(build, C.Db)
        L = Output.Output(__file__, C.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request getchromName(%s %s)" % (build, acc))

        self.__checkBuild(L, build, C.Db)
#        self.__checkChrom(L, D, name)

        result = D.get_chromName(acc)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing getchromName(%s %s)" % (build,
                     acc))

        del D,L,C
        return result
    #chromosomeName

    @soap(Mandatory.String, Mandatory.String, _returns = Array(Mandatory.String))
    def numberConversion(self, build, variant) :
        """Converts I{c.} to I{g.} notation or vice versa


        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg variant: The variant in either I{c.} or I{g.} notation, full HGVS
                      notation, including NM_ or NC_ accession number.
        @type variant: string

        @return: The variant(s) in either I{g.} or I{c.} notation.
        @rtype: list
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
            result = converter.chrom2c(variant, "list")
        else:
            result = [""]

        O.addMessage(__file__, -1, "INFO",
                     "Finished processing cTogConversion(%s %s)" % (
                     build, variant))
        return result
    #numberConversion

    @soap(Mandatory.String, _returns = CheckSyntaxOutput)
    def checkSyntax(self, variant):
        """Checks the syntax of a variant.

        @arg variant: The variant to check.
        @type variant: string

        @return: Object with fields:
                 - valid: A boolean indicating parse result (true for
                          succes, false in case of a parse error).
                 - messages: List of (error) messages as strings.
        @rtype: object
        """
        C = Config.Config() # Read the configuration file.
        O = Output.Output(__file__, C.Output)
        O.addMessage(__file__, -1, "INFO",
                     "Received request checkSyntax(%s)" % (variant))

        result = CheckSyntaxOutput()

        self.__checkVariant(O, variant)

        P = Parser.Nomenclatureparser(O)
        parsetree = P.parse(variant)
        del C, P
        result.valid = bool(parsetree)

        O.addMessage(__file__, -1, "INFO",
                     "Finished processing checkSyntax(%s)" % (variant))

        result.messages = O.getSoapMessages()

        del O
        return result
    #checkSyntax

    @soap(Mandatory.String, _returns = MutalyzerOutput)
    def runMutalyzer(self, variant) :
        """Todo: documentation.
        """
        C = Config.Config() # Read the configuration file.
        O = Output.Output(__file__, C.Output)
        O.addMessage(__file__, -1, "INFO",
                     "Received request runMutalyzer(%s)" % (variant))
        Mutalyzer.process(variant, C, O)

        result = MutalyzerOutput()

        # We force the results to strings here, because some results
        # may be of type Bio.Seq.Seq which soaplib doesn't like.
        #
        # todo: We might have to also do this elsewhere.

        result.original = str(O.getIndexedOutput("original", 0))
        result.mutated = str(O.getIndexedOutput("mutated", 0))

        result.origMRNA = str(O.getIndexedOutput("origMRNA", 0))
        result.mutatedMRNA = str(O.getIndexedOutput("mutatedMRNA", 0))

        result.origCDS = str(O.getIndexedOutput("origCDS", 0))
        result.newCDS = str(O.getIndexedOutput("newCDS", 0))

        result.origProtein = str(O.getIndexedOutput("oldprotein", 0))
        result.newProtein = str(O.getIndexedOutput("newprotein", 0))
        result.altProtein = str(O.getIndexedOutput("altProtein", 0))

        result.errors, result.warnings, result.summary = O.Summary()

        O.addMessage(__file__, -1, "INFO",
                     "Finished processing runMutalyzer(%s)" % (variant))

        result.messages = O.getSoapMessages()

        return result
    #runMutalyzer

    @soap(Mandatory.String, Mandatory.String, _returns = TranscriptNameInfo)
    def getGeneAndTranscript(self, genomicReference, transcriptReference) :
        """Todo: documentation.
        """
        C = Config.Config()
        O = Output.Output(__file__, C.Output)
        D = Db.Cache(C.Db)

        O.addMessage(__file__, -1, "INFO",
            "Received request getGeneAndTranscript(%s, %s)" % (genomicReference,
            transcriptReference))
        retriever = Retriever.GenBankRetriever(C.Retriever, O, D)
        record = retriever.loadrecord(genomicReference)

        ret = TranscriptNameInfo()
        for i in record.geneList :
            for j in i.transcriptList :
                if j.transcriptID == transcriptReference :
                    ret.transcriptName = "%s_v%s" % (i.name, j.name)
                    ret.productName = j.transcriptProduct
                #if

        O.addMessage(__file__, -1, "INFO",
            "Finished processing getGeneAndTranscript(%s, %s)" % (
            genomicReference, transcriptReference))
        return ret
    #getGeneAndTranscript

    @soap(Mandatory.String, _returns = Mandatory.String)
    def upLoadGenBankLocalFile(self, content) :
        pass
    #upLoadGenBankLocalFile

    @soap(Mandatory.String, _returns = Mandatory.String)
    def upLoadGenBankRemoteFile(self, url) :
        pass
    #upLoadGenBankLocalFile

    @soap(Mandatory.String, Mandatory.String, Mandatory.Integer,
        Mandatory.Integer, _returns = Mandatory.String)
    def sliceChromosomeByGene(self, geneSymbol, organism, upStream,
        downStream) :
        """
        Todo: documentation, error handling, argument checking, tests.
        """

        C = Config.Config()
        O = Output.Output(__file__, C.Output)
        D = Db.Cache(C.Db)
        retriever = Retriever.GenBankRetriever(C.Retriever, O, D)

        O.addMessage(__file__, -1, "INFO",
            "Received request sliceChromosomeByGene(%s, %s, %s, %s)" % (
            geneSymbol, organism, upStream, downStream))

        UD = retriever.retrievegene(geneSymbol, organism, upStream, downStream)

        O.addMessage(__file__, -1, "INFO",
            "Finished processing sliceChromosomeByGene(%s, %s, %s, %s)" % (
            geneSymbol, organism, upStream, downStream))

        return UD
    #sliceChromosomeByGene

    @soap(Mandatory.String, _returns = Mandatory.String)
    def sliceChromosome(self, chromAccNo, start, end, orientation) :
        pass
    #sliceChromosome
#MutalyzerService

# WSGI application for use with e.g. Apache/mod_wsgi
soap_application = Application([MutalyzerService],
                               'http://mutalyzer.nl/2.0/services', # namespace
                               'Mutalyzer')
application = wsgi.Application(soap_application)

# We can also use the built-in webserver by executing this file directly
if __name__ == '__main__':
    # Todo: Setting the working directory probably doesn't work
    from wsgiref.simple_server import make_server
    print 'Listening to http://localhost:8081/'
    print 'WDSL file is at http://localhost:8081/?wsdl'
    make_server('localhost', 8081, application).serve_forever()
