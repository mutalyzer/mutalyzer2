"""
Mutalyzer webservices.

@todo: Do we really use namespaces correctly?
@todo: For some reason, the server exposes its location including ?wsdl.
@todo: More thourough input checking. The @soap decorator does not do any
       kind of strictness checks on the input. For example, in
       transcriptInfo, the build argument must really be present. (Hint:
       use __checkBuild.)
"""


# WSGI applications should never print anything to stdout. We redirect to
# stderr, but eventually Mutalyzer should be fixed to never just 'print'
# anything.
# http://code.google.com/p/modwsgi/wiki/DebuggingTechniques
import sys
sys.stdout = sys.stderr

# Log exceptions to stdout
import logging; logging.basicConfig()

from soaplib.core import Application
from soaplib.core.service import soap
from soaplib.core.service import DefinitionBase
from soaplib.core.model.primitive import String, Integer, Boolean, DateTime
from soaplib.core.model.clazz import Array
from soaplib.core.model.exception import Fault
from soaplib.core.server import wsgi
import os
import socket
from operator import itemgetter, attrgetter

import mutalyzer
from mutalyzer.config import Config
from mutalyzer.output import Output
from mutalyzer.grammar import Grammar
from mutalyzer.sync import CacheSync
from mutalyzer import variantchecker
from mutalyzer import Db
from mutalyzer.mapping import Converter
from mutalyzer import Retriever
from mutalyzer import GenRecord
from mutalyzer.models import *


class MutalyzerService(DefinitionBase):
    """
    Mutalyzer webservices.

    These methods are made public via a SOAP interface.
    """
    def __init__(self, environ=None):
        self._config = Config()
        super(MutalyzerService, self).__init__(environ)
    #__init__

    def __checkBuild(self, L, build) :
        """
        Check if the build is supported (hg18 or hg19).

        Returns:
            - Nothing (but raises an EARG exception).

        @arg L: An output object for logging.
        @type L: object
        @arg build: The human genome build name that needs to be checked.
        @type build: string
        """

        if not build in self._config.Db.dbNames :
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

    @soap(Mandatory.String, Mandatory.String, Mandatory.Integer, Boolean,
        _returns = Array(Mandatory.String))
    def getTranscripts(self, build, chrom, pos, versions=False) :
        """
        Get all the transcripts that overlap with a chromosomal position.

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
        @type pos: int
        @kwarg versions: If set to True, also include transcript versions.
        @type versions: bool

        @return: A list of transcripts.
        @rtype: list
        """
        L = Output(__file__, self._config.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request getTranscripts(%s %s %s %s)" % (build,
                     chrom, pos, versions))

        self.__checkBuild(L, build)
        D = Db.Mapping(build, self._config.Db)

        self.__checkChrom(L, D, chrom)
        self.__checkPos(L, pos)

        ret = D.get_Transcripts(chrom, pos, pos, True)

        #filter out the accNo
        if versions:
            ret = [r[0] + '.' + str(r[-1]) for r in ret]
        else:
            ret = [r[0] for r in ret]

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing getTranscripts(%s %s %s %s)" % (build,
                     chrom, pos, versions))

        L.addMessage(__file__, -1, "INFO",
                     "We return %s" % ret)

        del D, L
        return ret
    #getTranscripts

    @soap(Mandatory.String, Mandatory.String, _returns = Array(Mandatory.String))
    def getTranscriptsByGeneName(self, build, name):
        """
        Todo: documentation.
        """
        L = Output(__file__, self._config.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request getTranscriptsByGene(%s %s)" % (build,
                     name))

        self.__checkBuild(L, build)
        D = Db.Mapping(build, self._config.Db)

        ret = D.get_TranscriptsByGeneName(name)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing getTranscriptsByGene(%s %s)" % (
                     build, name))

        return ret
    #getTranscriptsByGene

    @soap(Mandatory.String, Mandatory.String, Mandatory.Integer,
        Mandatory.Integer, Mandatory.Integer, _returns = Array(Mandatory.String))
    def getTranscriptsRange(self, build, chrom, pos1, pos2, method) :
        """
        Get all the transcripts that overlap with a range on a chromosome.

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
        L = Output(__file__, self._config.Output)

        L.addMessage(__file__, -1, "INFO",
            "Received request getTranscriptsRange(%s %s %s %s %s)" % (build,
            chrom, pos1, pos2, method))

        D = Db.Mapping(build, self._config.Db)
        self.__checkBuild(L, build)

        ret = D.get_Transcripts(chrom, pos1, pos2, method)

        #filter out the accNo
        ret = [r[0] for r in ret]

        L.addMessage(__file__, -1, "INFO",
            "Finished processing getTranscriptsRange(%s %s %s %s %s)" % (
            build, chrom, pos1, pos2, method))

        del D, L
        return ret
    #getTranscriptsRange

    @soap(Mandatory.String, Mandatory.String, _returns = Mandatory.String)
    def getGeneName(self, build, accno) :
        """
        Find the gene name associated with a transcript.

        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg accno: The identifier of a transcript.
        @type accno: string

        @return: The name of the associated gene.
        @rtype: string
        """
        L = Output(__file__, self._config.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request getGeneName(%s %s)" % (build, accno))

        D = Db.Mapping(build, self._config.Db)
        self.__checkBuild(L, build)

        ret = D.get_GeneName(accno.split('.')[0])

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing getGeneName(%s %s)" % (build, accno))

        del D, L
        return ret
    #getGeneName


    @soap(Mandatory.String, Mandatory.String, Mandatory.String,
        Mandatory.String, _returns = Mapping)
    def mappingInfo(self, LOVD_ver, build, accNo, variant) :
        """
        Search for an NM number in the MySQL database, if the version
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
        L = Output(__file__, self._config.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Reveived request mappingInfo(%s %s %s %s)" % (
                        LOVD_ver, build, accNo, variant))

        conv = Converter(build, self._config, L)
        result = conv.mainMapping(accNo, variant)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing mappingInfo(%s %s %s %s)" % (
                        LOVD_ver, build, accNo, variant))

        del L
        return result
    #mappingInfo

    @soap(Mandatory.String, Mandatory.String, Mandatory.String,
        _returns = Transcript)
    def transcriptInfo(self, LOVD_ver, build, accNo) :
        """
        Search for an NM number in the MySQL database, if the version
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
        O = Output(__file__, self._config.Output)

        O.addMessage(__file__, -1, "INFO",
                     "Received request transcriptInfo(%s %s %s)" % (LOVD_ver,
                     build, accNo))

        converter = Converter(build, self._config, O)
        T = converter.mainTranscript(accNo)

        O.addMessage(__file__, -1, "INFO",
                     "Finished processing transcriptInfo(%s %s %s)" % (
                     LOVD_ver, build, accNo))
        return T
    #transcriptInfo

    @soap(Mandatory.String, Mandatory.String, _returns = Mandatory.String)
    def chromAccession(self, build, name) :
        """
        Get the accession number of a chromosome, given a name.

        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg name: The name of a chromosome (e.g. chr1).
        @type name: string

        @return: The accession number of a chromosome.
        @rtype: string
        """
        D = Db.Mapping(build, self._config.Db)
        L = Output(__file__, self._config.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request chromAccession(%s %s)" % (build, name))

        self.__checkBuild(L, build)
        self.__checkChrom(L, D, name)

        result = D.chromAcc(name)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing chromAccession(%s %s)" % (build,
                     name))

        del D,L
        return result
    #chromAccession

    @soap(Mandatory.String, Mandatory.String, _returns = Mandatory.String)
    def chromosomeName(self, build, accNo) :
        """
        Get the name of a chromosome, given a chromosome accession number.

        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg accNo: The accession number of a chromosome (NC_...).
        @type accNo: string

        @return: The name of a chromosome.
        @rtype: string
        """
        D = Db.Mapping(build, self._config.Db)
        L = Output(__file__, self._config.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request chromName(%s %s)" % (build, accNo))

        self.__checkBuild(L, build)
#        self.__checkChrom(L, D, name)

        result = D.chromName(accNo)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing chromName(%s %s)" % (build,
                     accNo))

        del D,L
        return result
    #chromosomeName

    @soap(Mandatory.String, Mandatory.String, _returns = Mandatory.String)
    def getchromName(self, build, acc) :
        """
        Get the chromosome name, given a transcript identifier (NM number).

        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg acc: The NM accession number (version NOT included).
        @type acc: string

        @return: The name of a chromosome.
        @rtype: string
        """
        D = Db.Mapping(build, self._config.Db)
        L = Output(__file__, self._config.Output)

        L.addMessage(__file__, -1, "INFO",
                     "Received request getchromName(%s %s)" % (build, acc))

        self.__checkBuild(L, build)
#        self.__checkChrom(L, D, name)

        result = D.get_chromName(acc)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing getchromName(%s %s)" % (build,
                     acc))

        del D,L
        return result
    #chromosomeName

    @soap(Mandatory.String, Mandatory.String, _returns = Array(Mandatory.String))
    def numberConversion(self, build, variant) :
        """
        Converts I{c.} to I{g.} notation or vice versa


        @arg build: The human genome build (hg19 or hg18).
        @type build: string
        @arg variant: The variant in either I{c.} or I{g.} notation, full HGVS
                      notation, including NM_ or NC_ accession number.
        @type variant: string

        @return: The variant(s) in either I{g.} or I{c.} notation.
        @rtype: list
        """
        D = Db.Mapping(build, self._config.Db)
        O = Output(__file__, self._config.Output)
        O.addMessage(__file__, -1, "INFO",
                     "Received request cTogConversion(%s %s)" % (
                     build, variant))
        converter = Converter(build, self._config, O)
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
        """
        Checks the syntax of a variant.

        @arg variant: The variant to check.
        @type variant: string

        @return: Object with fields:
                 - valid: A boolean indicating parse result (true for
                          succes, false in case of a parse error).
                 - messages: List of (error) messages as strings.
        @rtype: object
        """
        output = Output(__file__, self._config.Output)
        output.addMessage(__file__, -1, "INFO",
                          "Received request checkSyntax(%s)" % (variant))

        result = CheckSyntaxOutput()

        self.__checkVariant(output, variant)

        grammar = Grammar(output)
        parsetree = grammar.parse(variant)
        result.valid = bool(parsetree)

        output.addMessage(__file__, -1, "INFO",
                          "Finished processing checkSyntax(%s)" % (variant))

        result.messages = []
        for message in output.getMessages():
            soap_message = SoapMessage()
            soap_message.errorcode = message.code
            soap_message.message = message.description
            result.messages.append(soap_message)

        return result
    #checkSyntax

    @soap(Mandatory.String, _returns = MutalyzerOutput)
    def runMutalyzer(self, variant) :
        """
        Todo: documentation.
        """
        O = Output(__file__, self._config.Output)
        O.addMessage(__file__, -1, "INFO",
                     "Received request runMutalyzer(%s)" % (variant))
        variantchecker.check_variant(variant, self._config, O)

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

        result.chromDescription = \
            O.getIndexedOutput("genomicChromDescription", 0)
        result.genomicDescription = \
            O.getIndexedOutput("genomicDescription", 0)
        result.transcriptDescriptions = O.getOutput("descriptions")
        result.proteinDescriptions = O.getOutput("protDescriptions")

        raw_variants = []
        for v in O.getOutput("visualisation"):
            r = RawVariant()
            r.description = v[0]
            r.visualisation = '%s\n%s' % (v[1], v[2])
            raw_variants.append(r)

        result.rawVariants = raw_variants

        result.errors, result.warnings, result.summary = O.Summary()

        O.addMessage(__file__, -1, "INFO",
                     "Finished processing runMutalyzer(%s)" % (variant))

        result.messages = []
        for message in O.getMessages():
            soap_message = SoapMessage()
            soap_message.errorcode = message.code
            soap_message.message = message.description
            result.messages.append(soap_message)

        return result
    #runMutalyzer

    @soap(Mandatory.String, Mandatory.String, _returns = TranscriptNameInfo)
    def getGeneAndTranscript(self, genomicReference, transcriptReference) :
        """
        Todo: documentation.
        """
        O = Output(__file__, self._config.Output)
        D = Db.Cache(self._config.Db)

        O.addMessage(__file__, -1, "INFO",
            "Received request getGeneAndTranscript(%s, %s)" % (genomicReference,
            transcriptReference))
        retriever = Retriever.GenBankRetriever(self._config.Retriever, O, D)
        record = retriever.loadrecord(genomicReference)

        GenRecordInstance = GenRecord.GenRecord(O, self._config.GenRecord)
        GenRecordInstance.record = record
        GenRecordInstance.checkRecord()

        ret = TranscriptNameInfo()
        for i in GenRecordInstance.record.geneList :
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

    @soap(Mandatory.String, String, _returns = Array(TranscriptInfo))
    def getTranscriptsAndInfo(self, genomicReference, geneName=None):
        """
        Given a genomic reference, return all its transcripts with their
        transcription/cds start/end sites and exons.

        @arg genomicReference: Name of a reference sequence.
        @type genomicReference: string

        @arg geneName: Name of gene to restrict returned transcripts to.
                       Default is to return all transcripts.
        @type geneName: string

        @return: Array of TranscriptInfo objects with fields:
                 - name
                 - id
                 - product
                 - cTransStart
                 - gTransStart
                 - cTransEnd
                 - gTransEnd
                 - sortableTransEnd
                 - cCDSStart
                 - gCDSStart
                 - cCDSStop
                 - gCDSStop
                 - locusTag
                 - linkMethod
                 - exons: Array of ExonInfo objects with fields:
                          - cStart
                          - gStart
                          - cStop
                          - gStop
                 - proteinTranscript: ProteinTranscript object with fields:
                                      - name
                                      - id
                                      - product
        """
        O = Output(__file__, self._config.Output)
        D = Db.Cache(self._config.Db)

        O.addMessage(__file__, -1, "INFO",
            "Received request getTranscriptsAndInfo(%s)" % genomicReference)
        retriever = Retriever.GenBankRetriever(self._config.Retriever, O, D)
        record = retriever.loadrecord(genomicReference)

        # Todo: If loadRecord failed (e.g. DTD missing), we should abort here.
        GenRecordInstance = GenRecord.GenRecord(O, self._config.GenRecord)
        GenRecordInstance.record = record
        GenRecordInstance.checkRecord()

        transcripts = []

        # The following loop is basically the same as building the legend in
        # the name checker web interface (website.Check).

        for gene in GenRecordInstance.record.geneList:
            # Only return transcripts for requested gene (if there was one)
            if geneName and gene.name != geneName:
                continue
            for transcript in sorted(gene.transcriptList,
                                     key=attrgetter('name')):

                # Exclude nameless transcripts
                if not transcript.name: continue

                t = TranscriptInfo()

                # Some raw info we don't use directly:
                # - transcript.CDS.location        CDS start and stop (g)
                # - transcript.CDS.positionList:   CDS splice sites (g) ?
                # - transcript.mRNA.location:      translation start and stop (g)
                # - transcript.mRNA.positionList:  splice sites (g)

                t.exons = []
                for i in range(0, transcript.CM.numberOfExons() * 2, 2):
                    exon = ExonInfo()
                    exon.gStart = transcript.CM.getSpliceSite(i)
                    exon.cStart = transcript.CM.g2c(exon.gStart)
                    exon.gStop = transcript.CM.getSpliceSite(i + 1)
                    exon.cStop = transcript.CM.g2c(exon.gStop)
                    t.exons.append(exon)

                # Beware that CM.info() gives a made-up value for trans_end,
                # which is sortable (no * notation). We therefore cannot use
                # it in our output and use the end position of the last exon
                # instead. The made-up value is still useful for sorting, so
                # we return it as sortableTransEnd.
                trans_start, sortable_trans_end, cds_stop = transcript.CM.info()
                cds_start = 1

                t.cTransEnd = str(t.exons[-1].cStop)
                t.gTransEnd = t.exons[-1].gStop
                t.sortableTransEnd = sortable_trans_end

                # Todo: If we have no CDS info, CM.info() gives trans_end as
                # value for cds_stop. This is an artifact to accomodate LOVD
                # stupidity an should probably be removed sometime.
                #if not transcript.CDS: cds_stop = None

                t.name = '%s_v%s' % (gene.name, transcript.name)
                t.id = transcript.transcriptID
                t.product = transcript.transcriptProduct
                t.cTransStart = str(trans_start)
                t.gTransStart = transcript.CM.x2g(trans_start, 0)
                t.cCDSStart = str(cds_start)
                t.gCDSStart = transcript.CM.x2g(cds_start, 0)
                t.cCDSStop = str(cds_stop)
                t.gCDSStop = transcript.CM.x2g(cds_stop, 0)
                t.locusTag = transcript.locusTag
                t.linkMethod = transcript.linkMethod

                t.proteinTranscript = None

                if transcript.translate:
                    p = ProteinTranscript()
                    p.name = '%s_i%s' % (gene.name, transcript.name)
                    p.id = transcript.proteinID
                    p.product = transcript.proteinProduct
                    t.proteinTranscript = p

                transcripts.append(t)

        O.addMessage(__file__, -1, "INFO",
            "Finished processing getTranscriptsAndInfo(%s)" % genomicReference)

        return transcripts
    #getTranscriptsAndInfo

    @soap(Mandatory.String, _returns = Mandatory.String)
    def upLoadGenBankLocalFile(self, content) :
        """
        Not implemented yet.
        """
        raise Fault('ENOTIMPLEMENTED', 'Not implemented yet')
    #upLoadGenBankLocalFile

    @soap(Mandatory.String, _returns = Mandatory.String)
    def upLoadGenBankRemoteFile(self, url) :
        """
        Not implemented yet.
        """
        raise Fault('ENOTIMPLEMENTED', 'Not implemented yet')
    #upLoadGenBankRemoteFile

    @soap(Mandatory.String, Mandatory.String, Mandatory.Integer,
        Mandatory.Integer, _returns = Mandatory.String)
    def sliceChromosomeByGene(self, geneSymbol, organism, upStream,
        downStream) :
        """
        Todo: documentation, error handling, argument checking, tests.
        """
        O = Output(__file__, self._config.Output)
        D = Db.Cache(self._config.Db)
        retriever = Retriever.GenBankRetriever(self._config.Retriever, O, D)

        O.addMessage(__file__, -1, "INFO",
            "Received request sliceChromosomeByGene(%s, %s, %s, %s)" % (
            geneSymbol, organism, upStream, downStream))

        UD = retriever.retrievegene(geneSymbol, organism, upStream, downStream)

        O.addMessage(__file__, -1, "INFO",
            "Finished processing sliceChromosomeByGene(%s, %s, %s, %s)" % (
            geneSymbol, organism, upStream, downStream))

        # Todo: use SOAP Fault object here (see Trac issue #41).
        if not UD:
            error = 'The request could not be completed\n' \
                    + '\n'.join(map(lambda m: str(m), O.getMessages()))
            raise Exception(error)

        return UD
    #sliceChromosomeByGene

    @soap(Mandatory.String, Mandatory.Integer, Mandatory.Integer,
        Mandatory.Integer, _returns = Mandatory.String)
    def sliceChromosome(self, chromAccNo, start, end, orientation) :
        """
        Todo: documentation, error handling, argument checking, tests.
        """
        O = Output(__file__, self._config.Output)
        D = Db.Cache(self._config.Db)
        retriever = Retriever.GenBankRetriever(self._config.Retriever, O, D)

        O.addMessage(__file__, -1, "INFO",
            "Received request sliceChromosome(%s, %s, %s, %s)" % (
            chromAccNo, start, end, orientation))

        UD = retriever.retrieveslice(chromAccNo, start, end, orientation)

        O.addMessage(__file__, -1, "INFO",
            "Finished processing sliceChromosome(%s, %s, %s, %s)" % (
            chromAccNo, start, end, orientation))

        return UD
    #sliceChromosome

    @soap(_returns = InfoOutput)
    def info(self):
        """
        Gives some static application information, such as the current running
        version.

        @return: Object with fields:
            - version: A string of the current running version.
            - versionParts: The parts of the current running version as a list
                of strings.
            - releaseDate: The release date for the running version as a
                string, or the empty string in case of a development version.
            - nomenclatureVersion: Version of the HGVS nomenclature used.
            - nomenclatureVersionParts: The parts of the HGVS nomenclature
                version as a list of strings.
            - serverName: The name of the server that is being queried.
            - contactEmail: The email address to contact for more information.
        @rtype: object
        """
        output = Output(__file__, self._config.Output)
        output.addMessage(__file__, -1, 'INFO', 'Received request info')

        result = InfoOutput()
        result.version = mutalyzer.__version__
        result.versionParts = mutalyzer.__version_info__
        if mutalyzer.RELEASE:
            result.releaseDate = mutalyzer.__date__
        else:
            result.releaseDate = ''
        result.nomenclatureVersion = mutalyzer.NOMENCLATURE_VERSION
        result.nomenclatureVersionParts = mutalyzer.NOMENCLATURE_VERSION_INFO
        result.serverName = socket.gethostname()
        result.contactEmail = mutalyzer.__contact__

        output.addMessage(__file__, -1, 'INFO', 'Finished processing info')
        return result
    #info

    @soap(DateTime, _returns = Array(CacheEntry))
    def getCache(self, created_since=None):
        """
        Get a list of entries from the local cache created since given date.

        This method is intended to be used by Mutalyzer itself to synchronize
        the cache between installations on different servers.
        """
        output = Output(__file__, self._config.Output)

        output.addMessage(__file__, -1, 'INFO',
                          'Received request getCache')

        database = Db.Cache(self._config.Db)
        sync = CacheSync(self._config.Retriever, output, database)

        cache = sync.local_cache(created_since)

        def cache_entry_to_soap(entry):
            e = CacheEntry()
            for attr in ('name', 'gi', 'hash', 'chromosomeName',
                         'chromosomeStart', 'chromosomeStop',
                         'chromosomeOrientation', 'url', 'created', 'cached'):
                setattr(e, attr, entry[attr])
            return e

        output.addMessage(__file__, -1, 'INFO',
                          'Finished processing getCache')

        return map(cache_entry_to_soap, cache)
    #getCache

    @soap(Mandatory.String, _returns = Array(Mandatory.String))
    def getdbSNPDescriptions(self, rs_id):
        """
        Lookup HGVS descriptions for a dbSNP rs identifier.

        @arg rs_id: The dbSNP rs identifier, e.g. 'rs9919552'.
        @type rs_id: string

        @return: List of HGVS descriptions.
        @rtype: list(string)
        """
        output = Output(__file__, self._config.Output)

        output.addMessage(__file__, -1, 'INFO',
                          'Received request getdbSNPDescription(%s)' % rs_id)

        retriever = Retriever.Retriever(self._config.Retriever, output, None)
        descriptions = retriever.snpConvert(rs_id)

        output.addMessage(__file__, -1, 'INFO',
                          'Finished processing getdbSNPDescription(%s)' % rs_id)

        # Todo: use SOAP Fault object here (see Trac issue #41).
        messages = output.getMessages()
        if messages:
            error = 'The request could not be completed\n' \
                    + '\n'.join(map(lambda m: str(m), output.getMessages()))
            raise Exception(error)

        return descriptions
    #getdbSNPDescriptions
#MutalyzerService


# WSGI application for use with e.g. Apache/mod_wsgi
soap_application = Application([MutalyzerService], mutalyzer.SOAP_NAMESPACE,
                               'Mutalyzer')
# Note: We would like to create the wsgi.Application instance only in the
# bin/mutalyzer-webservice.wsgi script, but unfortunately this breaks the
# get_wsdl method of soap_application which we use to generate API
# documentation in website.py.
application = wsgi.Application(soap_application)
