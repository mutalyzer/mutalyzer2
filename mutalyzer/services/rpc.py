"""
Mutalyzer RPC services.

@todo: More thourough input checking. The @soap decorator does not do any
       kind of strictness checks on the input. For example, in
       transcriptInfo, the build argument must really be present. (Hint:
       use __checkBuild.)
       We should use the built-in validator functionality from Spyne for
       this.
"""


from spyne.decorator import srpc
from spyne.service import ServiceBase
from spyne.model.primitive import String, Integer, Boolean, DateTime
from spyne.model.complex import Array
from spyne.model.fault import Fault
import os
import socket
import tempfile
from operator import itemgetter, attrgetter

import mutalyzer
from mutalyzer.config import settings
from mutalyzer.db import session
from mutalyzer.db.models import BatchJob, BatchQueueItem
from mutalyzer.output import Output
from mutalyzer.grammar import Grammar
from mutalyzer.sync import CacheSync
from mutalyzer import variantchecker
from mutalyzer import Db
from mutalyzer.mapping import Converter
from mutalyzer import File
from mutalyzer import Retriever
from mutalyzer import GenRecord
from mutalyzer import Scheduler
from mutalyzer.models import *
from mutalyzer import describe


def _checkBuild(L, build) :
    """
    Check if the build is supported (hg18, hg19, or mm10).

    Returns:
        - Nothing (but raises an EARG exception).

    @arg L: An output object for logging.
    @type L: object
    @arg build: The human genome build name that needs to be checked.
    @type build: string
    """

    if not build in settings.DB_NAMES:
        L.addMessage(__file__, 4, "EARG", "EARG %s" % build)
        raise Fault("EARG",
                    "The build argument (%s) was not a valid " \
                    "build name." % build)
    #if
#_checkBuild


def _checkChrom(L, D, chrom) :
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
        raise Fault("EARG", "The chrom argument (%s) was not a valid " \
            "chromosome name." % chrom)
    #if
#_checkChrom


def _checkPos(L, pos) :
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
        raise Fault("ERANGE", "The pos argument (%i) is out of range." % pos)
    #if
#_checkPos


def _checkVariant(L, variant) :
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
#_checkVariant


class MutalyzerService(ServiceBase):
    """
    Mutalyzer web services.

    These methods are made public via a SOAP interface.
    """
    def __init__(self, environ=None):
        super(MutalyzerService, self).__init__(environ)
    #__init__

    @srpc(Mandatory.ByteArray, String, String,  _returns=String)
    def submitBatchJob(data, process='NameChecker', argument=''):
        """
        Submit a batch job.

        Input and output file formats for batch jobs are explained on the
        website <https://mutalyzer.nl/batch>.

        On error an exception is raised:
          - detail: Human readable description of the error.
          - faultstring: A code to indicate the type of error.
              - EPARSE: The batch input could not be parsed.
              - EMAXSIZE: Input file exceeds maximum size.

        @arg data: Input file.
        @arg process: Optional type of the batch job, choose from: NameChecker
            (default), SyntaxChecker, PositionConverter, SnpConverter.
        @arg argument: Additional argument. Currently only used if batch_type
            is PositionConverter, denoting the human genome build.

        @return: Batch job identifier.
        """
        output = Output(__file__)

        counter = Db.Counter()
        counter.increment('batchjob', 'webservice')

        scheduler = Scheduler.Scheduler()
        file_instance = File.File(output)

        batch_types = ['NameChecker', 'SyntaxChecker', 'PositionConverter',
                       'SnpConverter']

        if process not in batch_types:
            raise Fault('EARG',
                        'The process argument must be one of %s.'
                        % ', '.join(batch_types))

        # Note that the max file size check below might be bogus, since Spyne
        # first checks the total request size, which by default has a maximum
        # of 2 megabytes.
        # In that case, a senv:Client.RequestTooLong faultstring is returned.

        # Todo: Set maximum request size by specifying the max_content_length
        #     argument for spyne.server.wsgi.WsgiApplication in all webservice
        #     instantiations.

        max_size = settings.MAX_FILE_SIZE

        batch_file = tempfile.TemporaryFile()
        size = 0
        try:
            for chunk in data:
                size += len(chunk)
                if size > max_size:
                    raise Fault('EMAXSIZE',
                                'Only files up to %s megabytes are accepted.' % (float(max_size) / 1048576))
                batch_file.write(chunk)
            batch_file.seek(0)
            job, columns = file_instance.parseBatchFile(batch_file)
        finally:
            try:
                batch_file.close()
            except IOError:
                pass

        if job is None:
            raise Fault('EPARSE', 'Could not parse input file, please check your file format.')

        result_id = scheduler.addJob('job@webservice', job, columns,
                                     'webservice', process, argument)
        return result_id

    @srpc(Mandatory.String, _returns=Integer)
    def monitorBatchJob(job_id):
        """
        Get the number of entries left for a batch job.

        Input and output file formats for batch jobs are explained on the
        website <https://mutalyzer.nl/batch>.

        @arg job_id: Batch job identifier.

        @return: Number of entries left.
        """
        return BatchQueueItem.query.join(BatchJob).filter_by(result_id=job_id).count()

    @srpc(Mandatory.String, _returns=ByteArray)
    def getBatchJob(job_id):
        """
        Get the result of a batch job.

        Input and output file formats for batch jobs are explained on the
        website <https://mutalyzer.nl/batch>.

        On error an exception is raised:
          - detail: Human readable description of the error.
          - faultstring: A code to indicate the type of error.
              - EBATCHNOTREADY: The batch job result is not yet ready.

        @arg job_id: Batch job identifier.

        @return: Batch job result file.
        """
        left = BatchQueueItem.query.join(BatchJob).filter_by(result_id=job_id).count()

        if left > 0:
            raise Fault('EBATCHNOTREADY', 'Batch job result is not yet ready.')

        filename = 'Results_%s.txt' % job_id
        handle = open(os.path.join(settings.CACHE_DIR, filename))
        return handle

    @srpc(Mandatory.String, Mandatory.String, Mandatory.Integer, Boolean,
        _returns=Array(Mandatory.String))
    def getTranscripts(build, chrom, pos, versions=False) :
        """
        Get all the transcripts that overlap with a chromosomal position.

        On error an exception is raised:
          - detail       ; Human readable description of the error.
          - faultstring: ; A code to indicate the type of error.
              - EARG   ; The argument was not valid.
              - ERANGE ; An invalid range was given.

        @arg build: The genome build (hg19, hg18, mm10).
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
        L = Output(__file__)

        L.addMessage(__file__, -1, "INFO",
            "Received request getTranscripts(%s %s %s %s)" % (build, chrom,
            pos, versions))

        _checkBuild(L, build)
        D = Db.Mapping(build)

        _checkChrom(L, D, chrom)
        _checkPos(L, pos)

        ret = D.get_Transcripts(chrom, pos, pos, True)

        #filter out the accNo
        if versions:
            ret = [r[0] + '.' + str(r[-1]) for r in ret]
        else:
            ret = [r[0] for r in ret]

        L.addMessage(__file__, -1, "INFO",
            "Finished processing getTranscripts(%s %s %s %s)" % (build, chrom,
            pos, versions))

        L.addMessage(__file__, -1, "INFO", "We return %s" % ret)

        del D, L
        return ret
    #getTranscripts

    @srpc(Mandatory.String, Mandatory.String, _returns=Array(Mandatory.String))
    def getTranscriptsByGeneName(build, name):
        """
        Todo: documentation.
        """
        L = Output(__file__)

        L.addMessage(__file__, -1, "INFO",
            "Received request getTranscriptsByGene(%s %s)" % (build, name))

        _checkBuild(L, build)
        D = Db.Mapping(build)

        ret = D.get_TranscriptsByGeneName(name)

        L.addMessage(__file__, -1, "INFO",
            "Finished processing getTranscriptsByGene(%s %s)" % (build, name))

        if ret :
            l = []
            for i in ret :
                l.append(i[0] + '.' + str(i[13]))
            return l

        return []
    #getTranscriptsByGene

    @srpc(Mandatory.String, Mandatory.String, Mandatory.Integer,
        Mandatory.Integer, Mandatory.Integer, _returns=Array(Mandatory.String))
    def getTranscriptsRange(build, chrom, pos1, pos2, method) :
        """
        Get all the transcripts that overlap with a range on a chromosome.

        @arg build: The genome build (hg19, hg18, mm10).
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
        L = Output(__file__)

        L.addMessage(__file__, -1, "INFO",
            "Received request getTranscriptsRange(%s %s %s %s %s)" % (build,
            chrom, pos1, pos2, method))

        D = Db.Mapping(build)
        _checkBuild(L, build)

        ret = D.get_Transcripts(chrom, pos1, pos2, method)

        #filter out the accNo
        ret = [r[0] for r in ret]

        L.addMessage(__file__, -1, "INFO",
            "Finished processing getTranscriptsRange(%s %s %s %s %s)" % (
            build, chrom, pos1, pos2, method))

        del D, L
        return ret
    #getTranscriptsRange

    @srpc(Mandatory.String, Mandatory.String, Mandatory.Integer,
        Mandatory.Integer, Mandatory.Integer,
        _returns=Array(TranscriptMappingInfo))
    def getTranscriptsMapping(build, chrom, pos1, pos2, method):
        """
        Get all the transcripts and their info that overlap with a range on a
        chromosome.

        @arg build: The genome build (hg19, hg18, mm10).
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

        @return: Array of TranscriptMappingInfo objects with fields:
                 - name
                 - version
                 - gene
                 - protein
                 - orientation
                 - start
                 - stop
                 - cds_start
                 - cds_stop
        """
        output = Output(__file__)
        output.addMessage(__file__, -1, 'INFO', 'Received request ' \
            'getTranscriptsRange(%s %s %s %s %s)' % (build, chrom, pos1, pos2,
            method))

        _checkBuild(output, build)

        database = Db.Mapping(build)
        transcripts = []

        for transcript in database.get_Transcripts(chrom, pos1, pos2, method):
            t = TranscriptMappingInfo()
            d = dict(zip(('transcript', 'selector', 'selector_version',
                'start', 'stop', 'cds_start', 'cds_stop', 'exon_starts',
                'exon_stops', 'gene', 'chromosome', 'orientation', 'protein',
                'version'), transcript))
            if d['orientation'] == '-':
                d['start'], d['stop'] = d['stop'], d['start']
                d['cds_start'], d['cds_stop'] = d['cds_stop'], d['cds_start']
            t.name = d['transcript']
            t.version = d['version']
            t.gene = d['gene']
            t.protein = d['protein']
            t.orientation = d['orientation']
            t.start = d['start']
            t.stop = d['stop']
            t.cds_start = d['cds_start']
            t.cds_stop = d['cds_stop']
            transcripts.append(t)

        output.addMessage(__file__, -1, 'INFO', 'Finished processing ' \
            'getTranscriptsRange(%s %s %s %s %s)' % (build, chrom, pos1, pos2,
            method))

        return transcripts
    #getTranscriptsMapping

    @srpc(Mandatory.String, Mandatory.String, _returns=Mandatory.String)
    def getGeneName(build, accno) :
        """
        Find the gene name associated with a transcript.

        @arg build: The genome build (hg19, hg18, mm10).
        @type build: string
        @arg accno: The identifier of a transcript.
        @type accno: string

        @return: The name of the associated gene.
        @rtype: string
        """
        L = Output(__file__)

        L.addMessage(__file__, -1, "INFO",
            "Received request getGeneName(%s %s)" % (build, accno))

        D = Db.Mapping(build)
        _checkBuild(L, build)

        ret = D.get_GeneName(accno.split('.')[0])

        L.addMessage(__file__, -1, "INFO",
            "Finished processing getGeneName(%s %s)" % (build, accno))

        del D, L
        return ret
    #getGeneName


    @srpc(Mandatory.String, Mandatory.String, Mandatory.String,
        Mandatory.String, _returns=Mapping)
    def mappingInfo(LOVD_ver, build, accNo, variant) :
        """
        Search for an NM number in the MySQL database, if the version
        number matches, get the start and end positions in a variant and
        translate these positions to I{g.} notation if the variant is in I{c.}
        notation and vice versa.

          - If no end position is present, the start position is assumed to be
            the end position.
          - If the version number is not found in the database, an error
            message is generated and a suggestion for an other version is
            given.
          - If the reference sequence is not found at all, an error is
            returned.
          - If no variant is present, an error is returned.
          - If the variant is not accepted by the nomenclature parser, a parse
            error will be printed.

        @arg LOVD_ver: The LOVD version.
        @type LOVD_ver: string
        @arg build: The genome build (hg19, hg18, mm10).
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
        L = Output(__file__)

        L.addMessage(__file__, -1, "INFO",
            "Reveived request mappingInfo(%s %s %s %s)" % (LOVD_ver, build,
            accNo, variant))

        conv = Converter(build, L)
        result = conv.mainMapping(accNo, variant)

        L.addMessage(__file__, -1, "INFO",
            "Finished processing mappingInfo(%s %s %s %s)" % (LOVD_ver, build,
            accNo, variant))

        del L
        return result
    #mappingInfo

    @srpc(Mandatory.String, Mandatory.String, Mandatory.String,
        _returns=Transcript)
    def transcriptInfo(LOVD_ver, build, accNo) :
        """
        Search for an NM number in the MySQL database, if the version
        number matches, the transcription start and end and CDS end
        in I{c.} notation is returned.

        @arg LOVD_ver: The LOVD version.
        @type LOVD_ver: string
        @arg build: The genome build (hg19, hg18, mm10).
        @type build: string
        @arg accNo: The NM accession number and version.
        @type accNo: string

        @return: Complex object:
          - trans_start  ; Transcription start in I{c.} notation.
          - trans_stop   ; Transcription stop in I{c.} notation.
          - CDS_stop     ; CDS stop in I{c.} notation.
        @rtype: object
        """
        O = Output(__file__)

        O.addMessage(__file__, -1, "INFO",
            "Received request transcriptInfo(%s %s %s)" % (LOVD_ver, build,
            accNo))

        converter = Converter(build, O)
        T = converter.mainTranscript(accNo)

        O.addMessage(__file__, -1, "INFO",
            "Finished processing transcriptInfo(%s %s %s)" % (LOVD_ver, build,
            accNo))
        return T
    #transcriptInfo

    @srpc(Mandatory.String, Mandatory.String, _returns=Mandatory.String)
    def chromAccession(build, name) :
        """
        Get the accession number of a chromosome, given a name.

        @arg build: The genome build (hg19, hg18, mm10).
        @type build: string
        @arg name: The name of a chromosome (e.g. chr1).
        @type name: string

        @return: The accession number of a chromosome.
        @rtype: string
        """
        D = Db.Mapping(build)
        L = Output(__file__)

        L.addMessage(__file__, -1, "INFO",
            "Received request chromAccession(%s %s)" % (build, name))

        _checkBuild(L, build)
        _checkChrom(L, D, name)

        result = D.chromAcc(name)

        L.addMessage(__file__, -1, "INFO",
            "Finished processing chromAccession(%s %s)" % (build, name))

        del D,L
        return result[0]
    #chromAccession

    @srpc(Mandatory.String, Mandatory.String, _returns=Mandatory.String)
    def chromosomeName(build, accNo) :
        """
        Get the name of a chromosome, given a chromosome accession number.

        @arg build: The genome build (hg19, hg18, mm10).
        @type build: string
        @arg accNo: The accession number of a chromosome (NC_...).
        @type accNo: string

        @return: The name of a chromosome.
        @rtype: string
        """
        D = Db.Mapping(build)
        L = Output(__file__)

        L.addMessage(__file__, -1, "INFO",
            "Received request chromName(%s %s)" % (build, accNo))

        _checkBuild(L, build)
#        self._checkChrom(L, D, name)

        result = D.chromName(accNo)

        L.addMessage(__file__, -1, "INFO",
            "Finished processing chromName(%s %s)" % (build, accNo))

        del D,L
        return result
    #chromosomeName

    @srpc(Mandatory.String, Mandatory.String, _returns=Mandatory.String)
    def getchromName(build, acc) :
        """
        Get the chromosome name, given a transcript identifier (NM number).

        @arg build: The genome build (hg19, hg18, mm10).
        @type build: string
        @arg acc: The NM accession number (version NOT included).
        @type acc: string

        @return: The name of a chromosome.
        @rtype: string
        """
        D = Db.Mapping(build)
        L = Output(__file__)

        L.addMessage(__file__, -1, "INFO",
            "Received request getchromName(%s %s)" % (build, acc))

        _checkBuild(L, build)
#        self._checkChrom(L, D, name)

        result = D.get_chromName(acc)

        L.addMessage(__file__, -1, "INFO",
            "Finished processing getchromName(%s %s)" % (build, acc))

        del D,L
        return result
    #chromosomeName

    @srpc(Mandatory.String, Mandatory.String, String,
        _returns=Array(Mandatory.String))
    def numberConversion(build, variant, gene=None):
        """
        Converts I{c.} to I{g.} notation or vice versa

        @arg build: The genome build (hg19, hg18, mm10).
        @type build: string
        @arg variant: The variant in either I{c.} or I{g.} notation, full HGVS
            notation, including NM_ or NC_ accession number.
        @type variant: string
        @kwarg gene: Optional gene name. If given, return variant descriptions
            on all transcripts for this gene.
        @type gene: string

        @return: The variant(s) in either I{g.} or I{c.} notation.
        @rtype: list
        """
        D = Db.Mapping(build)
        O = Output(__file__)
        O.addMessage(__file__, -1, "INFO",
            "Received request cTogConversion(%s %s)" % (build, variant))

        counter = Db.Counter()
        counter.increment('positionconvert', 'webservice')

        converter = Converter(build, O)
        variant = converter.correctChrVariant(variant)

        if "c." in variant or "n." in variant:
            result = [converter.c2chrom(variant)]
        elif "g." in variant or "m." in variant:
            result = converter.chrom2c(variant, "list", gene=gene)
        else:
            result = [""]

        O.addMessage(__file__, -1, "INFO",
            "Finished processing cTogConversion(%s %s)" % (build, variant))
        return result
    #numberConversion

    @srpc(Mandatory.String, _returns=CheckSyntaxOutput)
    def checkSyntax(variant):
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
        output = Output(__file__)
        output.addMessage(__file__, -1, "INFO",
            "Received request checkSyntax(%s)" % (variant))

        counter = Db.Counter()
        counter.increment('checksyntax', 'webservice')

        result = CheckSyntaxOutput()

        _checkVariant(output, variant)

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

    @srpc(Mandatory.String, _returns=MutalyzerOutput)
    def runMutalyzer(variant) :
        """
        Run the Mutalyzer name checker.

        @arg variant: The variant description to check.
        @type variant: string

        @return: Object with fields:
            - referenceId: Identifier of the reference sequence used.
            - sourceId: Identifier of the reference sequence source, e.g. the
                chromosomal accession number and version in case referenceId
                is a  UD reference created as a chromosomal slice.
            - sourceAccession: Accession number of the reference sequence
                source (only for genbank references).
            - sourceVersion: Version number of the reference sequence source
                (only for genbank references).
            - sourceGi: GI number of the reference sequence source (only for
                genbank references).
            - molecule: Molecular type of the reference sequence.
            - original: Original sequence.
            - mutated: Mutated sequence.
            - origMRNA: Original transcript sequence.
            - mutatedMRNA: Mutated transcript sequence.
            - origCDS: Original CDS.
            - newCDS: Mutated CDS.
            - origProtein: Original protein sequence.
            - newProtein: Mutated protein sequence.
            - altProtein: Alternative mutated protein sequence.
            - errors: Number of errors.
            - warnings: Number of warnings.
            - summary: Summary of messages.
            - chromDescription: Chromosomal description.
            - genomicDescription: Genomic description.
            - transcriptDescriptions: List of transcript descriptions.
            - proteinDescriptions: List of protein descriptions.
            - rawVariants: List of raw variants where each raw variant is
                represented by an object with fields:
                - description: Description of the raw variant.
                - visualisation: ASCII visualisation of the raw variant.
            - exons: If a transcript is selected, array of ExonInfo objects
                for each exon in the selected transcript with fields:
                - cStart
                - gStart
                - cStop
                - gStop
            - messages: List of (error) messages.
        """
        O = Output(__file__)
        O.addMessage(__file__, -1, "INFO",
            "Received request runMutalyzer(%s)" % (variant))

        counter = Db.Counter()
        counter.increment('namecheck', 'webservice')

        variantchecker.check_variant(variant, O)

        result = MutalyzerOutput()

        result.referenceId = O.getIndexedOutput('reference_id', 0)
        result.sourceId = O.getIndexedOutput('source_id', 0)
        result.sourceAccession = O.getIndexedOutput('source_accession', 0)
        result.sourceVersion = O.getIndexedOutput('source_version', 0)
        result.sourceGi = O.getIndexedOutput('source_gi', 0)
        result.molecule = O.getIndexedOutput('molecule', 0)

        # We force the results to strings here, because some results
        # may be of type Bio.Seq.Seq which spyne doesn't like.
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

        if O.getIndexedOutput('hasTranscriptInfo', 0, False):
            result.exons = []
            for e in O.getOutput('exonInfo'):
                exon = ExonInfo()
                exon.gStart, exon.gStop, exon.cStart, exon.cStop = e
                result.exons.append(exon)

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

    @srpc(Mandatory.String, Mandatory.String, _returns=TranscriptNameInfo)
    def getGeneAndTranscript(genomicReference, transcriptReference) :
        """
        Todo: documentation.
        """
        O = Output(__file__)
        D = Db.Cache()

        O.addMessage(__file__, -1, "INFO",
            "Received request getGeneAndTranscript(%s, %s)" % (
            genomicReference, transcriptReference))
        retriever = Retriever.GenBankRetriever(O, D)
        record = retriever.loadrecord(genomicReference)

        GenRecordInstance = GenRecord.GenRecord(O)
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

    @srpc(Mandatory.String, String, _returns=Array(TranscriptInfo))
    def getTranscriptsAndInfo(genomicReference, geneName=None):
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
                 - chromTransStart
                 - cTransEnd
                 - gTransEnd
                 - chromTransEnd
                 - sortableTransEnd
                 - cCDSStart
                 - gCDSStart
                 - chromCDSStart
                 - cCDSStop
                 - gCDSStop
                 - chromCDSStop
                 - locusTag
                 - linkMethod
                 - exons: Array of ExonInfo objects with fields:
                          - cStart
                          - gStart
                          - chromStart
                          - cStop
                          - gStop
                          - chromStop
                 - proteinTranscript: ProteinTranscript object with fields:
                                      - name
                                      - id
                                      - product
        """
        O = Output(__file__)
        D = Db.Cache()

        O.addMessage(__file__, -1, "INFO",
            "Received request getTranscriptsAndInfo(%s, %s)" % (
            genomicReference, geneName))
        retriever = Retriever.GenBankRetriever(O, D)
        record = retriever.loadrecord(genomicReference)

        # Todo: If loadRecord failed (e.g. DTD missing), we should abort here.
        GenRecordInstance = GenRecord.GenRecord(O)
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
                # - transcript.mRNA.location:      translation start and stop
                #                                  (g)
                # - transcript.mRNA.positionList:  splice sites (g)

                t.exons = []
                for i in range(0, transcript.CM.numberOfExons() * 2, 2):
                    exon = ExonInfo()
                    exon.gStart = transcript.CM.getSpliceSite(i)
                    exon.cStart = transcript.CM.g2c(exon.gStart)
                    exon.chromStart = GenRecordInstance.record.toChromPos(
                        exon.gStart)
                    exon.gStop = transcript.CM.getSpliceSite(i + 1)
                    exon.cStop = transcript.CM.g2c(exon.gStop)
                    exon.chromStop = GenRecordInstance.record.toChromPos(
                        exon.gStop)
                    t.exons.append(exon)

                # Beware that CM.info() gives a made-up value for trans_end,
                # which is sortable (no * notation). We therefore cannot use
                # it in our output and use the end position of the last exon
                # instead. The made-up value is still useful for sorting, so
                # we return it as sortableTransEnd.
                trans_start, sortable_trans_end, cds_stop = \
                    transcript.CM.info()
                cds_start = 1

                t.cTransEnd = str(t.exons[-1].cStop)
                t.gTransEnd = t.exons[-1].gStop
                t.chromTransEnd = GenRecordInstance.record.toChromPos(
                    t.gTransEnd)
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
                t.chromTransStart = GenRecordInstance.record.toChromPos(
                    t.gTransStart)
                t.cCDSStart = str(cds_start)
                t.gCDSStart = transcript.CM.x2g(cds_start, 0)
                t.chromCDSStart = GenRecordInstance.record.toChromPos(
                    t.gCDSStart)
                t.cCDSStop = str(cds_stop)
                t.gCDSStop = transcript.CM.x2g(cds_stop, 0)
                t.chromCDSStop = GenRecordInstance.record.toChromPos(t.gCDSStop)
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

    @srpc(Mandatory.String, _returns=Mandatory.String)
    def upLoadGenBankLocalFile(content) :
        """
        Not implemented yet.
        """
        raise Fault('ENOTIMPLEMENTED', 'Not implemented yet')
    #upLoadGenBankLocalFile

    @srpc(Mandatory.String, _returns=Mandatory.String)
    def upLoadGenBankRemoteFile(url) :
        """
        Not implemented yet.
        """
        raise Fault('ENOTIMPLEMENTED', 'Not implemented yet')
    #upLoadGenBankRemoteFile

    @srpc(Mandatory.String, Mandatory.String, Mandatory.Integer,
        Mandatory.Integer, _returns=Mandatory.String)
    def sliceChromosomeByGene(geneSymbol, organism, upStream,
        downStream) :
        """
        Todo: documentation, error handling, argument checking, tests.
        """
        O = Output(__file__)
        D = Db.Cache()
        retriever = Retriever.GenBankRetriever(O, D)

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

    @srpc(Mandatory.String, Mandatory.Integer, Mandatory.Integer,
        Mandatory.Integer, _returns=Mandatory.String)
    def sliceChromosome(chromAccNo, start, end, orientation) :
        """
        Todo: documentation, error handling, argument checking, tests.

        @arg orientation: Orientation of the slice. 1 for forward, 2 for
            reverse complement.
        @type orientation: integer
        """
        O = Output(__file__)
        D = Db.Cache()
        retriever = Retriever.GenBankRetriever(O, D)

        O.addMessage(__file__, -1, "INFO",
            "Received request sliceChromosome(%s, %s, %s, %s)" % (
            chromAccNo, start, end, orientation))

        UD = retriever.retrieveslice(chromAccNo, start, end, orientation)

        O.addMessage(__file__, -1, "INFO",
            "Finished processing sliceChromosome(%s, %s, %s, %s)" % (
            chromAccNo, start, end, orientation))

        return UD
    #sliceChromosome

    @srpc(_returns=InfoOutput)
    def info():
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
        output = Output(__file__)
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

    @srpc(_returns=Mandatory.String)
    def ping():
        """
        Simple function to test the interface.

        @return: Always the value 'pong'.
        @rtype: string
        """
        return 'pong'
    #ping

    @srpc(Mandatory.String, Mandatory.String, _returns=Allele)
    def descriptionExtract(reference, observed):
        """
        Extract the HGVS variant description from a reference sequence and an
        observed sequence.

        Note that this only works on DNA sequences for now.
        """
        output = Output(__file__)

        output.addMessage(__file__, -1, 'INFO',
            'Received request descriptionExtract')

        result = Allele()
        result.allele = describe.describe(reference, observed)
        result.description = describe.alleleDescription(result.allele)

        output.addMessage(__file__, -1, 'INFO',
            'Finished processing descriptionExtract')

        return result
    #descriptionExtract

    @srpc(DateTime, _returns=Array(CacheEntry))
    def getCache(created_since=None):
        """
        Get a list of entries from the local cache created since given date.

        This method is intended to be used by Mutalyzer itself to synchronize
        the cache between installations on different servers.
        """
        output = Output(__file__)

        output.addMessage(__file__, -1, 'INFO', 'Received request getCache')

        database = Db.Cache()
        sync = CacheSync(output, database)

        cache = sync.local_cache(created_since)

        def cache_entry_to_soap(entry):
            e = CacheEntry()
            for attr in ('name', 'gi', 'hash', 'chromosomeName',
                'chromosomeStart', 'chromosomeStop', 'chromosomeOrientation',
                'url', 'created', 'cached'):
                setattr(e, attr, entry[attr])
            return e

        output.addMessage(__file__, -1, 'INFO', 'Finished processing getCache')

        return map(cache_entry_to_soap, cache)
    #getCache

    @srpc(Mandatory.String, _returns=Array(Mandatory.String))
    def getdbSNPDescriptions(rs_id):
        """
        Lookup HGVS descriptions for a dbSNP rs identifier.

        @arg rs_id: The dbSNP rs identifier, e.g. 'rs9919552'.
        @type rs_id: string

        @return: List of HGVS descriptions.
        @rtype: list(string)
        """
        output = Output(__file__)

        output.addMessage(__file__, -1, 'INFO',
            'Received request getdbSNPDescription(%s)' % rs_id)

        counter = Db.Counter()
        counter.increment('snpconvert', 'webservice')

        retriever = Retriever.Retriever(output, None)
        descriptions = retriever.snpConvert(rs_id)

        output.addMessage(__file__, -1, 'INFO',
            'Finished processing getdbSNPDescription(%s)' % rs_id)

        # Todo: use SOAP Fault object here (see Trac issue #41).
        messages = output.getMessages()
        if messages:
            error = 'The request could not be completed\n' + \
                '\n'.join(map(lambda m: str(m), output.getMessages()))
            raise Exception(error)

        return descriptions
    #getdbSNPDescriptions
#MutalyzerService


# Close database session at end of each call.
def _shutdown_session(ctx):
    session.remove()
MutalyzerService.event_manager.add_listener('method_return_object',
                                            _shutdown_session)
MutalyzerService.event_manager.add_listener('method_exception_object',
                                            _shutdown_session)
