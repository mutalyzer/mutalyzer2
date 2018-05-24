"""
Mutalyzer RPC services.

@todo: More thourough input checking. The @soap decorator does not do any
       kind of strictness checks on the input. For example, in
       transcriptInfo, the build argument must really be present.
       We should use the built-in validator functionality from Spyne for
       this.
"""


from __future__ import unicode_literals

import binning
from datetime import datetime
from spyne.decorator import rpc, srpc
from spyne.service import ServiceBase
from spyne.model.primitive import Integer, Boolean, DateTime, Unicode
from spyne.model.complex import Array
from spyne.model.fault import Fault
import io
import os
import socket
from operator import attrgetter
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.sql import func

import extractor

import mutalyzer
from mutalyzer.config import settings
from mutalyzer.db import session
from mutalyzer.db import session as sessiongb
from mutalyzer.db.models import (Assembly, Chromosome, BatchJob,
                                 BatchQueueItem, TranscriptMapping)
from mutalyzer.output import Output
from mutalyzer.grammar import Grammar
from mutalyzer.sync import CacheSync
from mutalyzer import announce
from mutalyzer import ncbi
from mutalyzer import stats
from mutalyzer import variantchecker
from mutalyzer.mapping import Converter
from mutalyzer import File
from mutalyzer import Retriever
from mutalyzer import GenRecord
from mutalyzer import Scheduler
from mutalyzer.models import *
from mutalyzer.nc_db import get_entire_nc_record


def create_rpc_fault(output):
    """
    Create an RPC Fault exception from the error message in `output` with the
    highest error level.

    If there are no error messages, construct a generic Fault.
    """
    try:
        message = sorted(output.getMessages(),
                         key=attrgetter('level'),
                         reverse=True)[0]
        return Fault(unicode(message.code), unicode(message.description))
    except IndexError:
        return Fault('ESERVER', 'The request could not be completed')


class MutalyzerService(ServiceBase):
    """
    Mutalyzer web services.

    These methods are made public via a SOAP interface.
    """
    def __init__(self, environ=None):
        super(MutalyzerService, self).__init__(environ)
    #__init__

    @rpc(Mandatory.ByteArray, Unicode, Unicode, Unicode, _returns=Unicode)
    def submitBatchJob(ctx, data, process='NameChecker', argument='', email=None):
        """
        Submit a batch job.

        Input and output file formats for batch jobs are explained on the
        website <https://mutalyzer.nl/batch-jobs>.

        Batch jobs are processed using round-robin scheduling grouped by email
        address (or client IP address if no email address is specified). Per
        email address, jobs are processed sequentially in order of submission.
        This means you will not see any progress on this job until all your
        earlier jobs have finished.

        On error an exception is raised:
          - detail: Human readable description of the error.
          - faultstring: A code to indicate the type of error.
              - EPARSE: The batch input could not be parsed.
              - EMAXSIZE: Input file exceeds maximum size.

        @arg data: Input file (base64 encoded).
        @arg process: Optional type of the batch job, choose from: NameChecker
            (default), SyntaxChecker, PositionConverter, SnpConverter.
        @arg argument: Additional argument. Currently only used if batch_type
            is PositionConverter, denoting the human genome build.
        @arg email: Optional email address. Notification of job completion
            will be sent to this address.

        @return: Batch job identifier.
        """
        output = Output(__file__)

        stats.increment_counter('batch-job/webservice')

        scheduler = Scheduler.Scheduler()
        file_instance = File.File(output)

        batch_types = {'NameChecker': 'name-checker',
                       'SyntaxChecker': 'syntax-checker',
                       'PositionConverter': 'position-converter',
                       'SnpConverter': 'snp-converter'}

        if process not in batch_types:
            raise Fault('EARG',
                        'The process argument must be one of %s.'
                        % ', '.join(batch_types))

        # The Python type for `data` should be a sequence of `str` objects,
        # but it seems we sometimes just get one `str` object. Perhaps only in
        # the unit tests, but let's fix that anyway.
        if isinstance(data, str):
            data = [data]

        # Note that the max file size check below might be bogus, since Spyne
        # first checks the total request size, which by default has a maximum
        # of 2 megabytes.
        # In that case, a senv:Client.RequestTooLong faultstring is returned.

        # Todo: Set maximum request size by specifying the max_content_length
        #     argument for spyne.server.wsgi.WsgiApplication in all webservice
        #     instantiations.
        if sum(len(s) for s in data) > settings.MAX_FILE_SIZE:
            raise Fault('EMAXSIZE',
                        'Only files up to %d megabytes are accepted.'
                        % (settings.MAX_FILE_SIZE // 1048576))

        batch_file = io.BytesIO()
        for d in data:
            batch_file.write(d)

        job, columns = file_instance.parseBatchFile(batch_file)
        batch_file.close()

        if job is None:
            raise Fault('EPARSE', 'Could not parse input file, please check your file format.')

        if not email:
            # If no email address is specified, we create a fake one based on
            # the caller's IP address. This makes sure the scheduler processes
            # jobs grouped by user.
            try:
                address = unicode(ctx.transport.req_env['REMOTE_ADDR'])
            except (AttributeError, KeyError):
                address = 'localhost'
            email = '%s@webservice.mutalyzer' % address

        result_id = scheduler.addJob(email, job, columns,
                                     batch_types[process], argument)
        return result_id

    @srpc(Mandatory.Unicode, _returns=Integer)
    def monitorBatchJob(job_id):
        """
        Get the number of entries left for a batch job.

        Input and output file formats for batch jobs are explained on the
        website <https://mutalyzer.nl/batch-jobs>.

        @arg job_id: Batch job identifier.

        @return: Number of entries left.
        """
        return BatchQueueItem.query.join(BatchJob).filter_by(result_id=job_id).count()

    @srpc(Mandatory.Unicode, _returns=ByteArray)
    def getBatchJob(job_id):
        """
        Get the result of a batch job.

        Input and output file formats for batch jobs are explained on the
        website <https://mutalyzer.nl/batch-jobs>.

        On error an exception is raised:
          - detail: Human readable description of the error.
          - faultstring: A code to indicate the type of error.
              - EBATCHNOTREADY: The batch job result is not yet ready.

        @arg job_id: Batch job identifier.

        @return: Batch job result file (UTF-8, base64 encoded).
        """
        left = BatchQueueItem.query.join(BatchJob).filter_by(result_id=job_id).count()

        if left > 0:
            raise Fault('EBATCHNOTREADY', 'Batch job result is not yet ready.')

        filename = 'batch-job-%s.txt' % job_id
        handle = open(os.path.join(settings.CACHE_DIR, filename), 'rb')
        return handle

    @srpc(Mandatory.Unicode, Mandatory.Unicode, Mandatory.Integer, Boolean,
        _returns=Array(Mandatory.Unicode))
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
        @arg pos: A position on the chromosome (one-based).
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

        try:
            assembly = Assembly.by_name_or_alias(build)
        except NoResultFound:
            L.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                            "build name." % build)

        try:
            chromosome = assembly.chromosomes.filter_by(name=chrom).one()
        except NoResultFound:
            L.addMessage(__file__, 4, "EARG", "EARG %s" % chrom)
            raise Fault("EARG", "The chrom argument (%s) was not a valid " \
                            "chromosome name." % chrom)

        pos = max(min(pos, binning.MAX_POSITION + 1), 1)
        bins = binning.overlapping_bins(pos - 1, pos)
        mappings = chromosome.transcript_mappings.filter(
            TranscriptMapping.bin.in_(bins),
            TranscriptMapping.start <= pos,
            TranscriptMapping.stop >= pos
        ).order_by(
            TranscriptMapping.start,
            TranscriptMapping.stop,
            TranscriptMapping.gene,
            TranscriptMapping.accession,
            TranscriptMapping.version,
            TranscriptMapping.transcript)

        L.addMessage(__file__, -1, "INFO",
                     "Finished processing getTranscripts(%s %s %s %s)"
                     % (build, chrom, pos, versions))

        return [mapping.get_reference(include_version=versions)
                for mapping in mappings]
    #getTranscripts

    @srpc(Mandatory.Unicode, Mandatory.Unicode, _returns=Array(Mandatory.Unicode))
    def getTranscriptsByGeneName(build, name):
        """
        Todo: documentation.
        """
        L = Output(__file__)

        L.addMessage(__file__, -1, "INFO",
            "Received request getTranscriptsByGene(%s %s)" % (build, name))

        try:
            assembly = Assembly.by_name_or_alias(build)
        except NoResultFound:
            L.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                            "build name." % build)

        mappings = TranscriptMapping.query \
            .filter(TranscriptMapping.chromosome.has(assembly=assembly),
                    TranscriptMapping.gene == name)

        L.addMessage(__file__, -1, "INFO",
            "Finished processing getTranscriptsByGene(%s %s)" % (build, name))

        return [mapping.reference for mapping in mappings]
    #getTranscriptsByGeneName

    @srpc(Mandatory.Unicode, Mandatory.Unicode, Mandatory.Integer,
          Mandatory.Integer, Mandatory.Integer, Boolean,
          _returns=Array(Mandatory.Unicode))
    def getTranscriptsRange(build, chrom, pos1, pos2, method, versions=False):
        """
        Get all the transcripts that overlap with a range on a chromosome.

        The range should be provided as one-based, inclusive positions.

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
        @kwarg versions: If set to True, also include transcript versions.
        @type versions: bool

        @return: A list of transcripts.
        @rtype: list
        """
        L = Output(__file__)

        L.addMessage(__file__, -1, "INFO",
            "Received request getTranscriptsRange(%s %s %s %s %s)" % (build,
            chrom, pos1, pos2, method))

        if pos1 > pos2:
            L.addMessage(__file__, 4, 'EARG',
                         'Invalid range: %d-%d' % (pos1, pos2))
            raise Fault('EARG',
                        'Invalid range (%d-%d) with start position greater '
                        'than stop position.' % (pos1, pos2))

        if pos1 < 1 or pos2 > binning.MAX_POSITION + 1:
            L.addMessage(__file__, 4, 'EARG',
                         'Invalid range: %d-%d' % (pos1, pos2))
            raise Fault('EARG',
                        'Invalid range (%d-%d) exceeding chromosome bounds.'
                        % (pos1, pos2))

        try:
            assembly = Assembly.by_name_or_alias(build)
        except NoResultFound:
            L.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                            "build name." % build)

        try:
            chromosome = assembly.chromosomes.filter_by(name=chrom).one()
        except NoResultFound:
            L.addMessage(__file__, 4, "EARG", "EARG %s" % chrom)
            raise Fault("EARG", "The chrom argument (%s) was not a valid " \
                            "chromosome name." % chrom)

        if method:
            bins = binning.overlapping_bins(pos1 - 1, pos2)
            range_filter = (TranscriptMapping.bin.in_(bins),
                            TranscriptMapping.start <= pos2,
                            TranscriptMapping.stop >= pos1)
        else:
            bins = binning.contained_bins(pos1 - 1, pos2)
            range_filter = (TranscriptMapping.bin.in_(bins),
                            TranscriptMapping.start >= pos1,
                            TranscriptMapping.stop <= pos2)

        mappings = chromosome.transcript_mappings.filter(*range_filter).order_by(
            TranscriptMapping.start,
            TranscriptMapping.stop,
            TranscriptMapping.gene,
            TranscriptMapping.accession,
            TranscriptMapping.version,
            TranscriptMapping.transcript)

        L.addMessage(__file__, -1, "INFO",
            "Finished processing getTranscriptsRange(%s %s %s %s %s)" % (
            build, chrom, pos1, pos2, method))

        return [mapping.get_reference(include_version=versions)
                for mapping in mappings]
    #getTranscriptsRange

    @srpc(Mandatory.Unicode, Mandatory.Unicode, Mandatory.Integer,
        Mandatory.Integer, Mandatory.Integer,
        _returns=Array(TranscriptMappingInfo))
    def getTranscriptsMapping(build, chrom, pos1, pos2, method):
        """
        Get all the transcripts and their info that overlap with a range on a
        chromosome.

        The range should be provided as one-based, inclusive positions.

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
                 - transcript
                 - name
                 - version
                 - gene
                 - orientation
                 - start
                 - stop
                 - cds_start
                 - cds_stop
        All returned ranges are one-based, inclusive, and in gene
        orientation.
        """
        output = Output(__file__)
        output.addMessage(__file__, -1, 'INFO', 'Received request ' \
            'getTranscriptsMapping(%s %s %s %s %s)' % (build, chrom, pos1, pos2,
            method))

        if pos1 > pos2:
            output.addMessage(__file__, 4, 'EARG',
                         'Invalid range: %d-%d' % (pos1, pos2))
            raise Fault('EARG',
                        'Invalid range (%d-%d) with start position greater '
                        'than stop position.' % (pos1, pos2))

        if pos1 < 1 or pos2 > binning.MAX_POSITION + 1:
            L.addMessage(__file__, 4, 'EARG',
                         'Invalid range: %d-%d' % (pos1, pos2))
            raise Fault('EARG',
                        'Invalid range (%d-%d) exceeding chromosome bounds.'
                        % (pos1, pos2))

        try:
            assembly = Assembly.by_name_or_alias(build)
        except NoResultFound:
            output.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                            "build name." % build)

        try:
            chromosome = assembly.chromosomes.filter_by(name=chrom).one()
        except NoResultFound:
            output.addMessage(__file__, 4, "EARG", "EARG %s" % chrom)
            raise Fault("EARG", "The chrom argument (%s) was not a valid " \
                            "chromosome name." % chrom)

        if method:
            bins = binning.overlapping_bins(pos1 - 1, pos2)
            range_filter = (TranscriptMapping.bin.in_(bins),
                            TranscriptMapping.start <= pos2,
                            TranscriptMapping.stop >= pos1)
        else:
            bins = binning.contained_bins(pos1 - 1, pos2)
            range_filter = (TranscriptMapping.bin.in_(bins),
                            TranscriptMapping.start >= pos1,
                            TranscriptMapping.stop <= pos2)

        mappings = chromosome.transcript_mappings.filter(*range_filter).order_by(
            TranscriptMapping.start,
            TranscriptMapping.stop,
            TranscriptMapping.gene,
            TranscriptMapping.accession,
            TranscriptMapping.version,
            TranscriptMapping.transcript)

        transcripts = []

        for mapping in mappings:
            t = TranscriptMappingInfo()
            t.transcript = mapping.reference
            t.name = mapping.accession
            t.version = mapping.version
            t.gene = mapping.gene
            t.orientation = '-' if mapping.orientation == 'reverse' else '+'
            if mapping.orientation == 'reverse':
                t.start, t.stop = mapping.stop, mapping.start
            else:
                t.start, t.stop = mapping.start, mapping.stop
            if mapping.orientation == 'reverse':
                t.cds_start, t.cds_stop = mapping.cds_stop, mapping.cds_start
            else:
                t.cds_start, t.cds_stop = mapping.cds_start, mapping.cds_stop
            transcripts.append(t)

        output.addMessage(__file__, -1, 'INFO', 'Finished processing ' \
            'getTranscriptsMapping(%s %s %s %s %s)' % (build, chrom, pos1, pos2,
            method))

        return transcripts
    #getTranscriptsMapping

    @srpc(Mandatory.Unicode, Mandatory.Unicode, _returns=Mandatory.Unicode)
    def getGeneName(build, accno) :
        """
        Find the gene name associated with a transcript.

        @arg build: The genome build (hg19, hg18, mm10).
        @type build: string
        @arg accno: The identifier of a transcript, with optional version.
        @type accno: string

        @return: The name of the associated gene.
        @rtype: string
        """
        L = Output(__file__)

        L.addMessage(__file__, -1, "INFO",
            "Received request getGeneName(%s %s)" % (build, accno))

        try:
            assembly = Assembly.by_name_or_alias(build)
        except NoResultFound:
            L.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                            "build name." % build)

        # Accept an optional accession version.
        try:
            accno, version = accno.split('.')
            version = int(version)
        except ValueError:
            version = None

        mapping = TranscriptMapping.query.filter(
            TranscriptMapping.chromosome.has(assembly=assembly),
            TranscriptMapping.accession == accno
        )
        if version:
            mapping = mapping.filter_by(version=version)
        mapping = mapping.order_by(TranscriptMapping.version.desc()).first()

        L.addMessage(__file__, -1, "INFO",
            "Finished processing getGeneName(%s %s)" % (build, accno))

        if not mapping:
            L.addMessage(__file__, 4, "ENOTFOUND", "ENOTFOUND %s %s %s" % (build, accno, version))
            raise Fault("ENOTFOUND",
                        "Transcript %s%s not found for build %s." % (
                            accno, '.%d' % version if version else '', build))
        return mapping.gene
    #getGeneName

    @srpc(Mandatory.Unicode, Mandatory.Unicode, Mandatory.Unicode,
        Mandatory.Unicode, _returns=Mapping)
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
        @arg accNo: The NM accession number and version or LRG.
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

        try:
            assembly = Assembly.by_name_or_alias(build)
        except NoResultFound:
            L.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                            "build name." % build)

        conv = Converter(assembly, L)
        result = conv.mainMapping(accNo, variant)

        L.addMessage(__file__, -1, "INFO",
            "Finished processing mappingInfo(%s %s %s %s)" % (LOVD_ver, build,
            accNo, variant))

        del L
        return result
    #mappingInfo

    @srpc(Mandatory.Unicode, Mandatory.Unicode, Mandatory.Unicode,
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
        @arg accNo: The NM accession number and version or LRG.
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

        try:
            assembly = Assembly.by_name_or_alias(build)
        except NoResultFound:
            O.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                            "build name." % build)

        converter = Converter(assembly, O)
        T = converter.mainTranscript(accNo)

        O.addMessage(__file__, -1, "INFO",
            "Finished processing transcriptInfo(%s %s %s)" % (LOVD_ver, build,
            accNo))
        return T
    #transcriptInfo

    @srpc(Mandatory.Unicode, Mandatory.Unicode, _returns=Mandatory.Unicode)
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
        L = Output(__file__)
        L.addMessage(__file__, -1, "INFO",
            "Received request chromAccession(%s %s)" % (build, name))

        try:
            assembly = Assembly.by_name_or_alias(build)
        except NoResultFound:
            L.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                            "build name." % build)

        try:
            chromosome = assembly.chromosomes.filter_by(name=name).one()
        except NoResultFound:
            L.addMessage(__file__, 4, "EARG", "EARG %s" % name)
            raise Fault("EARG", "The name argument (%s) was not a valid " \
                            "chromosome name." % name)

        L.addMessage(__file__, -1, "INFO",
            "Finished processing chromAccession(%s %s)" % (build, name))

        return chromosome.accession
    #chromAccession

    @srpc(Mandatory.Unicode, Mandatory.Unicode, _returns=Mandatory.Unicode)
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
        L = Output(__file__)
        L.addMessage(__file__, -1, "INFO",
            "Received request chromName(%s %s)" % (build, accNo))

        try:
            assembly = Assembly.by_name_or_alias(build)
        except NoResultFound:
            L.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                            "build name." % build)

        try:
            chromosome = assembly.chromosomes.filter_by(accession=accNo).one()
        except NoResultFound:
            L.addMessage(__file__, 4, "EARG", "EARG %s" % accNo)
            raise Fault("EARG", "The accNo argument (%s) was not a valid " \
                            "chromosome accession." % accNo)

        L.addMessage(__file__, -1, "INFO",
            "Finished processing chromName(%s %s)" % (build, accNo))

        return chromosome.name
    #chromosomeName

    @srpc(Mandatory.Unicode, Mandatory.Unicode, _returns=Mandatory.Unicode)
    def getchromName(build, acc) :
        """
        Get the chromosome name, given a transcript identifier (NM number).

        @arg build: The genome build (hg19, hg18, mm10).
        @type build: string
        @arg acc: The NM accession number (version NOT included) or LRG.
        @type acc: string

        @return: The name of a chromosome.
        @rtype: string
        """
        L = Output(__file__)

        L.addMessage(__file__, -1, "INFO",
            "Received request getchromName(%s %s)" % (build, acc))

        try:
            assembly = Assembly.by_name_or_alias(build)
        except NoResultFound:
            L.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                            "build name." % build)

        mapping = TranscriptMapping.query \
            .filter(TranscriptMapping.chromosome.has(assembly=assembly),
                    TranscriptMapping.accession == acc) \
            .first()

        L.addMessage(__file__, -1, "INFO",
            "Finished processing getchromName(%s %s)" % (build, acc))

        return mapping.chromosome.name
    #chromosomeName

    @srpc(Mandatory.Unicode, Mandatory.Unicode, Unicode,
          _returns=Array(Mandatory.Unicode))
    def numberConversion(build, variant, gene=None):
        """
        Converts I{c.} to I{g.} notation or vice versa

        @arg build: The genome build (hg19, hg18, mm10).
        @type build: string
        @arg variant: The variant in either I{c.} or I{g.} notation, full HGVS
            notation, including NM_, NC_, or LRG_ accession number.
        @type variant: string
        @kwarg gene: Optional gene name. If given, return variant descriptions
            on all transcripts for this gene.
        @type gene: string

        @return: The variant(s) in either I{g.} or I{c.} notation.
        @rtype: list
        """
        O = Output(__file__)
        O.addMessage(__file__, -1, "INFO",
            "Received request cTogConversion(%s %s)" % (build, variant))

        stats.increment_counter('position-converter/webservice')

        try:
            assembly = Assembly.by_name_or_alias(build)
        except NoResultFound:
            O.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                            "build name." % build)

        converter = Converter(assembly, O)
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

    @srpc(Mandatory.Unicode, _returns=CheckSyntaxOutput)
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

        stats.increment_counter('syntax-checker/webservice')

        if not variant :
            output.addMessage(__file__, 4, "EARG", "EARG no variant")
            raise Fault("EARG", "The variant argument is not provided.")

        result = CheckSyntaxOutput()

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

    @srpc(Mandatory.Unicode, _returns=MutalyzerOutput)
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
            - legend: Array with name and information per transcript variant
                and protein isoform. Each entry is an array with the following
                fields:
                - name
                - id
                - locusTag
                - product
                - linkMethod
            - messages: List of (error) messages.
        """
        O = Output(__file__)
        O.addMessage(__file__, -1, "INFO",
            "Received request runMutalyzer(%s)" % (variant))

        stats.increment_counter('name-checker/webservice')

        variantchecker.check_variant(variant, O)

        result = MutalyzerOutput()

        result.referenceId = O.getIndexedOutput('reference_id', 0)
        result.sourceId = O.getIndexedOutput('source_id', 0)
        result.sourceAccession = O.getIndexedOutput('source_accession', 0)
        result.sourceVersion = O.getIndexedOutput('source_version', 0)
        result.molecule = O.getIndexedOutput('molecule', 0)

        result.original = O.getIndexedOutput("original", 0)
        result.mutated = O.getIndexedOutput("mutated", 0)

        result.origMRNA = O.getIndexedOutput("origMRNA", 0)
        result.mutatedMRNA = O.getIndexedOutput("mutatedMRNA", 0)

        result.origCDS = O.getIndexedOutput("origCDS", 0)
        result.newCDS = O.getIndexedOutput("newCDS", 0)

        result.origProtein = O.getIndexedOutput("oldProtein", 0)
        result.newProtein = O.getIndexedOutput("newProtein", 0)
        result.altProtein = O.getIndexedOutput("altProtein", 0)

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

        result.legend = []
        for name, id, locusTag, product, linkMethod in O.getOutput('legends'):
            record = LegendRecord()
            record.name = name
            record.id = id
            record.locusTag = locusTag
            record.product = product
            record.linkMethod = linkMethod
            result.legend.append(record)

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

    @srpc(Mandatory.Unicode, RequestExtras, _returns=MutalyzerOutput)
    def runMutalyzerLight(variant, extras) :
        """
        Run the Mutalyzer name checker.

        @arg variant: The variant description to check.
        @type variant: string

        @arg extras: Additional fields to be included in the response.
        @type extras: RequestExtras


        @return: Default response object contains the following fields:
            - referenceId: Identifier of the reference sequence used.
            - sourceId: Identifier of the reference sequence source, e.g. the
                chromosomal accession number and version in case referenceId
                is a  UD reference created as a chromosomal slice.
            - sourceAccession: Accession number of the reference sequence
                source (only for genbank references).
            - sourceVersion: Version number of the reference sequence source
                (only for genbank references).
            - molecule: Molecular type of the reference sequence.
            - errors: Number of errors.
            - warnings: Number of warnings.
            - summary: Summary of messages.
            - chromDescription: Chromosomal description.
            - genomicDescription: Genomic description.
            - transcriptDescriptions: List of transcript descriptions.
            - proteinDescriptions: List of protein descriptions.
            - exons: If a transcript is selected, array of ExonInfo objects
                for each exon in the selected transcript with fields:
                - cStart
                - gStart
                - cStop
                - gStop
            - legend: Array with name and information per transcript variant
                and protein isoform. Each entry is an array with the following
                fields:
                - name
                - id
                - locusTag
                - product
                - linkMethod
            - messages: List of (error) messages.
            If extras is utilized the following fields are included, according
            to their selection:
                - original: Original sequence.
                - mutated: Mutated sequence.
                - varDetails:
                   - reference_file: reference file
                   - ref: reference bases ('.' for insertions)
                   - alt: alternative bases ('.' for deletions)
                   - start: start position
                   - stop: stop position
                   - info: information about the conversion process failure
        """
        def check_param(obj, attr=None):
            """
            Checks whether `attr` is present in the `obj`. Added specifically
            for the `extras` parameter.
            :param obj: Call parameter.
            :param attr: Attribute to search for in obj.
            :return:
            """
            if not obj:
                return False
            if obj and not attr:
                return True
            if hasattr(obj, attr) and getattr(obj, attr):
                return True
            else:
                return False

        def get_variant_details():
            """
            Extracts the variant description details from the parse tree and
            the sequence.

            Currently accepts only one variant description and three mutation
            types ('del', 'ins', and 'subst') with their corresponding
            locations as natural numbers only ('c.10-3_10-5del' and other
            positions such as '*3' are not accepted). Thus, it is mainly for
            genomic references with 'g.' coordinates and no specific annotated
            segment.

            Example of accepted variant descriptions and their output terms:
            - input: 'NC_000014.8:g.94844865G>A'
                - output terms:
                    - 'reference_file': 'NC_000014.8'
                    - 'start': '94844865'
                    - 'stop': '94844865'
                    - 'ref': 'G'
                    - 'alt': 'A'
                    - 'operation': 'subst'
            - input: 'NC_000011.9:g.47353833_47353857del'
                - output terms:
                    - 'reference_file': 'NC_000011.9'
                    - 'start': '47353833'
                    - 'stop': '47353857'
                    - 'ref': 'AGGGAAGCCATCCAGGCTGAGAGGG'
                    - 'alt': '.'
                    - 'operation': 'del'
            - input: 'NC_000008.10:g.10480387_10480388dup'
                - output terms:
                    - 'reference_file': 'NC_000008.10'
                    - 'start': '10480387'
                    - 'stop': '10480388'
                    - 'ref': '.'
                    - 'alt': 'GC'
                    - 'operation': 'dup'
            - input: 'NC_000008.10:g.10480387_10480388insA'
                - output terms:
                    - 'reference_file': 'NC_000008.10'
                    - 'start': '10480388'
                    - 'stop': '10480388'
                    - 'ref': '.'
                    - 'alt': 'A'
                    - 'operation': 'ins'
            - input: 'NC_000014.8:g.19400000delinsGT'
                - output terms:
                    - 'reference_file': 'NC_000014.8'
                    - 'stop': '19400000'
                    - 'start': '19400000'
                    - 'ref': 'A'
                    - 'alt': 'GT'
                    - 'operation': 'delins'

            For 'del' variants the sequence is required.
            For 'ins' and 'subst' variants the sequence is not required with
            the 'ref' and 'alt' terms being obtained directly from the HGVS
            description.
            """
            varDetails = VarDetails()
            varDetails.info = []
            output = Output(__file__)
            grammar = Grammar(output)
            if result.genomicDescription:
                variant_tree = grammar.parse(result.genomicDescription)
            else:
                varDetails.info.append("No genomic description generated.")
                return varDetails

            # Only one variant accepted for the moment
            if variant_tree and not variant_tree.SingleAlleleVarSet:
                # Extract the only variant
                variant = variant_tree.RawVar

                # Should be NCBI and not LRG
                if variant_tree.RefSeqAcc and variant_tree.Version:
                    record_id = variant_tree.RefSeqAcc + '.' + variant_tree.Version
                    varDetails.reference_file = record_id
                else:
                    varDetails.info.append("Only ncbi references with accession "
                                           "and version are accepted.")
                    return varDetails

                # Extract the variant type
                try:
                    mutation_type = variant.MutationType
                except AttributeError:
                    varDetails.info.append("No operation (mutation) type found.")
                    return varDetails

                # Extract the the positions
                varDetails.start = varDetails.stop = variant.StartLoc.PtLoc.Main
                if variant.EndLoc:
                    varDetails.stop = variant.EndLoc.PtLoc.Main
                if abs(int(varDetails.stop) - int(varDetails.start)) > 524288:
                    varDetails.info.append("Too long (> 524288 bases) sequence change.")
                    return varDetails

                # Extract the 'ref' and 'alt' terms
                if mutation_type == 'subst':
                    varDetails.ref = unicode(variant.Arg1)
                    varDetails.alt = unicode(variant.Arg2)
                    varDetails.operation = 'subst'
                elif mutation_type == 'del':
                    varDetails.ref = unicode(O.getIndexedOutput("original", 0)\
                        [int(varDetails.start)-1:int(varDetails.stop)])
                    varDetails.alt = '.'
                    varDetails.operation = 'del'
                elif mutation_type == 'dup':
                    varDetails.ref = '.'
                    varDetails.alt = unicode(O.getIndexedOutput("original", 0)\
                        [int(varDetails.start)-1:int(varDetails.stop)])
                    varDetails.operation = 'dup'
                elif mutation_type == 'ins':
                    varDetails.ref = '.'
                    varDetails.alt = unicode(variant.Seq.Sequence)
                    varDetails.operation = 'ins'
                elif mutation_type == 'delins':
                    varDetails.ref = unicode(O.getIndexedOutput("original", 0)\
                        [int(varDetails.start)-1:int(varDetails.stop)])
                    varDetails.alt = unicode(variant.Seq.Sequence)
                    varDetails.operation = 'delins'
                else:
                    varDetails.info.append("Conversion not performed since "
                                           "'%s' operation is not supported."
                                           % mutation_type)
                    return varDetails
            else:
                varDetails.info.append("Conversion not performed since multiple "
                                       "variants are present in the description.")
                return varDetails

            return varDetails

        O = Output(__file__)
        O.addMessage(__file__, -1, "INFO",
            "Received request runMutalyzerLight(%s)" % (variant))

        stats.increment_counter('name-checker/webservice')

        if not check_param(extras, 'original') and not check_param(extras,'varDetails'):
            O.addOutput('add_original_sequence_to_output', 'not')
        if not check_param(extras, 'mutated'):
            O.addOutput('add_mutated_sequence_to_output', 'not')
        variantchecker.check_variant(variant, O)

        result = MutalyzerOutput()

        result.referenceId = O.getIndexedOutput('reference_id', 0)
        result.sourceId = O.getIndexedOutput('source_id', 0)
        result.sourceAccession = O.getIndexedOutput('source_accession', 0)
        result.sourceVersion = O.getIndexedOutput('source_version', 0)
        result.molecule = O.getIndexedOutput('molecule', 0)

        if check_param(extras, 'original'):
            result.original = O.getIndexedOutput("original", 0)
        if check_param(extras, 'mutated'):
            result.mutated = O.getIndexedOutput("mutated", 0)

        result.chromDescription = \
            O.getIndexedOutput("genomicChromDescription", 0)
        result.genomicDescription = \
            O.getIndexedOutput("genomicDescription", 0)

        if check_param(extras, 'varDetails'):
            result.varDetails = get_variant_details()

        result.transcriptDescriptions = O.getOutput("descriptions")
        result.proteinDescriptions = O.getOutput("protDescriptions")

        if O.getIndexedOutput('hasTranscriptInfo', 0, False):
            result.exons = []
            for e in O.getOutput('exonInfo'):
                exon = ExonInfo()
                exon.gStart, exon.gStop, exon.cStart, exon.cStop = e
                result.exons.append(exon)

        result.legend = []
        for name, id, locusTag, product, linkMethod in O.getOutput('legends'):
            record = LegendRecord()
            record.name = name
            record.id = id
            record.locusTag = locusTag
            record.product = product
            record.linkMethod = linkMethod
            result.legend.append(record)

        result.errors, result.warnings, result.summary = O.Summary()

        O.addMessage(__file__, -1, "INFO",
            "Finished processing runMutalyzerLight(%s)" % (variant))

        result.messages = []
        for message in O.getMessages():
            soap_message = SoapMessage()
            soap_message.errorcode = message.code
            soap_message.message = message.description
            result.messages.append(soap_message)

        return result
    #runMutalyzerLight

    @srpc(Mandatory.Unicode, Mandatory.Unicode, _returns=TranscriptNameInfo)
    def getGeneAndTranscript(genomicReference, transcriptReference) :
        """
        Todo: documentation.
        """
        O = Output(__file__)

        O.addMessage(__file__, -1, "INFO",
            "Received request getGeneAndTranscript(%s, %s)" % (
            genomicReference, transcriptReference))
        retriever = Retriever.GenBankRetriever(O)
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

    @srpc(Mandatory.Unicode, Unicode, _returns=Array(TranscriptInfo))
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

        O.addMessage(__file__, -1, "INFO",
            "Received request getTranscriptsAndInfo(%s, %s)" % (
            genomicReference, geneName))

        # We try first to
        if 'NC' in genomicReference:
            record = get_entire_nc_record(genomicReference, geneName=geneName)
        else:
            record = None

        if record is None:
            retriever = Retriever.GenBankRetriever(O)
            record = retriever.loadrecord(genomicReference)

        if record is None:
            raise Fault("EARG",
                        "Unknown reference file: %s" % genomicReference)

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

                t.cTransEnd = unicode(t.exons[-1].cStop)
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
                t.cTransStart = unicode(trans_start)
                t.gTransStart = transcript.CM.x2g(trans_start, 0)
                t.chromTransStart = GenRecordInstance.record.toChromPos(
                    t.gTransStart)
                t.cCDSStart = unicode(cds_start)
                t.gCDSStart = transcript.CM.x2g(cds_start, 0)
                t.chromCDSStart = GenRecordInstance.record.toChromPos(
                    t.gCDSStart)
                t.cCDSStop = unicode(cds_stop)
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

    @srpc(Mandatory.ByteArray, _returns=Mandatory.Unicode)
    def uploadGenBankLocalFile(data):
        """
        Upload a genbank file.

        @arg data: Genbank file (UTF-8, base64 encoded).
        @return: UD accession number for the uploaded genbank file.
        """
        output = Output(__file__)
        retriever = Retriever.GenBankRetriever(output)

        output.addMessage(__file__, -1, 'INFO',
                          'Received request uploadGenBankLocalFile()')

        # The Python type for `data` should be a sequence of `str` objects,
        # but it seems we sometimes just get one `str` object. Perhaps only in
        # the unit tests, but let's fix that anyway.
        if isinstance(data, str):
            data = [data]

        # Note that the max file size check below might be bogus, since Spyne
        # first checks the total request size, which by default has a maximum
        # of 2 megabytes.
        # In that case, a senv:Client.RequestTooLong faultstring is returned.

        # Todo: Set maximum request size by specifying the max_content_length
        #     argument for spyne.server.wsgi.WsgiApplication in all webservice
        #     instantiations.
        if sum(len(s) for s in data) > settings.MAX_FILE_SIZE:
            raise Fault('EMAXSIZE',
                        'Only files up to %d megabytes are accepted.'
                        % (settings.MAX_FILE_SIZE // 1048576))

        ud = retriever.uploadrecord(b''.join(data))

        output.addMessage(__file__, -1, 'INFO',
                          'Finished processing uploadGenBankLocalFile()')

        if not ud:
            raise create_rpc_fault(output)

        return ud
    #upLoadGenBankLocalFile

    @srpc(Mandatory.Unicode, _returns=Mandatory.Unicode)
    def uploadGenBankRemoteFile(url) :
        """
        Not implemented yet.
        """
        raise Fault('ENOTIMPLEMENTED', 'Not implemented yet')
    #upLoadGenBankRemoteFile

    @srpc(Mandatory.Unicode, Mandatory.Unicode, Mandatory.Integer,
        Mandatory.Integer, _returns=Mandatory.Unicode)
    def sliceChromosomeByGene(geneSymbol, organism, upStream,
        downStream) :
        """
        Retrieve part of the reference genome for a (HGNC) gene symbol.

        @arg geneSymbol: Gene symbol.
        @type geneSymbol: string
        @arg organism: Organism name without spaces.
        @type organism: string
        @arg upStream: Number of 5' flanking bases to include.
        @type upStream: integer
        @arg downStream: Number of 3' flanking bases to include.
        @type upStream: integer

        This uses the NCBI Entrez search engine and is therefore based on the
        current Entrez assembly for the given organism.

        @return: UD accession number for created slice.
        """
        # Todo: error handling, argument checking, tests.
        O = Output(__file__)
        retriever = Retriever.GenBankRetriever(O)

        O.addMessage(__file__, -1, "INFO",
            "Received request sliceChromosomeByGene(%s, %s, %s, %s)" % (
            geneSymbol, organism, upStream, downStream))

        UD = retriever.retrievegene(geneSymbol, organism, upStream, downStream)

        O.addMessage(__file__, -1, "INFO",
            "Finished processing sliceChromosomeByGene(%s, %s, %s, %s)" % (
            geneSymbol, organism, upStream, downStream))

        if not UD:
            raise create_rpc_fault(O)

        return UD
    #sliceChromosomeByGene

    @srpc(Mandatory.Unicode, Mandatory.Integer, Mandatory.Integer,
        Mandatory.Integer, _returns=Mandatory.Unicode)
    def sliceChromosome(chromAccNo, start, end, orientation) :
        """
        Retrieve a range of a chromosome by accession number.

        @arg chromAccNo: Chromosome or contig by accession number.
        @type chromAccNo: string
        @arg start: Start position (one-based, inclusive, in reference
          orientation).
        @type start: integer
        @arg end: End position (one-based, inclusive, in reference
          orientation).
        @type end: integer
        @arg orientation: Orientation of the slice. 1 for forward, 2 for
          reverse complement.
        @type orientation: integer
        """
        # Todo: error handling, argument checking, tests.
        O = Output(__file__)
        retriever = Retriever.GenBankRetriever(O)

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
            - announcement: Announcement body text.
            - announcementUrl: URL to go with the announcement for further
                information.
        @rtype: object
        """
        output = Output(__file__)
        output.addMessage(__file__, -1, 'INFO', 'Received request info')

        result = InfoOutput()
        result.version = mutalyzer.__version__
        result.versionParts = mutalyzer.__version_info__
        if mutalyzer.__version_info__[-1] == 'dev':
            result.releaseDate = ''
        else:
            result.releaseDate = mutalyzer.__date__
        result.nomenclatureVersion = mutalyzer.NOMENCLATURE_VERSION
        result.nomenclatureVersionParts = mutalyzer.NOMENCLATURE_VERSION_INFO
        result.serverName = socket.gethostname()
        result.contactEmail = mutalyzer.__contact__

        announcement = announce.get_announcement()
        if announcement:
            result.announcement = announcement['body']
            result.announcementUrl = announcement['url']

        output.addMessage(__file__, -1, 'INFO', 'Finished processing info')
        return result
    #info

    @srpc(_returns=Mandatory.Unicode)
    def ping():
        """
        Simple function to test the interface.

        @return: Always the value 'pong'.
        @rtype: string
        """
        return 'pong'
    #ping

    @srpc(Mandatory.Unicode, Mandatory.Unicode, _returns=Allele)
    def descriptionExtract(reference, observed):
        """
        Extract the HGVS variant description from a reference sequence and an
        observed sequence.

        Note that this only works on DNA sequences for now.
        """
        output = Output(__file__)

        output.addMessage(__file__, -1, 'INFO',
            'Received request descriptionExtract')

        stats.increment_counter('description-extractor/webservice')

        if (len(reference) > settings.EXTRACTOR_MAX_INPUT_LENGTH or
            len(observed) > settings.EXTRACTOR_MAX_INPUT_LENGTH):
            raise Fault('EMAXSIZE',
                        'Input sequences are restricted to {} bp.'
                        .format(settings.EXTRACTOR_MAX_INPUT_LENGTH))

        allele = extractor.describe_dna(reference, observed)

        result = Allele()
        result.allele = []
        for variant in allele:
            raw_var = RawVar()
            raw_var.start = variant.start
            raw_var.start_offset = variant.start_offset
            raw_var.end = variant.end
            raw_var.end_offset = variant.end_offset
            raw_var.sample_start = variant.sample_start
            raw_var.sample_start_offset = variant.sample_start_offset
            raw_var.sample_end = variant.sample_end
            raw_var.sample_end_offset = variant.sample_end_offset
            raw_var.type = variant.type
            raw_var.deleted = unicode(variant.deleted)
            raw_var.inserted = unicode(variant.inserted)
            raw_var.weight = variant.weight()
            raw_var.shift = variant.shift
            raw_var.description = unicode(variant)
            result.allele.append(raw_var)
        result.description = unicode(allele)

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

        sync = CacheSync(output)

        cache = sync.local_cache(created_since)

        def cache_entry_to_soap(entry):
            e = CacheEntry()
            for attr in ('name', 'source', 'source_data', 'hash', 'created',
                         'cached'):
                setattr(e, attr, entry[attr])
            return e

        output.addMessage(__file__, -1, 'INFO', 'Finished processing getCache')

        return map(cache_entry_to_soap, cache)
    #getCache

    @srpc(Mandatory.Unicode, _returns=Array(Mandatory.Unicode))
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

        stats.increment_counter('snp-converter/webservice')

        descriptions = ncbi.rsid_to_descriptions(rs_id, output)

        output.addMessage(__file__, -1, 'INFO',
            'Finished processing getdbSNPDescription(%s)' % rs_id)

        if output.getMessages():
            raise create_rpc_fault(output)

        return descriptions
    #getdbSNPDescriptions

    @srpc(Mandatory.Unicode, Unicode, _returns=GeneLocation)
    def getGeneLocation(gene, build=None):
        """
        Get the location of a gene on the given genome build (assembly), using
        the system's transcript mapping database.

        @arg gene: Gene symbol.
        @type gene: string
        @arg build: Genome build (assembly) by name or alias. If omited,
          the system's default assembly is used.
        @type build: string

        @return: Object with the following fields:
          - gene: Gene symbol.
          - start: Gene start position (one-based, inclusive, in chromosomal
              orientation). If multiple transcripts for the gene are known,
              this contains the lowest start position.
          - stop: Gene stop position (one-based, inclusive, in chromosomal
              orientation). If multiple transcripts for the gene are known,
              this contains the highest stop position.
          - orientation: Gene orientation, either 'forward' or 'reverse'.
          - chromosome_name: Gene chromosome by name (e.g., 'chrX').
          - chromosome_accession: Gene chromosome by accession (e.g.,
              'NC_000023.10').
          - assembly_name: Selected genome build (assembly) by name (e.g.,
              'GRCh37').
          - assembly_alias: Selected genome build (assembly) by alias (e.g.,
              'hg19').
        """
        output = Output(__file__)

        output.addMessage(__file__, -1, 'INFO',
                          'Received request getGeneLocation(%s, %s)'
                          % (gene, build))

        try:
            assembly = Assembly.by_name_or_alias(build or
                                                 settings.DEFAULT_ASSEMBLY)
        except NoResultFound:
            output.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault("EARG",
                        "The build argument (%s) was not a valid " \
                        "build name." % build)

        # From all the transcripts for this gene, get the lowest start
        # position and highest stop position. For integrity, we group by
        # chromosome and orientation.
        # Order by chromosome name for disambiguation, as is done in
        # Convertor._get_mapping()
        mapping = \
            session.query(func.min(TranscriptMapping.start),
                          func.max(TranscriptMapping.stop),
                          TranscriptMapping.orientation,
                          Chromosome.name,
                          Chromosome.accession) \
                   .filter(TranscriptMapping.chromosome.has(assembly=assembly),
                           TranscriptMapping.gene == gene) \
                   .join(TranscriptMapping.chromosome) \
                   .group_by(Chromosome.id,
                             TranscriptMapping.orientation) \
                   .order_by(Chromosome.name.asc()) \
                   .first()

        if not mapping:
            output.addMessage(__file__, 4, "EARG", "EARG %s" % gene)
            raise Fault("EARG",
                        "No location was found for gene %s." % gene)

        start, stop, orientation, chromosome_name, chromosome_accession \
            = mapping

        result = GeneLocation()
        result.gene = gene
        result.start = start
        result.stop = stop
        result.orientation = orientation
        result.chromosome_name = chromosome_name
        result.chromosome_accession = chromosome_accession
        result.assembly_name = assembly.name
        result.assembly_alias = assembly.alias

        output.addMessage(__file__, -1, 'INFO',
                          'Finished processing getGeneLocation(%s %s)'
                          % (gene, build))

        return result
    #getGeneLocation
#MutalyzerService


    @srpc(Mandatory.Unicode, _returns=Array(TranscriptInChromosome))
    def mapTranscriptToChromosomes(transcript):
        """
        Searches for a transcript reference (NM) in the chromosomal database
        and returns the chromosome accession.version numbers list. Example:
            transcript input:
                NM_000267.3
            output:
                [
                    {
                        "transcript": "NM_000267.3",
                        "chromosome": "NC_000017.10"
                    },
                    {
                        "transcript": "NM_000267.3",
                        "chromosome": "NC_000017.11"
                    }
                ]

        @arg transcript: The transcript 'accession[.version]'.
        @type transcript: string

        @return: Complex object:
          - transcript  ; Transcription accesion.version.
          - chromosome   ; Chromosome accesioin.version.
        @rtype: object
        """
        O = Output(__file__)

        O.addMessage(__file__, -1, "INFO",
                     "Received request mapTranscriptToChromosomes(%s)" % transcript)

        from mutalyzer.dbgb.models import Transcript, Reference

        if '.' not in transcript:
            acc, ver = transcript, None
        else:
            acc, ver = transcript.split('.')

        result = []

        if ver is not None:
            transcripts = Transcript.query.\
                filter_by(transcript_accession=acc, transcript_version=ver).all()
        else:
            transcripts = Transcript.query.\
                filter_by(transcript_accession=acc).all()
        for transcript in transcripts:
            reference = Reference.query.\
                filter_by(id=transcript.reference_id).first()
            transcript_output = TranscriptInChromosome()
            transcript_output.transcript = transcript.transcript_accession\
                                           + '.' + transcript.transcript_version
            transcript_output.chromosome = reference.accession\
                                           + '.' + reference.version
            result.append(transcript_output)

        sessiongb.remove()

        O.addMessage(__file__, -1, "INFO",
                     "Finished processing mapTranscriptToChromosomes(%s)" % transcript)
        return result
    # mapTranscriptToChromosomes


# Close database session at end of each call.
def _shutdown_session(ctx):
    session.remove()
MutalyzerService.event_manager.add_listener('method_return_object',
                                            _shutdown_session)
MutalyzerService.event_manager.add_listener('method_exception_object',
                                            _shutdown_session)
