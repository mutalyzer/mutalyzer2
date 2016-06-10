"""
Models backed by SQL using SQLAlchemy.
"""


from __future__ import unicode_literals

from datetime import datetime
import sqlite3
import uuid

import binning
from sqlalchemy import event, or_
from sqlalchemy import (Boolean, Column, DateTime, Enum, ForeignKey, Index,
                        Integer, String, Text, TypeDecorator)
from sqlalchemy.engine import Engine
from sqlalchemy.orm import backref, relationship

from mutalyzer import db


BATCH_JOB_TYPES = ('name-checker',
                   'syntax-checker',
                   'position-converter',
                   'snp-converter')


@event.listens_for(Engine, 'connect')
def set_sqlite_pragma(dbapi_connection, connection_record):
    """
    We use foreign keys (and ``ON DELETE CASCADE`` on some of these), but in
    SQLite these are only enforced if ``PRAGMA foreign_keys=ON`` is executed
    on all connections before use.

    [1] http://docs.sqlalchemy.org/en/latest/dialects/sqlite.html#foreign-key-support
    """
    if isinstance(dbapi_connection, sqlite3.Connection):
        cursor = dbapi_connection.cursor()
        cursor.execute('PRAGMA foreign_keys=ON')
        cursor.close()


class Positions(TypeDecorator):
    """
    Represents an immutable list of integers as a concatentation of their
    string serializations.

    Adapted from the `Marshal JSON Strings
    <http://docs.sqlalchemy.org/en/latest/core/types.html#marshal-json-strings>`_
    example in the SQLAlchemy documentation.
    """
    impl = Text

    def process_bind_param(self, value, dialect):
        if value is not None:
            value = ','.join(unicode(i) for i in value)
        return value

    def process_result_value(self, value, dialect):
        if value is None:
            return None
        if value == '':
            return []
        return [int(s) for s in value.split(',')]


class BatchJob(db.Base):
    """
    Batch job.
    """
    __tablename__ = 'batch_jobs'
    __table_args__ = {'mysql_engine': 'InnoDB', 'mysql_charset': 'utf8'}

    id = Column(Integer, primary_key=True)

    #: Email address of user who submitted the job.
    email = Column(String(200))

    #: Type of batch job.
    job_type = Column(Enum(*BATCH_JOB_TYPES, name='job_type'), nullable=False)

    #: Optional argument (currently only used when `job_type` is
    #: ``PositionConverter``, where it denotes the assembly).
    argument = Column(String(20))

    #: Identifier to use in the job result filename and thus the URL for
    #: downloading it. We don't use the auto-incrementing `id` field for this,
    #: since it can be guessed by any user.
    result_id = Column(String(50), nullable=False, index=True, unique=True)

    #: Date and time of creation.
    added = Column(DateTime)

    def __init__(self, job_type, email=None, argument=None):
        self.job_type = job_type
        self.email = email
        self.argument = argument
        self.result_id = unicode(uuid.uuid4())
        self.added = datetime.now()

    def __repr__(self):
        return '<BatchJob %r result_id=%r email=%r>' \
            % (self.job_type, self.result_id, self.email)


class BatchQueueItem(db.Base):
    """
    The smallest unit of work within a batch job.
    """
    __tablename__ = 'batch_queue_items'
    __table_args__ = {'mysql_engine': 'InnoDB', 'mysql_charset': 'utf8'}

    id = Column(Integer, primary_key=True)
    batch_job_id = Column(Integer,
                          ForeignKey('batch_jobs.id', ondelete='CASCADE'),
                          nullable=False)

    #: Input line from a batch job.
    item = Column(String(200), nullable=False)

    #: A list of flags set for this item. A flag consists of either an ``A``,
    #: ``S`` or ``C`` followed by a digit, which refers to the reason of
    #: alteration/skip. We simply store the concatenation of these flags.
    flags = Column(String(20), nullable=False)

    #: The :class:`BatchJob` for this item.
    batch_job = relationship(
        BatchJob,
        backref=backref('batch_queue_items', lazy='dynamic',
                        cascade='all, delete-orphan',
                        passive_deletes=True))

    def __init__(self, batch_job, item, flags=None):
        self.batch_job = batch_job
        self.item = item
        self.flags = flags or ''

    def __repr__(self):
        return '<BatchQueueItem %r flags=%r>' % (self.item, self.flags)


Index('batch_queue_item_with_batch_job',
      BatchQueueItem.batch_job_id, BatchQueueItem.id)


class Reference(db.Base):
    """
    Cached information about a reference sequence.
    """
    __tablename__ = 'references'
    __table_args__ = {'mysql_engine': 'InnoDB', 'mysql_charset': 'utf8'}

    id = Column(Integer, primary_key=True)

    #: Accession number for this reference, including the version number if
    #: applicable (e.g., ``AL449423.14``, ``NM_000059.3``,
    #: ``UD_138781341344``).
    accession = Column(String(20), nullable=False, index=True, unique=True)

    #: MD5 checksum of the reference file.
    checksum = Column(String(32), nullable=False, index=True, unique=True)

    #: Source of the reference, needed to re-fetch the reference file. Meaning
    #: of the allowed values:
    #: - ``ncbi``: NCBI nuccore database.
    #: - ``ncbi_slice``: NCBI nuccore database by slicing. `source_data` is a
    #:   string of the form ``<accession>:<start>:<stop>:<orientation>`` with
    #:   one-based, inclusive start and stop position (in reference
    #:   orientation) and ``forwad`` or ``reverse`` for ``<orientation>``.
    #: - ``lrg``: Official LRG repository.
    #: - ``url``: Custom HTTP/HTTPS/FTP URL. `source_data` contains the URL.
    #: - `upload``: Custom file upload. There's no way to re-fetch this file.
    source = Column(Enum('ncbi', 'ncbi_slice', 'lrg', 'url', 'upload',
                         name='reference_source'),
                    nullable=False)

    #: Additional data needed to re-fetch the reference file. The data depends
    #: on the value of `source` and must be serialized as a string.
    source_data = Column(String(255))

    #: The corresponding GI number, if available.
    geninfo_identifier = Column(String(13), index=True, unique=True)

    #: Date and time of creation.
    added = Column(DateTime)

    def __init__(self, accession, checksum, source, source_data=None,
                 geninfo_identifier=None):
        self.accession = accession
        self.checksum = checksum
        self.source = source
        self.source_data = source_data
        self.geninfo_identifier = geninfo_identifier
        self.added = datetime.now()

    def __repr__(self):
        return '<Reference %r>' % self.accession


Index('reference_source_data',
      Reference.source, Reference.source_data)


class Assembly(db.Base):
    """
    Genome assembly.
    """
    __tablename__ = 'assemblies'
    __table_args__ = {'mysql_engine': 'InnoDB', 'mysql_charset': 'utf8'}

    id = Column(Integer, primary_key=True)

    #: Genome Reference Consortium name (e.g., ``GRCh37``, ``GRCm38``).
    name = Column(String(30), nullable=False, unique=True)

    #: Short alias (e.g., ``hg19``, ``mm10``).
    alias = Column(String(10), unique=True)

    #: NCBI taxonomy identifier (e.g., 9606, 10090).
    taxonomy_id = Column(Integer, nullable=False)

    #: NCBI taxonomy name (e.g., ``Homo sapiens``, ``Mus musculus``).
    taxonomy_common_name = Column(String(50), nullable=False)

    #: Criteria to order assemblies by in user-visible lists.
    order_by_criteria = taxonomy_common_name.asc(), alias.asc()

    def __init__(self, name, taxonomy_id, taxonomy_common_name, alias=None):
        self.name = name
        self.taxonomy_id = taxonomy_id
        self.taxonomy_common_name = taxonomy_common_name
        self.alias = alias

    def __repr__(self):
        return '<Assembly %r taxonomy=%r>' % (self.name,
                                              self.taxonomy_common_name)

    @classmethod
    def by_name_or_alias(cls, name_or_alias):
        """
        Returns an :class:`Assembly` instance for the given `name_or_alias` if
        it exists, or raises :exc:`sqlalchemy.orm.exc.NoResultFound`
        otherwise.
        """
        return cls.query.filter(or_(cls.name == name_or_alias,
                                    cls.alias == name_or_alias)).one()


class Chromosome(db.Base):
    """
    Chromosome name and accession number in a genome assembly.
    """
    __tablename__ = 'chromosomes'
    __table_args__ = {'mysql_engine': 'InnoDB', 'mysql_charset': 'utf8'}

    id = Column(Integer, primary_key=True)
    assembly_id = Column(Integer,
                         ForeignKey('assemblies.id', ondelete='CASCADE'),
                         nullable=False)

    #: Chromosome name (e.g., ``chr1``, ``chrM``).
    name = Column(String(30), nullable=False)

    #: Chromosome accession number (e.g., ``NC_000001.10``, ``NC_012920.1``).
    accession = Column(String(30), nullable=False)

    # Todo: Use constant definitions for this type.
    #: Enclosing organelle. Used to differentiate between references requiring
    #: ``m.`` and ``g.`` descriptions.
    organelle = Column(Enum('nucleus', 'mitochondrion',
                            name='organelle'))

    #: The :class:`Assembly` this chromosome is in.
    assembly = relationship(
        Assembly,
        lazy='joined',
        innerjoin=True,
        backref=backref('chromosomes', lazy='dynamic',
                        cascade='all, delete-orphan',
                        passive_deletes=True))

    def __init__(self, assembly, name, accession, organelle):
        self.assembly = assembly
        self.name = name
        self.accession = accession
        self.organelle = organelle

    def __repr__(self):
        return '<Chromosome %r accession=%r>' % (self.name, self.accession)


Index('chromosome_name',
      Chromosome.assembly_id, Chromosome.name,
      unique=True)
Index('chromosome_accession',
      Chromosome.assembly_id, Chromosome.accession,
      unique=True)


class TranscriptMapping(db.Base):
    """
    Mapping of a gene transcript on a chromosome.

    .. note:: All positions are ordered according to the chromosome,
       irrespective of transcript orientation. For example, `start` is always
       smaller than `stop`.
    """
    __tablename__ = 'transcript_mappings'
    __table_args__ = {'mysql_engine': 'InnoDB', 'mysql_charset': 'utf8'}

    id = Column(Integer, primary_key=True)
    chromosome_id = Column(Integer,
                           ForeignKey('chromosomes.id', ondelete='CASCADE'),
                           nullable=False)

    #: Type of reference sequence containing the transcript.
    reference_type = Column(Enum('refseq', 'lrg', name='reference_type'),
                            nullable=False)

    #: Accession number for the transcript reference, not including the
    #: version number (e.g., ``NM_000020``, ``NC_012920``, ``LRG_14``).
    accession = Column(String(20), nullable=False)

    #: Accession version (e.g., 3, 2). Not applicaple when `reference_type` is
    #: ``lrg``.
    version = Column(Integer)

    #: Gene symbol (e.g., ``DMD``, ``PSEN1``, ``TRNS1``).
    gene = Column(String(30), nullable=False)

    #: Transcript number in reference (e.g., 1, 2).
    transcript = Column(Integer, nullable=False)

    # Todo: Use constants and some utilities for conversion.
    #: The orientation of the transcript on the chromosome.
    orientation = Column(Enum('forward', 'reverse', name='orientation'),
                         nullable=False)

    #: The start position of the transcript on the chromosome (one-based,
    #: inclusive, in chromosomal orientation).
    start = Column(Integer, nullable=False)

    #: The stop position of the transcript on the chromosome (one-based,
    #: inclusive, in chromosomal orientation).
    stop = Column(Integer, nullable=False)

    #: Bin index that can be used for faster range-based queries. See the
    #: `interval binning documentation <http://interval-binning.readthedocs.org/>`_
    #: for more information.
    bin = Column(Integer, nullable=False, index=True)

    #: The CDS start position of the transcript on the chromosome (one-based,
    #: inclusive, in chromosomal orientation).
    cds_start = Column(Integer)

    #: The CDS stop position of the transcript on the chromosome (one-based,
    #: inclusive).
    cds_stop = Column(Integer)

    #: The exon start positions of the transcript on the chromosome
    #: (one-based, inclusive, in chromosomal orientation).
    exon_starts = Column(Positions, nullable=False)

    #: The exon start positions of the transcript on the chromosome
    #: (one-based, inclusive, in chromosomal orientation).
    exon_stops = Column(Positions, nullable=False)

    #: If `False`, variant descriptions can use just the accession number
    #: without gene and transcript selector (e.g., ``NM_000020:c.1del``,
    #: ``LRG_1:c.1del``). If `True`, gene and transcript selection is
    #: necessary (e.g., ``NC_012920(TRNI_v001):c.1del``, ``LRG_1t1:c.1del``).
    select_transcript = Column(Boolean, nullable=False)

    #: Source of the mapping.
    source = Column(Enum('ucsc', 'ncbi', 'ebi', 'reference', name='source'),
                    nullable=False)

    #: The :class:`Assembly` this chromosome is in.
    chromosome = relationship(
        Chromosome,
        lazy='joined',
        innerjoin=True,
        backref=backref('transcript_mappings', lazy='dynamic',
                        cascade='all, delete-orphan',
                        passive_deletes=True))

    def __init__(self, chromosome, reference_type, accession, gene,
                 orientation, start, stop, exon_starts, exon_stops, source,
                 transcript=1, cds=None, select_transcript=False,
                 version=None):
        self.chromosome = chromosome
        self.reference_type = reference_type
        self.accession = accession
        self.gene = gene
        self.orientation = orientation
        self.start = start
        self.stop = stop
        self.exon_starts = exon_starts
        self.exon_stops = exon_stops
        self.source = source
        self.transcript = transcript
        self.cds = cds
        self.select_transcript = select_transcript
        self.version = version
        self.bin = binning.assign_bin(self.start - 1, self.stop)

    def __repr__(self):
        return ('<TranscriptMapping accession=%r version=%r gene=%r '
                'transcript=%r>'
                % (self.accession, self.version, self.gene, self.transcript))

    @classmethod
    def create_or_update(cls, chromosome, reference_type, accession, gene,
                         orientation, start, stop, exon_starts, exon_stops,
                         source, transcript=1, cds=None,
                         select_transcript=False, version=None):
        """
        Returns an :class:`TranscriptMapping` instance with the given values.
        If a row with a duplicate key already exists, it is deleted first.

        .. note:: Unfortunately, SQLAlchemy does not have `ON DUPLICATE KEY
            UPDATE` functionality, which would performance-wise be the best
            way to update transcript mappings. Even if it had, PostgreSQL does
            not. This class method implements an alternative albeit at the
            cost of querying the table for an existing entry on each update.

            From PostgreSQL 9.5 we do have ``INSERT ... ON CONFLICT DO UPDATE``
            which we might want to use in the future.

            http://stackoverflow.com/questions/17267417/how-do-i-do-an-upsert-merge-insert-on-duplicate-update-in-postgresql
        """
        # Actually we should do a `lock transcript_mappings in exclusive mode;`
        # to prevent concurrent reads to miss this entry.
        cls.query.filter_by(
            chromosome=chromosome, accession=accession, version=version,
            gene=gene, transcript=transcript
        ).delete()
        return cls(chromosome, reference_type, accession, gene, orientation,
                   start, stop, exon_starts, exon_stops, source,
                   transcript=transcript, cds=cds,
                   select_transcript=select_transcript, version=version)

    @property
    def coding(self):
        """
        Set to `True` iff the transcript is coding.
        """
        return self.cds_start is not None and self.cds_stop is not None

    @property
    def cds(self):
        """
        Tuple of CDS start and stop positions on the chromosome, or `None` if
        the transcript is non-coding.
        """
        if self.coding:
            return self.cds_start, self.cds_stop
        return None

    @cds.setter
    def cds(self, cds):
        self.cds_start, self.cds_stop = cds or (None, None)

    def get_reference(self, include_version=True):
        """
        Get fully qualified reference for this transcript.

        You would usually want to use the simpler :attr:`reference` property
        instead, except if the accession number may not include version number
        (which we consider bad practice).
        """
        if include_version and self.version:
            accession = '%s.%i' % (self.accession, self.version)
        else:
            accession = self.accession

        if self.select_transcript:
            if self.reference_type == 'lrg':
                selector = 't%d' % self.transcript
            elif self.transcript:
                selector = '(%s_v%.3i)' % (self.gene, self.transcript)
            else:
                selector = '(%s)' % self.gene
        else:
            selector = ''

        return '%s%s' % (accession, selector)

    @property
    def reference(self):
        """
        Fully qualified reference for this transcript.
        """
        return self.get_reference()


Index('transcript_mapping_transcript',
      TranscriptMapping.accession, TranscriptMapping.version,
      TranscriptMapping.gene, TranscriptMapping.transcript,
      TranscriptMapping.chromosome_id,
      unique=True)
