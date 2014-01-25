"""
Models backed by SQL using SQLAlchemy.
"""


from datetime import datetime
import sqlite3
import uuid

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
            value = ','.join(str(i) for i in value)
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

    #: URL for downloading the job result file. This would usually be a view
    #: on the Mutalyzer website.
    download_url = Column(String(200))

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

    def __init__(self, job_type, email=None, download_url=None,
                 argument=None):
        self.job_type = job_type
        self.email = email
        self.download_url = download_url
        self.argument = argument
        self.result_id = str(uuid.uuid4())
        self.added = datetime.now()

    def __repr__(self):
        return '<BatchJob %r filename=%r email=%r>' \
            % (self.job_type, self.filename, self.email)


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

    #: The corresponding GI number, if available.
    geninfo_identifier = Column(String(13), index=True, unique=True)

    #: The accession number from which we took a slice, if available.
    slice_accession = Column(String(20))

    #: The start position on the accession number from which we took a slice,
    #: if available.
    slice_start = Column(Integer)

    #: The stop position on the accession number from which we took a slice,
    #: if available.
    slice_stop = Column(Integer)

    #: The orientation on the accession number from which we took a slice, if
    #: available.
    slice_orientation = Column(Enum('forward', 'reverse',
                                    name='slice_orentation'))

    #: The URL from which the reference file was downloaded, if available.
    download_url = Column(String(255), index=True, unique=True)

    #: Date and time of creation.
    added = Column(DateTime)

    def __init__(self, accession, checksum, geninfo_identifier=None,
                 slice_accession=None, slice_start=None, slice_stop=None,
                 slice_orientation=None, download_url=None):
        self.accession = accession
        self.checksum = checksum
        self.geninfo_identifier = geninfo_identifier
        self.slice_accession = slice_accession
        self.slice_start = slice_start
        self.slice_stop = slice_stop
        self.slice_orientation = slice_orientation
        self.download_url = download_url
        self.added = datetime.now()

    def __repr__(self):
        return '<Reference %r>' % self.accession


Index('reference_slice',
      Reference.slice_accession, Reference.slice_start, Reference.slice_stop,
      Reference.slice_orientation,
      unique=True)


# Todo: Perhaps it is a better fit to implement this with Redis.
class TranscriptProteinLink(db.Base):
    """
    Cached link between a transcript and protein reference.
    """
    __tablename__ = 'transcript_protein_links'
    __table_args__ = {'mysql_engine': 'InnoDB', 'mysql_charset': 'utf8'}

    id = Column(Integer, primary_key=True)

    #: Accession number for the transcript, not including the version number
    #: (e.g., ``NM_018195`, ``XM_005270562``, ``NR_015380``).
    transcript_accession = Column(String(20), nullable=False, index=True,
                                  unique=True)

    #: Accession number for the protein, not including the version number
    #: (e.g., ``NP_060665``, ``XP_005258635``). If `NULL`, the record states
    #: that no protein is linked to the transcript by the NCBI.
    protein_accession = Column(String(20), index=True)

    #: Date and time of creation.
    added = Column(DateTime)

    def __init__(self, transcript_accession, protein_accession=None):
        self.transcript_accession = transcript_accession
        self.protein_accession = protein_accession
        self.added = datetime.now()

    def __repr__(self):
        return '<TranscriptProteinLink transcript=%r protein=%r>' \
            % (self.transcript_accession, self.protein_accession)


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

    #: Type of organelle.
    organelle_type = Column(Enum('chromosome', 'mitochondrion',
                                 name='organelle_type'))

    #: The :class:`Assembly` this chromosome is in.
    assembly = relationship(
        Assembly,
        lazy='joined',
        innerjoin=True,
        backref=backref('chromosomes', lazy='dynamic',
                        cascade='all, delete-orphan',
                        passive_deletes=True))

    def __init__(self, assembly, name, accession, organelle_type):
        self.assembly = assembly
        self.name = name
        self.accession = accession
        self.organelle_type = organelle_type

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
    orientation = Column(Enum('forward', 'reverse', name='orentation'),
                         nullable=False)

    #: The start position of the transcript on the chromosome.
    start = Column(Integer, nullable=False)

    #: The stop position of the transcript on the chromosome.
    stop = Column(Integer, nullable=False)

    #: The CDS start position of the transcript on the chromosome.
    cds_start = Column(Integer)

    #: The CDS stop position of the transcript on the chromosome.
    cds_stop = Column(Integer)

    #: The exon start positions of the transcript on the chromosome.
    exon_starts = Column(Positions, nullable=False)

    #: The exon start positions of the transcript on the chromosome.
    exon_stops = Column(Positions, nullable=False)

    #: If `False`, variant descriptions can use just the accession number
    #: without gene and transcript selector (e.g., ``NM_000020:c.1del``,
    #: ``LRG_1:c.1del``). If `True`, gene and transcript selection is
    #: necessary (e.g., ``NC_012920(TRNI_v001):c.1del``, ``LRG_1_t1:c.1del``).
    select_transcript = Column(Boolean, nullable=False)

    #: Source of the mapping.
    source = Column(Enum('ucsc', 'ncbi', 'reference', name='source'),
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
        Returns an existing :class:`TranscriptMapping` instance with the given
        values for `chromosome`, `accession`, `version`, `gene`, and
        `transcript` if it exists, or a new instance otherwise.

        .. note:: Unfortunately, SQLAlchemy does not have `ON DUPLICATE KEY
            UPDATE` functionality, which would performance-wise be the best
            way to update transcript mappings. This class method implements an
            alternative albeit at the cost of querying the table for an
            existing entry on each update.
        """
        instance = cls.query.filter_by(
            chromosome=chromosome, accession=accession, version=version,
            gene=gene, transcript=transcript).first()
        if instance is None:
            instance = cls(chromosome, reference_type, accession, gene,
                           orientation, start, stop, exon_starts, exon_stops,
                           source, transcript=transcript, cds=cds,
                           select_transcript=select_transcript,
                           version=version)
        else:
            instance.reference_type = reference_type
            instance.orientation = orientation
            instance.start = start
            instance.stop = stop
            instance.exon_starts = exon_starts
            instance.exon_stops = exon_stops
            instance.source = source
            instance.cds = cds
            instance.select_transcript = select_transcript
        return instance

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


Index('transcript_mapping_transcript',
      TranscriptMapping.accession, TranscriptMapping.version,
      TranscriptMapping.gene, TranscriptMapping.transcript,
      TranscriptMapping.chromosome_id,
      unique=True)


def create_all():
    db.Base.metadata.drop_all(db.session.get_bind())
    db.Base.metadata.create_all(db.session.get_bind())
    db.session.commit()

    # Todo: Use alembic.

    # if using alembic:
    #from alembic.config import Config
    #from alembic import command
    #alembic_cfg = Config("alembic.ini")
    #command.stamp(alembic_cfg, "head")
