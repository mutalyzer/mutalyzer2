"""
Models backed by SQL using SQLAlchemy.
"""


from datetime import datetime
import sqlite3
import uuid

from sqlalchemy import (event, Column, Boolean, BigInteger, DateTime, ForeignKey,
                        Integer, Numeric, String, Table, Text, Index, Enum,
                        UniqueConstraint)
from sqlalchemy.engine import Engine
from sqlalchemy.orm import relationship, backref

from mutalyzer import db


BATCH_JOB_TYPES = ('NameChecker', 'SyntaxChecker', 'PositionConverter',
                   'SnpConverter')


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


class BatchJob(db.Base):
    """
    Batch job.
    """
    __tablename__ = 'batch_jobs'
    __table_args__ = {'mysql_engine': 'InnoDB', 'mysql_charset': 'utf8'}

    id = Column(Integer, primary_key=True)

    #: Email address of user who submitted the job.
    email = Column(String(200))

    #: Prefix of URL for downloading the job result file. This would usually
    #: be the Mutalyzer website base URL.
    download_url_prefix = Column(String(200))

    #: Type of batch job.
    job_type = Column(Enum(*BATCH_JOB_TYPES, name='job_type'), nullable=False)

    #: Optional argument (currently only used with job type PositionConverter,
    #: where it denotes the assembly).
    argument = Column(String(20))

    #: Identifier to use in the job result filename and thus the URL for
    #: downloading it. We don't use the auto-incrementing `id` field for this,
    #: since it can be guessed by any user.
    result_id = Column(String(50), nullable=False, index=True, unique=True)

    #: Date and time of creation.
    added = Column(DateTime)

    def __init__(self, job_type, email=None, download_url_prefix=None,
                 argument=None):
        self.job_type = job_type
        self.email = email
        self.download_url_prefix = download_url_prefix
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

    #: A list of flags set for this item. A flag consists of either an A, S or
    #: C followed by a digit, which refers to the reason of alteration/skip.
    #: We simply store the concatenation of these flags.
    flags = Column(String(20), nullable=False)

    #: The :class:`BatchJob` for this item.
    batch_job = relationship(
        BatchJob,
        backref=backref('batch_jobs', lazy='dynamic',
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
    #: applicable (e.g., AL449423.14, NM_000059.3, UD_138781341344).
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
    #: (e.g., NM_018195, XM_005270562, NR_015380).
    transcript_accession = Column(String(20), nullable=False, index=True,
                                  unique=True)

    #: Accession number for the protein, not including the version number
    #: (e.g., NP_060665, XP_005258635). If `NULL`, the record states that no
    #: protein is linked to the transcript by the NCBI.
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


def create_all():
    db.Base.metadata.drop_all(db.session.get_bind())
    db.Base.metadata.create_all(db.session.get_bind())

    # Todo: Use alembic.

    # if using alembic:
    #from alembic.config import Config
    #from alembic import command
    #alembic_cfg = Config("alembic.ini")
    #command.stamp(alembic_cfg, "head")
