"""
Queries on database models.
"""


# Todo: Some (all) of these are probably better defined as class methods on
#   the models they work with.


from __future__ import unicode_literals

from datetime import datetime, timedelta

from sqlalchemy import and_, or_
import sqlalchemy.exc

from mutalyzer.config import settings
from mutalyzer.db import session
from mutalyzer.db.models import BatchQueueItem, TranscriptProteinLink


def pop_batch_queue_item(batch_job):
    """
    Get the next batch queue item for the given batch job. Return its fields
    as a tuple `item`, `flags` and remove it from the database.

    If no batch queue item could be found for this batch job, return `None`.

    .. note:: Originally, finding the next batch queue item was done using a
        more complicated query::

            SELECT QueueID, Input, Flags
            FROM BatchQueue
            WHERE QueueID = (
                SELECT MIN(QueueID)
                FROM BatchQueue
                GROUP BY JobID
                HAVING JobID = {batch_job.id}
            );

        However, I couldn't see any significant performance difference in my
        latest benchmarks, so we stick with the more obvious query for now.
    """
    batch_queue_item = BatchQueueItem.query \
        .filter_by(batch_job=batch_job) \
        .order_by(BatchQueueItem.id.asc()) \
        .first()
    if batch_queue_item is None:
        return None

    item, flags = batch_queue_item.item, batch_queue_item.flags

    session.delete(batch_queue_item)
    session.commit()

    return item, flags


def get_transcript_protein_link(accession, reverse=False):
    """
    Get a cached link between a transcript and a protein that is not expired
    according to the configuration settings `PROTEIN_LINK_EXPIRATION` and
    `NEGATIVE_PROTEIN_LINK_EXPIRATION`.

    :arg str accession: Accession number (without version number) to lookup
      link for.
    :arg bool reverse: If `True`, `accession` is assumed to be a protein
      accession number, otherwise `accession` is assumed to be a transcript
      accession number.

    Note that the link may be negative, i.e., the knowledge that no link
    exists can also be cached. In that case, the `protein_accession` field of
    the resulting `TranscriptProteinLink` object is `None`.

    Returns `None` if no link (positive or negative) is found.
    """
    link_datetime = datetime.now() - \
        timedelta(seconds=settings.PROTEIN_LINK_EXPIRATION)
    negative_link_datetime = datetime.now() - \
        timedelta(seconds=settings.NEGATIVE_PROTEIN_LINK_EXPIRATION)

    # Query column must have `accession`, other column has the value we're
    # probably interested in.
    query_column = TranscriptProteinLink.transcript_accession
    other_column = TranscriptProteinLink.protein_accession

    if reverse:
        # Lookup by protein accession instead of transcript accession.
        query_column, other_column = other_column, query_column

    return TranscriptProteinLink.query.filter(
        query_column == accession,
        or_(and_(other_column.isnot(None),
                 TranscriptProteinLink.added >= link_datetime),
            and_(other_column.is_(None),
                 TranscriptProteinLink.added >= negative_link_datetime))
    ).first()


def update_transcript_protein_link(transcript_accession=None,
                                   protein_accession=None):
    """
    Update cached link between a transcript and a protein, or create it if it
    doesn't exist yet.

    :arg str transcript_accession: Transcript accession number (without
      version number).
    :arg str protein_accession: Protein accession number (without version
      number).

    At least one of `transcript_accession` or `protein_accession` must be not
    `None`.
    """
    if transcript_accession is None and protein_accession is None:
        raise ValueError('Link must have a transcript or protein')

    # Filter clauses to find links for either of the given accession numbers.
    clauses = []
    if transcript_accession is not None:
        clauses.append(TranscriptProteinLink.transcript_accession ==
                       transcript_accession)
    if protein_accession is not None:
        clauses.append(TranscriptProteinLink.protein_accession ==
                       protein_accession)

    # Delete any related existing links.
    TranscriptProteinLink.query.filter(or_(*clauses)).delete()
    session.commit()

    # There is a race condition here between deleting old links and adding the
    # new one. It's extremely unlikely to go wrong, and we can safely ignore
    # it anyway.
    link = TranscriptProteinLink(transcript_accession, protein_accession)
    try:
        session.add(link)
        session.commit()
    except sqlalchemy.exc.IntegrityError:
        session.rollback()
