"""
Queries on database models.
"""


# Todo: Some (all) of these are probably better defined as class methods on
#   the models they work with.


from __future__ import unicode_literals

from datetime import datetime, timedelta

from sqlalchemy import and_, or_

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


def get_transcript_protein_link(transcript_accession):
    """
    Get a cached link between a transcript and a protein that is not expired
    according to the configuration settings `PROTEIN_LINK_EXPIRATION` and
    `NEGATIVE_PROTEIN_LINK_EXPIRATION`.

    Note that the link may be negative, i.e., the knowledge that no link
    exists can also be cached. In that case, the `protein_accession` field of
    the resulting `TranscriptProteinLink` object is `None`.

    Returns `None` if no link (positive or negative) is found.
    """
    link_datetime = datetime.now() - \
        timedelta(seconds=settings.PROTEIN_LINK_EXPIRATION)
    negative_link_datetime = datetime.now() - \
        timedelta(seconds=settings.NEGATIVE_PROTEIN_LINK_EXPIRATION)

    return TranscriptProteinLink.query \
        .filter_by(transcript_accession=transcript_accession) \
        .filter(or_(
          and_(TranscriptProteinLink.protein_accession != None,
               TranscriptProteinLink.added >= link_datetime),
          and_(TranscriptProteinLink.protein_accession == None,
               TranscriptProteinLink.added >= negative_link_datetime))) \
        .first()


def update_transcript_protein_link(transcript_accession,
                                   protein_accession=None):
    """
    Update cached link between a transcript and a protein, or create it if it
    doesn't exist yet.
    """
    link = TranscriptProteinLink.query \
        .filter_by(transcript_accession=transcript_accession) \
        .first()

    if link is not None:
        link.protein_accession = protein_accession
        link.added = datetime.now()
    else:
        link = TranscriptProteinLink(transcript_accession, protein_accession)
        session.add(link)

    session.commit()
