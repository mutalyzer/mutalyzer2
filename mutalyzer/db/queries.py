"""
Queries on database models.
"""


from mutalyzer.db import session
from mutalyzer.db.models import BatchQueueItem


def pop_batch_queue_item(batch_job):
    """
    Get the next batch queue item for the given batch job. Return its fields
    as a tuple `item`, `flags` and remove it from the database.

    If not batch queue item could be found for this batch job, return `None`.

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
