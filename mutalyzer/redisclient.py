"""
Lazy global interface to Redis.

Redis connections can safely be shared among threads, so we can keep this very
simple and just use one global connection pool as created by `StrictRedis`.

.. note:: Redis is only used for non-essential features in Mutalyzer and
    therefore a Redis server is not a hard requirement to run Mutalyzer.

    If the `REDIS_URI` configuration setting is `None`, we silently
    instantiate a mock interface to Redis.

.. todo:: We currently use Redis for storing stat counters, but there are many
    opportunities to use it for caching. For example, which version numbers
    are available for a certain accession number, which is a costly operation
    and something we implemented quite a hack for to optimize in the batch
    name checker (the whole alter thing with batchflags).
"""


import redis

from mutalyzer.config import settings
from mutalyzer import util


class LazyClient(util.LazyObject):
    """
    A lazy proxy for a `StrictRedis` object.
    """
    def _setup(self):
        """
        Instantiate the Redis interface. This called the first time Redis is
        used.
        """
        if settings.REDIS_URI is None:
            import mockredis
            self._wrapped = mockredis.MockRedis(strict=True)
        else:
            self._wrapped = redis.StrictRedis.from_url(settings.REDIS_URI)


client = LazyClient()
