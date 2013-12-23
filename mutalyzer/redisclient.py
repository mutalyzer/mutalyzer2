"""
Lazy global interface to Redis.

Redis connections can safely be shared among threads, so we can keep this very
simple and just use one global connection pool as created by `StrictRedis`.

.. note:: Redis is only used for non-essential features in Mutalyzer and
    therefore a Redis server is not a hard requirement to run Mutalyzer.

    If the `REDIS_URI` configuration setting is `None`, we silently
    instantiate a mock interface to Redis.
"""


import redis

from mutalyzer import settings
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


client = LazyCLient()
