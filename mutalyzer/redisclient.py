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


from __future__ import unicode_literals

import redis

from mutalyzer.config import settings
from mutalyzer import util


class LazyClient(util.LazyObject):
    """
    A lazy proxy for a `StrictRedis` object.
    """
    def _setup(self):
        """
        Instantiate the Redis interface. This is called the first time Redis
        is used.
        """
        self.configure(settings.REDIS_URI)

    def configure(self, uri):
        """
        Configure the client with a new connectionpool for the given URI.
        """
        if uri is None:
            import mockredis
            self._wrapped = mockredis.MockRedis(strict=True)
        else:
            self._wrapped = redis.StrictRedis.from_url(uri,
                                                       decode_responses=True,
                                                       encoding='utf-8')


def configure_client(uri):
    """
    (Re)configure the client with the given URI.
    """
    global client
    client.configure(uri)


# Reconfigure the client if configuration is updated.
settings.on_update(configure_client, 'REDIS_URI')


#: Global :class:`LazyClient` instance. Use this for all communication with
#: Redis.
client = LazyClient()
