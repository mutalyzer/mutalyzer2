"""
Test configuration.
"""


from __future__ import unicode_literals

from fixtures import *  # noqa


DEFAULT_DATABASE_URIS = ['sqlite://']
DEFAULT_REDIS_URI = None


def pytest_addoption(parser):
    parser.addoption(
        '--database-uri', metavar='URI', dest='database_uris', default=[],
        action='append',
        help='Database connection, multiple allowed (default: in-memory '
        'SQLite database)')
    parser.addoption(
        '--redis-uri', metavar='URI', dest='redis_uri',
        default=DEFAULT_REDIS_URI,
        help='Redis connection (default: mock Redis server)')


def pytest_generate_tests(metafunc):
    if 'database_uri' in metafunc.fixturenames:
        metafunc.parametrize(
            'database_uri',
            metafunc.config.option.database_uris or DEFAULT_DATABASE_URIS)
