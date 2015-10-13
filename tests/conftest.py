"""
Test configuration.
"""


from __future__ import unicode_literals

from disable_network import turn_off_network, turn_on_network
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
    parser.addoption(
        '--allow-network', dest='allow_network', action='store_true',
        help='Allow non-localhost network access')


def pytest_generate_tests(metafunc):
    if 'database_uri' in metafunc.fixturenames:
        metafunc.parametrize(
            'database_uri',
            metafunc.config.option.database_uris or DEFAULT_DATABASE_URIS)


def pytest_configure(config):
    if not config.getoption('allow_network'):
        turn_off_network(verbose=config.option.verbose)


def pytest_unconfigure():
    # This is a no-op if network was never disabled.
    turn_on_network()
