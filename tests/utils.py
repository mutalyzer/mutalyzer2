"""
Utilities for unit tests.
"""


from __future__ import unicode_literals

from functools import wraps
import os
import shutil
import tempfile

from mutalyzer.config import settings
from mutalyzer.redisclient import client as redis
from mutalyzer import db


class TestEnvironment(object):
    """
    Configure Mutalyzer for unit tests. All storage is transient and isolated.
    """
    def __init__(self, fixtures=None):
        fixtures = fixtures or []

        self.cache_dir = tempfile.mkdtemp()

        log_handle, self.log_file = tempfile.mkstemp()
        os.close(log_handle)

        database_uri = os.getenv('MUTALYZER_TEST_DATABASE_URI', 'sqlite://')
        redis_uri = os.getenv('MUTALYZER_TEST_REDIS_URI', None)

        settings.configure({'DEBUG':        False,
                            'TESTING':      True,
                            'CACHE_DIR':    self.cache_dir,
                            'REDIS_URI':    redis_uri,
                            'DATABASE_URI': database_uri,
                            'LOG_FILE':     self.log_file})

        # Mutalyzer create tables automatically if we're using an SQLite
        # in-memory database.
        if database_uri != 'sqlite://':
            db.Base.metadata.drop_all(db.session.get_bind())
            db.Base.metadata.create_all(db.session.get_bind())

        if redis_uri is not None:
            redis.flushdb()

        for fixture in fixtures:
            fixture()

    def destroy(self):
        """
        Destroy all storage defined in the current environment.
        """
        db.session.remove()

        shutil.rmtree(self.cache_dir)
        os.unlink(self.log_file)


class MutalyzerTest(object):
    """
    Test class providing an isolated test environment for each test.
    """
    fixtures = ()

    def setup(self):
        self.environment = TestEnvironment(fixtures=self.fixtures)

    def teardown(self):
        self.environment.destroy()


def fix(*fixtures):
    """
    Decorator for a unit test setting up the specified fixtures.
    """
    def decorator(f):
        @wraps(f)
        def fixed_f(*args, **kwargs):
            for fixture in fixtures:
                fixture()
            return f(*args, **kwargs)
        return fixed_f
    return decorator
