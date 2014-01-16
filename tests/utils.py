"""
Utilities for unit tests.
"""


from functools import wraps
import os
import shutil
import tempfile

from mutalyzer.config import settings


class TestEnvironment(object):
    """
    Configure Mutalyzer for unit tests. All storage is transient and isolated.
    """
    def __init__(self, fixtures=None):
        fixtures = fixtures or []

        self.cache_dir = tempfile.mkdtemp()

        log_handle, self.log_file = tempfile.mkstemp()
        os.close(log_handle)

        settings.configure({'DEBUG':        False,
                            'TESTING':      True,
                            'CACHE_DIR':    self.cache_dir,
                            'REDIS_URI':    None,
                            'DATABASE_URI': 'sqlite://',
                            'LOG_FILE':     self.log_file})

        for fixture in fixtures:
            fixture()

    def destroy(self):
        """
        Destroy all storage defined in the current environment.
        """
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
