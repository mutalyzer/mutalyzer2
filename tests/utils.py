import os
import shutil
import tempfile

from mutalyzer.config import settings
from mutalyzer.db import models


def create_test_environment(database=False):
    """
    Configure Mutalyzer for unit tests. All storage is transient and isolated.
    """
    log_handle, log_filename = tempfile.mkstemp()
    os.close(log_handle)

    settings.configure(dict(
            DEBUG        = True,
            TESTING      = True,
            CACHE_DIR    = tempfile.mkdtemp(),
            DATABASE_URI = 'sqlite://',
            LOG_FILE     = log_filename))

    if database:
        models.create_all()


def destroy_environment():
    """
    Destroy all storage defined in the current environment.
    """
    shutil.rmtree(settings.CACHE_DIR)
    os.unlink(settings.LOG_FILE)
