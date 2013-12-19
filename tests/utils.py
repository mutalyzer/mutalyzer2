import os
import tempfile


log_handle, log_filename = tempfile.mkstemp()
os.close(log_handle)


TEST_SETTINGS = dict(
    DEBUG     = True,
    TESTING   = True,
    CACHE_DIR = tempfile.mkdtemp(),
    LOG_FILE  = log_filename
)
