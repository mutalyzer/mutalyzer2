"""
Default Mutalyzer settings. Override these with settings in the module
pointed-to by the `MUTALYZER_SETTINGS` environment variable.
"""


# Use Mutalyzer in debug mode.
DEBUG = True

# This address is used in contact information on the website, as sender in
# batch job notifications, and with retrieval of records at the NCBI using
# Entrez.
EMAIL = 'mutalyzer@humgen.nl'

# The cache directory. Used to store uploaded and downloaded files (e.g.,
# reference files from NCBI or user) and batch job results.
import tempfile
CACHE_DIR = tempfile.mkdtemp()

# Maximum size of the cache directory (in bytes).
MAX_CACHE_SIZE = 50 * 1048576 # 50 MB

# Maximum size for uploaded and downloaded files (in bytes).
MAX_FILE_SIZE = 10 * 1048576 # 10 MB

# Host name for local MySQL databases.
MYSQL_HOST = 'localhost'

# User for local MySQL databases.
MYSQL_USER = 'mutalyzer'

# Local MySQL database name.
MYSQL_DATABASE = 'mutalyzer'

# Available databases with mapping information.
DB_NAMES = ['hg18', 'hg19', 'mm10']

# Default database for mapping information.
DEFAULT_DB = 'hg19'

# Name and location of the log file.
import os
import tempfile
log_handle, log_filename = tempfile.mkstemp()
os.close(log_handle)
LOG_FILE = log_filename

# Level of logged messages.
LOG_LEVEL = 3

# Level of output messages.
OUTPUT_LEVEL = 1

# Format of time prefix for log messages. Can be anything that is accepted as
# the format argument of time.strftime.
# http://docs.python.org/2/library/time.html#time.strftime
LOG_TIME_FORMAT = "%Y-%m-%d %H:%M:%S"

# Prefix URL from where LRG files are fetched.
LRG_PREFIX_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/lrgex/'

# Allow for this fraction of errors in batch jobs.
BATCH_JOBS_ERROR_THRESHOLD = 0.05

# Number of days a cached transcript->protein link from the NCBI is considered
# valid.
PROTEIN_LINK_LIFETIME = 30

# Number of days a cached nonexisting transcript->protein link from the NCBI
# is considered valid.
PROTEIN_LINK_NONE_LIFETIME = 5

# Is Piwik enabled?
PIWIK = False

# Base URL for the Piwik server.
PIWIK_BASE_URL = 'https://piwik.example.com'

# Piwik site ID for Mutalyzer.
PIWIK_SITE_ID = 1
