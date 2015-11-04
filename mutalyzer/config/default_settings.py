"""
Default Mutalyzer settings. Override these with settings in the module
pointed-to by the `MUTALYZER_SETTINGS` environment variable.
"""


from __future__ import unicode_literals


# Use Mutalyzer in debug mode.
DEBUG = False

# We are running unit tests.
TESTING = False

# This address is used in contact information on the website, as sender in
# batch job notifications, and with retrieval of records at the NCBI using
# Entrez.
EMAIL = 'mutalyzer@humgen.nl'

# The cache directory. Used to store uploaded and downloaded files (e.g.,
# reference files from NCBI or user) and batch job results.
CACHE_DIR = '/tmp'

# Maximum size for uploaded and downloaded files (in bytes).
MAX_FILE_SIZE = 10 * 1048576 # 10 MB

# Maximum sequence length for description extractor (in bases).
EXTRACTOR_MAX_INPUT_LENGTH = 50 * 1000 # 50 Kbp

# The WSGI application runs behind a reverse proxy (e.g., nginx using
# proxy_pass). This needs to be set if the application is mapped to a URL
# other than / or a different HTTP scheme is used by the reverse proxy.
# http://flask.pocoo.org/snippets/35/
REVERSE_PROXIED = False

# Redis connection URI (can be any redis-py connection URI). Set to `None` to
# silently use a mock Redis. Redis is only used for non-essential features.
REDIS_URI = None

# Database connection URI (can be any SQLAlchemy connection URI).
DATABASE_URI = 'sqlite://'

# Default genome assembly (by name or alias).
DEFAULT_ASSEMBLY = 'hg19'

# Name and location of the log file.
LOG_FILE = '/tmp/mutalyzer.log'

# Level of logged messages.
LOG_LEVEL = 3

# Level of output messages.
OUTPUT_LEVEL = 1

# Format of time prefix for log messages. Can be anything that is accepted as
# the format argument of time.strftime.
# http://docs.python.org/2/library/time.html#time.strftime
LOG_TIME_FORMAT = "%Y-%m-%d %H:%M:%S"

# Prefix URL from where LRG files are fetched.
#LRG_PREFIX_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/lrgex/'
LRG_PREFIX_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/lrgex/SCHEMA_1_7_ARCHIVE/'

# Allow for this fraction of errors in batch jobs.
BATCH_JOBS_ERROR_THRESHOLD = 0.05

# Cache expiration time for negative transcript<->protein links from the NCBI
# (in seconds).
NEGATIVE_LINK_CACHE_EXPIRATION = 60 * 60 * 24 * 30

# URL to the website root (without trailing slash). Used for generating
# download links in the batch scheduler.
WEBSITE_ROOT_URL = None

# URL to the SOAP webservice WSDL document. Used to build the WSDL document
# and for linking to it from the documentation page on the website.
SOAP_WSDL_URL = None

# URL to the HTTP/RPC+JSON webservice root (without trailing slash). Used for
# linking to it from the documentation page on the website.
JSON_ROOT_URL = None

# Is Piwik enabled?
PIWIK = False

# Base URL for the Piwik server.
PIWIK_BASE_URL = 'https://piwik.example.com'

# Piwik site ID for Mutalyzer.
PIWIK_SITE_ID = 1

# Enable the Werkzeug reloader for the website.
# This is disabled by default due to a bug with using the reloader in
# combination with `python -m mutalyzer.entrypoints.website`.
# https://github.com/mitsuhiko/werkzeug/issues/461#issuecomment-139369694
USE_RELOADER = False
