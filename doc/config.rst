.. highlight:: none

.. _config:

Configuration
=============

This section describes how to configure Mutalyzer and includes a list of all
available configuration settings.

Mutalyzer looks for its configuration in the file specified by the
``MUTALYZER_SETTINGS`` environment variable. Make sure to always have this
environment variable set when invoking any component of Mutalyzer. One way of
doing this is by exporting it::

    $ export MUTALYZER_SETTINGS=~/mutalyzer/settings.py

If you like, you can add this command to your ``~/.bashrc`` to have it
executed every time you open a shell.

Another way is by prefixing your invocations with
``MUTALYZER_SETTINGS=...``. For example::

    $ MUTALYZER_SETTINGS=~/mutalyzer/settings.py mutalyzer-website


Example configuration
---------------------

If you followed the steps in :ref:`install`, this is a standard configuration
file that will work for you:

.. code-block:: python

    REDIS_URI = 'redis://localhost'
    DATABASE_URI = 'postgresql://mutalyzer:*****@localhost/mutalyzer'

This is not yet a minimal configuration. In fact, you can run Mutalyzer
without a configuration file since the default configuration works out of the
box. The default configuration uses a mock Redis instance, an in-memory
database backend and a temporary directory for cache and log files, so it is
not recommended for anything more than playing around.

The next section describes all available configuration settings.


Configuration settings
----------------------

Note that the configuration file is interpreted as a Python module, so you can
use arbitrary Python expressions as configuration values, or even import other
modules in it.

Unsetting a configuration setting is done by using the value `None`. If no
default value is mentioned for any configuration setting below it means it is
not set by default.

.. _config-email:

EMAIL
  The email address used in contact information on the website and sent with
  NCBI Entrez calls.

  `Default value:` ``mutalyzer@humgen.nl``

BATCH_NOTIFICATION_EMAIL
  The email address used as sender in batch job notifications. If set to
  `None`, the value of :ref:`EMAIL <config-email>` will be used.

  `Default value:` `None`

.. _config-debug:

DEBUG
  If set to `True`, Mutalyzer runs in debug mode and will show more
  information with errors.

  `Default value:` `False`

.. _config-cache-dir:

CACHE_DIR
  The cache directory which is used to store uploaded and downloaded files
  such as reference files from the NCBI and batch job results.

  `Default value:` ``/tmp``


User input settings
^^^^^^^^^^^^^^^^^^^

MAX_FILE_SIZE
  Maximum size for uploaded and downloaded files (in bytes).

  `Default value:` `10 * 1048576` (10 MB)

EXTRACTOR_MAX_INPUT_LENGTH
  Maximum sequence length for description extractor (in bases).

  `Default value:` `50 * 1000` (50 Kbp)

BATCH_JOBS_ERROR_THRESHOLD
  Allow for this fraction of errors in batch jobs.

  `Default value:` `0.05`


Database settings
^^^^^^^^^^^^^^^^^

DATABASE_URI
  SQLAlchemy database connection URI specifying the database used to store
  users, samples, variants, etcetera.

  ================   =============================================================
  Database system    Example URI
  ================   =============================================================
  PostgreSQL         ``postgresql://mutalyzer:*****@localhost/mutalyzer``
  MySQL              ``mysql://mutalyzer:*****@localhost/mutalyzer?charset=utf8``
  SQLite             ``sqlite:////tmp/mutalyzer.db``
  ================   =============================================================

  See the SQLAlchemy documentation on
  `Engine Configuration
  <http://docs.sqlalchemy.org/en/latest/core/engines.html>`_ for more
  information.

  `Default value:` ``sqlite://`` (in-memory SQLite database)

REDIS_URI
  Redis connection URI (can be any `redis-py
  <https://github.com/andymccurdy/redis-py>`_ connection URI). Set to `None`
  to silently use a mock Redis. Redis is only used for non-essential
  features such as caching of external resources.

  `Default value:` `None`


Settings for output and logging
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All Mutalyzer messages come with a level which can be one of:

======  ========  ======================================================
Level   Alias     Meaning
======  ========  ======================================================
-1      Log       Specifically log a message.
0       Debug     Debug information.
1       Info      Info.
2       Warning   Regular warnings.
3       Error     Serious errors that can be compensated for.
4       Fatal     Errors that are not recoverable.
5       Off       Can be used as a log/output level to turn off output.
======  ========  ======================================================

LOG_FILE
  Name and location of the log file.

  `Default value:` ``/tmp/mutalyzer.log``

LOG_LEVEL
  Level of logged messages.

  `Default value:` `3`

OUTPUT_LEVEL
  Level of output messages.

  `Default value:` `1`

LOG_TIME_FORMAT
  Format of time prefix for log messages. Can be anything that is accepted as
  the format argument of `time.strftime
  <http://docs.python.org/2/library/time.html#time.strftime>`_.

  `Default value:` ``%Y-%m-%d %H:%M:%S``


Website settings
^^^^^^^^^^^^^^^^

.. _config-reverse-proxied:

REVERSE_PROXIED
  If set to `True`, the WSGI application runs behind a reverse proxy (e.g.,
  nginx using ``proxy_pass``). This needs to be set if the application is
  mapped to a URL other than / or a different HTTP scheme is used by the
  reverse proxy.

  `Default value:` `False`

.. _config-website-root-url:

WEBSITE_ROOT_URL
  URL to the website root (without trailing slash). Used for generating
  download links in the batch scheduler.

  `Default value:` `None`

.. _config-soap-wsdl-url:

SOAP_WSDL_URL
  URL to the SOAP webservice WSDL document. Used to build the WSDL document
  and for linking to it from the documentation page on the website.

  `Default value:` `None`

.. _config-json-root-url:

JSON_ROOT_URL
  URL to the HTTP/RPC+JSON webservice root (without trailing slash). Used for
  linking to it from the documentation page on the website.

  `Default value:` `None`


Piwik settings
^^^^^^^^^^^^^^

`Piwik <http://piwik.org/>`_ is an Open Source analytics platform. Mutalyzer
has built-in support for visitor tracking with Piwik.

PIWIK
  If set to `True`, Piwik is enabled and some Javascript tracking code is
  included in every Mutalyzer website page.

  `Default value:` `False`

PIWIK_BASE_URL
  Base URL for the Piwik server.

  `Default value:` ``https://piwik.example.com``

PIWIK_SITE_ID
  Piwik site ID for Mutalyzer.

  `Default value:` `1`


Miscellaneous settings
^^^^^^^^^^^^^^^^^^^^^^

LRG_PREFIX_URL
  Prefix URL from where LRG files are fetched.

  `Default value:` ``ftp://ftp.ebi.ac.uk/pub/databases/lrgex/SCHEMA_1_7_ARCHIVE/``

DEFAULT_ASSEMBLY
  Default genome assembly (by name or alias).

  `Default value:` ``hg19``

NEGATIVE_LINK_CACHE_EXPIRATION
  Cache expiration time for negative transcript<->protein links from the NCBI
  (in seconds).

  `Default value:` `60 * 60 * 24 * 30` (30 days)

USE_RELOADER
  Enable the `Werkzeug reloader
  <http://werkzeug.pocoo.org/docs/0.10/serving/#reloader>`_ for the website.

  This is disabled by default due to `a bug with using the reloader
  <https://github.com/mitsuhiko/werkzeug/issues/461#issuecomment-139369694>`_
  in combination with ``python -m mutalyzer.entrypoints.website``.

  `Default value:` `False`
