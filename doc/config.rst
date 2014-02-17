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

.. note:: Todo: This section is incomplete.

Note that the configuration file is interpreted as a Python module, so you can
use arbitrary Python expressions as configuration values, or even import other
modules in it.

Unsetting a configuration setting is done by using the value `None`. If no
default value is mentioned for any configuration setting below it means it is
not set by default.


Database settings
^^^^^^^^^^^^^^^^^

DATABASE_URI
  SQLAlchemy database connection URI specifying the database used to store
  users, samples, variants, etcetera.

  ================   =====================================================
  Database system    Example URI
  ================   =====================================================
  PostgreSQL         ``postgresql://mutalyzer:*****@localhost/mutalyzer``
  MySQL              ``mysql://mutalyzer:*****@localhost/mutalyzer``
  SQLite             ``sqlite:////tmp/mutalyzer.db``
  ================   =====================================================

  See the SQLAlchemy documentation on
  `Engine Configuration
  <http://docs.sqlalchemy.org/en/latest/core/engines.html>`_ for more
  information.

  `Default value:` ``sqlite://`` (in-memory SQLite database)


Miscellaneous settings
^^^^^^^^^^^^^^^^^^^^^^

TESTING
  If set to `True`, Mutalyzer assumes to be running its unit tests. This is
  done automatically in the provided test suite, so you should never have to
  change this setting.

  `Default value:` `False`
