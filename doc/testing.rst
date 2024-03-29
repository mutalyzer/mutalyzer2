.. highlight:: none

.. _testing:

Testing
=======


Unit tests
----------

We use `pytest`_ for the unit tests. To run them, just type ``py.test`` from
the Mutalyzer source directory.

.. note:: The Mutalyzer package must be installed before running the unit
          tests.

By default, the tests use an in-memory SQLite database. This can be customized
with the ``--database-uri`` command line argument and a valid `SQLAlchemy
connection URI
<http://docs.sqlalchemy.org/en/rel_1_0/core/engines.html#database-urls>`_
(obviously, the contents of this database will be lost). For example, to use
an SQLite database on the filesystem::

    $ py.test --database-uri sqlite:////tmp/mutalyzer.sql

Or, using `pg_virtualenv
<https://alioth.debian.org/scm/loggerhead/pkg-postgresql/postgresql-common/trunk/view/head:/pg_virtualenv>`_
(included with the Debian PostgreSQL packages), to run the tests with
PostgreSQL::

    $ pg_virtualenv bash -c 'py.test --database-uri postgres://${PGUSER}:${PGPASSWORD}@${PGHOST}:${PGPORT}/${PGDATABASE}'

Multiple ``--database-uri`` arguments are allowed. Tests using the database
will be run once for every database specified.

Similarly, ``--redis-uri`` (only one allowed) specifies a Redis server to use
for testing. If unspecified, a mock Redis server is used.

Tests are `run automatically on Travis CI
<https://travis-ci.org/mutalyzer/mutalyzer2>`_ with SQLite, PostgreSQL, and
MySQL, for each pull request and push on GitHub.


Testing the web services
------------------------

To ease testing the web services during development, some simple web service
client scripts are included in the Mutalyzer source tree::

    $ cd extras/soap-tools
    $ ./info.py
    Version: 2.0.31
    Version parts: 2, 0, 31
    Release date: 21 August 2019
    Nomenclature version: 2.0
    Nomenclature version parts: 2, 0
    Server name: res-muta-app01
    Contact e-mail: info@mutalyzer.nl
    $

They simply call one of the web service functions and print the result. You
may have to change the server location defined at the top of these scripts.

.. note:: One of the scripts, ``run_batch_job.py``, provides an easy way to
          run a batch job from the command line. Some additional notes are
          available for `running this on a Windows machine
          <https://gist.github.com/jfjlaros/482fe9f0397e554ed29f>`_.


.. _pytest: http://pytest.org/
