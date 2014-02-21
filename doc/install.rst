.. highlight:: none

.. _install:

Installation
============

Mutalyzer depends on a database server, `Python`_ 2.7, and several Python
packages. `Redis`_ is a soft dependency. This section walks you through
installing Mutalyzer with Redis and using `PostgreSQL`_ as database server,
which is the recommended setup.

.. note:: All operating system specific instructions assume installation on a
   `Debian`_ 7 *wheezy* system. You'll have to figure out the necessary
   adjustements yourself if you're on another system.

The following steps will get Mutalyzer running on your system with the
recommended setup:

* :ref:`install-postgresql`
* :ref:`install-redis`
* :ref:`install-virtualenv`
* :ref:`install-setup`

At the bottom of this page some :ref:`alternative setups
<install-alternatives>` are documented.


.. _install-quick:

If you're in a hurry
--------------------

The impatient can run Mutalyzer without a database server and more such
nonsense with the following steps::

    $ pip install -r requirements.txt
    $ MUTALYZER_SETTINGS=/dev/null python -m mutalyzer.entrypoints.website

This starts the website frontend on the reported port using an in-memory
SQLite database.


.. _install-postgresql:

Database server: PostgreSQL
---------------------------

Install `PostgreSQL`_ and add a user for Mutalyzer. Create a database
(e.g. ``mutalyzer``) owned by the new user. For example::

    $ sudo apt-get install postgresql
    $ sudo -u postgres createuser --superuser $USER
    $ createuser --pwprompt --encrypted --no-adduser --no-createdb --no-createrole mutalyzer
    $ createdb --encoding=UNICODE --owner=mutalyzer mutalyzer

Also install some development libraries needed for building the ``psycopg2``
Python package later and add the package to the list of requirements::

    $ sudo apt-get install libpq-dev
    $ echo psycopg2 >> requirements.txt

This will make sure the Python PostgreSQL database adapter gets installed in
the :ref:`install-virtualenv` section.

.. seealso::

   :ref:`install-mysql`
     Alternatively, MySQL can be used as database server.

   :ref:`install-sqlite`
     Alternatively, SQLite can be used as database server.

   `Dialects -- SQLAlchemy documentation <http://docs.sqlalchemy.org/en/latest/dialects/index.html>`_
     In theory, any database supported by SQLAlchemy could work.


.. _install-redis:

Redis
-----

Mutalyzer uses Redis for non-critical fast storage such as statistics::

    $ sudo apt-get install redis-server

.. note:: Redis is a soft dependency, meaning that Mutalyzer will run without
    it (but may lack some non-essential features).


.. _install-virtualenv:

Python virtual environment
--------------------------

It is recommended to run Mutalyzer from a Python virtual environment, using
`virtualenv`_. Installing virtualenv and creating virtual environments is not
covered here.

Assuming you created and activated a virtual environment for Mutalyzer,
install all required Python packages::

    $ sudo apt-get install python-dev libmysqlclient-dev libxml2-dev libxslt-dev swig
    $ pip install -r requirements.txt

Install Mutalyzer::

    $ python setup.py install

.. note:: If you're planning on modifying the Mutalyzer source code, it might
    be convenient to install Mutalyzer in *development mode*::

        $ python setup.py develop

    Instead of copying the source code to the installation directory, this
    only links from the installation directory to the source code such that
    any changes you make to it are directly available in the environment.

Now might be a good time to run the unit tests::

    $ py.test

.. seealso::

   `virtualenv`_
     ``virtualenv`` is a tool to create isolated Python environments.

   `virtualenvwrapper`_
     ``virtualenvwrapper`` is a set of extensions to the ``virtualenv``
     tool. The extensions include wrappers for creating and deleting virtual
     environments and otherwise managing your development workflow.


.. _install-setup:

Mutalyzer setup
---------------

Mutalyzer looks for its configuration in the file specified by the
``MUTALYZER_SETTINGS`` environment variable. First create the file with your
configuration settings, for example::

    $ export MUTALYZER_SETTINGS=~/mutalyzer/settings.py
    $ cat > $MUTALYZER_SETTINGS
    REDIS_URI = 'redis://localhost'
    DATABASE_URI = 'postgresql://mutalyzer:*****@localhost/mutalyzer'

A script is included to setup the database::

    $ mutalyzer-admin setup-database --alembic-config migrations/alembic.ini

You can now proceed to :ref:`run`.

.. seealso::

   :ref:`config`
     For more information on the available configuration settings.


.. _install-alternatives:

Alternative setups
------------------

The remainder of this page documents some alternatives to the recommended
setup documented above.


.. _install-mysql:

Database server: MySQL
^^^^^^^^^^^^^^^^^^^^^^

Install `MySQL`_ and create a database (e.g. ``mutalyzer``) with all privileges
for the Mutalyzer user. For example::

    $ sudo apt-get install mysql-server
    $ mysql -h localhost -u root -p
    > create database mutalyzer;
    > grant all privileges on mutalyzer.* to mutalyzer@localhost identified by '*****';

The Python MySQL database adapter is a hard dependency regardless of your
choice of database server, so it'll get installed in the
:ref:`install-virtualenv` section.

In the :ref:`install-setup` section, make sure to use a MySQL database URI in
the Mutalyzer settings file, e.g.:

.. code-block:: python

    DATABASE_URI = 'mysql://mutalyzer:*****@localhost/mutalyzer'

.. seealso::

   :ref:`install-postgresql`
     The recommended setup uses PostgreSQL as database server.


.. _install-sqlite:

Database server: SQLite
^^^^^^^^^^^^^^^^^^^^^^^

You probably already have all you need for using `SQLite`_, so this section
consists of zero steps.

Just note that in the :ref:`install-setup` section, you should use an SQLite
database URI in the Mutalyzer settings file, e.g.:

.. code-block:: python

    DATABASE_URI = 'sqlite:////tmp/mutalyzer.db'

.. seealso::

   :ref:`install-postgresql`
     The recommended setup uses PostgreSQL as database server.


.. _Debian: http://www.debian.org/
.. _MySQL: http://www.mysql.com/
.. _PostgreSQL: http://www.postgresql.org/
.. _Python: http://python.org/
.. _Redis: http://redis.io/
.. _SQLite: http://www.sqlite.org/
.. _virtualenv: http://www.virtualenv.org/
.. _virtualenvwrapper: http://www.doughellmann.com/docs/virtualenvwrapper/
