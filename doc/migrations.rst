.. highlight:: none

.. _migrations:

Database migrations
===================

Managing database migrations is done using `Alembic`_, a database migration tool
for usage with SQLAlchemy.


Supported database systems
--------------------------

We try to take care to fully support migrations on PostgreSQL, MySQL, and
SQLite. Since our production deployments run PostgreSQL, migrations may be
optimized only for that system.


Zero-downtime migrations
------------------------

To support zero-downtime database migrations, all migrations must keep
compatibility with the existing codebase. The assumption is that deployments
always run database migrations first and update the code second.

This means some database changes may have to be broken down into several steps
and completed over several deployments. For example, deleting a column from a
table can be done as follows:

1. Remove the code that uses the column.
2. Deploy:

   1. Run migrations (nothing to be done).
   2. Update the code.

3. Create a migration that removes the column.
4. Deploy:

   1. Run migrations (column is removed).
   2. Update the code (nothing to be done).

On the other hand, adding a column is easier. Creating a migration that adds
the column and adding code to use it can safely be done at the same time if we
assume migrations are always run before updating the code.

To support smooth deployments, ideally a release is made between migrations
and code updates that should be deployed separately. This way also external
users can always safely upgrade one release at a time.

.. note:: Our automated `deployment of Mutalyzer with Ansible
          <https://github.com/mutalyzer/ansible-role-mutalyzer>`_ runs
          database migrations before updating the code.


Setting up Alembic
------------------

Before you can use Alembic, the database needs to be stamped with the current
revision. This is done automatically when using the ``-c`` argument with the
``mutalyzer-admin setup-database`` command (as is recommended in
:ref:`install-setup`)::

    $ mutalyzer-admin setup-database --alembic-config migrations/alembic.ini

An existing database can also be stamped manually using Alembic::

    $ alembic -c migrations/alembic.ini stamp head


Running migrations
------------------

Upgrading an existing database to the latest revision is done as follows::

    $ alembic -c migrations/alembic.ini upgrade head

Downgrades are explicitely unsupported and some migrations may not have
downgrades implemented at all.


Creating migrations
-------------------

To create a new migration, first update the SQLAlchemy models in the Mutalyzer
source code. Then, generate a migration with Alembic::

    $ alembic -c migrations/alembic.ini revision --autogenerate -m 'Some descriptive comment'

Template code for the upgrade and downgrade paths are written to a new file in
the ``migrations/versions`` directory. Alembic is smart enough to generate the
complete code for simple migrations, but please always have a look at the
generated code.

For more complex changes to the database schema, you'll have to add some code
manually. The same goes for any data migrations you might want to
include. Consult the Alembic documentation and existing migrations for some
common patterns.


.. _Alembic: http://alembic.readthedocs.org/
