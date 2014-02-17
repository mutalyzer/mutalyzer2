.. highlight:: none

.. _upgrade:

Upgrading
=========

Before upgrading Mutalyzer, stop any currently running instances. Then, update
your copy of the source code (using for example ``git pull`` on an existing
git clone).

Make sure to install any new requirements::

    $ pip install -r requirements.txt

Now install the new version::

    $ python setup.py install

Managing database migrations is done using `Alembic`_. This command will move
your database to the latest schema::

    $ alembic -c migrations/alembic.ini upgrade head


.. _Alembic: http://alembic.readthedocs.org/
