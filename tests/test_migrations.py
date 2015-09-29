"""
Test database migrations.
"""


from __future__ import unicode_literals

import os

import alembic.autogenerate
import alembic.command
import alembic.config
from alembic.migration import MigrationContext
from sqlalchemy import create_engine

from mutalyzer import db


def test_migrations():
    """
    Run all migrations and assert the result is up to date with the model
    definitions.

    We don't use `utils.MutalyzerTest` here, or `mutalyzer.db.session` in any
    way for that matter, since it will bootstrap the database schema.
    """
    database_uri = os.getenv('MUTALYZER_TEST_DATABASE_URI', 'sqlite://')

    alembic_config = alembic.config.Config('migrations/alembic.ini')
    engine = create_engine(database_uri)

    with engine.begin() as connection:
        # http://alembic.readthedocs.org/en/latest/cookbook.html#sharing-a-connection-with-a-series-of-migration-commands-and-environments
        alembic_config.attributes['connection'] = connection

        if database_uri != 'sqlite://':
            db.Base.metadata.drop_all(connection)

        # TODO: We might also want to test this with data. An implementation
        #   might look like this:
        #
        #   1. Create initial schema by running the first migration.
        #   2. Add some database content.
        #   3. Run the remaining migrations.

        alembic.command.upgrade(alembic_config, 'head')

        context = MigrationContext.configure(connection)
        assert not alembic.autogenerate.compare_metadata(
            context, db.Base.metadata)

    engine.dispose()
