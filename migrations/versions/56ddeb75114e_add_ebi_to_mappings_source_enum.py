"""Add EBI to mappings source enum

Revision ID: 56ddeb75114e
Revises: 1ed411f9fdfa
Create Date: 2016-02-11 04:21:40.658981

"""

from __future__ import unicode_literals

# revision identifiers, used by Alembic.
revision = '56ddeb75114e'
down_revision = u'1ed411f9fdfa'

from alembic import op
import sqlalchemy as sa


def upgrade():
    context = op.get_context()

    if context.bind.dialect.name == 'postgresql':
        # In PostgreSQL < 9.1 there was no ALTER TYPE for enums, so it would
        # have been something like:
        #
        #     ALTER TABLE foo ALTER COLUMN bar TYPE new_type USING bar::text::new_type;
        #
        # However, all my installations are PostgreSQL >= 9.1 and I think the
        # USING syntax is PostgreSQL-specific, so let's ignore that. It would
        # also come with all the hassle of moving old column values into the
        # new column.
        if context.bind.dialect.server_version_info >= (9, 3):
            op.execute('COMMIT')
            op.execute("ALTER TYPE source ADD VALUE IF NOT EXISTS 'ebi'")
            return
        if context.bind.dialect.server_version_info >= (9, 1):
            op.execute('COMMIT')
            op.execute("ALTER TYPE source ADD VALUE 'ebi'")
            return

    elif context.bind.dialect.name == 'sqlite':
        # SQLite doesn't support altering columns, so we have to wrap this in
        # a batch operation.
        with op.batch_alter_table('transcript_mappings') as batch_op:
            batch_op.alter_column(
                'source', nullable=False,
                type_=sa.Enum('ucsc', 'ncbi', 'ebi', 'reference', name='source')
            )
        return

    elif context.bind.dialect.name == 'mysql':
        # In MySQL we can simply alter the column.
        op.alter_column(
            'transcript_mappings', 'source', nullable=False,
            type_=sa.Enum('ucsc', 'ncbi', 'ebi', 'reference', name='source'),
            existing_type=sa.Enum('ucsc', 'ncbi', 'reference', name='source')
        )
        return

    raise Exception('Sorry, only PostgreSQL >= 9.1, SQLite, and MySQL are supported by this migration')


def downgrade():
    raise Exception('Downgrade not supported by this migration')
