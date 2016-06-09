"""Fix data migrations for Reference.source

Revision ID: 8ed4dce84496
Revises: 33c9465d5a60
Create Date: 2016-06-09 22:38:01.933037

"""

from __future__ import unicode_literals

# revision identifiers, used by Alembic.
revision = '8ed4dce84496'
down_revision = u'33c9465d5a60'

from alembic import op
import sqlalchemy as sa
from sqlalchemy import sql


def upgrade():
    # We repeat the data migration from migration c92f30c03b62 for sliced
    # references to correct the error in its original release. Slice start
    # position was used for both start and stop positions in the new column.
    connection = op.get_bind()

    # Inline table definition we can use in this migration.
    references = sql.table(
        'references',
        sql.column('id', sa.Integer()),
        sql.column('source', sa.Enum('ncbi', 'ncbi_slice', 'lrg', 'url', 'upload', name='reference_source')),
        sql.column('source_data', sa.String(255)),
        sql.column('slice_accession', sa.String(20)),
        sql.column('slice_start', sa.Integer()),
        sql.column('slice_stop', sa.Integer()),
        sql.column('slice_orientation', sa.Enum('forward', 'reverse', name='slice_orientation')))

    # Get all rows.
    result = connection.execute(
        references.select().with_only_columns([
            references.c.id,
            references.c.source,
            references.c.source_data,
            references.c.slice_accession,
            references.c.slice_start,
            references.c.slice_stop,
            references.c.slice_orientation
        ]).where(references.c.source == 'ncbi_slice'))

    # Generate parameter values for the UPDATE query below.
    def update_params(r):
        data = r.source_data
        if r.slice_accession:
            data = '{}:{}:{}:{}'.format(r.slice_accession, r.slice_start, r.slice_stop, r.slice_orientation)
        return {'r_id': r.id, 'r_source_data': data}

    # Process a few rows at a time, since they will be read in memory.
    while True:
        chunk = result.fetchmany(1000)
        if not chunk:
            break

        # Populate `source_data` based on existing column values.
        statement = references.update().where(
            references.c.id == sql.bindparam('r_id')
        ).values({'source_data': sql.bindparam('r_source_data')})

        # Execute UPDATE query for fetched rows.
        connection.execute(statement, [update_params(r) for r in chunk])


def downgrade():
    pass
