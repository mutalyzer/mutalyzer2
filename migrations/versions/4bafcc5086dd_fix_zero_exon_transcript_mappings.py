"""Fix zero-exon transcript mappings

Revision ID: 4bafcc5086dd
Revises: 2e062969eb54
Create Date: 2015-07-20 16:16:01.602964

"""

from __future__ import unicode_literals

# revision identifiers, used by Alembic.
revision = '4bafcc5086dd'
down_revision = u'2e062969eb54'

from alembic import op
from sqlalchemy import sql
import sqlalchemy as sa


def upgrade():
    transcript_mappings = sql.table('transcript_mappings',
                                    sql.column('start', sa.Integer()),
                                    sql.column('stop', sa.Integer()),
                                    sql.column('exon_starts', sa.Text()),
                                    sql.column('exon_stops', sa.Text()))
    # https://alembic.readthedocs.org/en/latest/ops.html#alembic.operations.Operations.execute
    op.execute(transcript_mappings
               .update()
               .where(transcript_mappings.c.exon_starts == op.inline_literal(''))
               .values({'exon_starts': transcript_mappings.c.start,
                        'exon_stops': transcript_mappings.c.stop}))


def downgrade():
    # We cannot reliably downgrade this migration.
    pass
