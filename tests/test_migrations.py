"""
Test database migrations.
"""


from __future__ import unicode_literals

from datetime import datetime

import alembic.autogenerate
import alembic.command
import alembic.config
from alembic.migration import MigrationContext
import sqlalchemy as sa
from sqlalchemy import create_engine, sql

from mutalyzer import db


def test_migrations(database_uri):
    """
    Run all migrations and assert the result is up to date with the model
    definitions.
    """
    alembic_config = alembic.config.Config('migrations/alembic.ini')
    engine = create_engine(database_uri)

    with engine.begin() as connection:
        # http://alembic.readthedocs.org/en/latest/cookbook.html#sharing-a-connection-with-a-series-of-migration-commands-and-environments
        alembic_config.attributes['connection'] = connection

        if database_uri != 'sqlite://':
            db.Base.metadata.drop_all(connection)

        # Create initial schema by running the first migration.
        alembic.command.upgrade(alembic_config, 'ea660b66f26')

        # Add some database content to run the migrations on.
        add_database_content(connection)

        # Run the remaining migrations.
        alembic.command.upgrade(alembic_config, 'head')

        context = MigrationContext.configure(connection)
        assert not alembic.autogenerate.compare_metadata(
            context, db.Base.metadata)

    engine.dispose()


def add_database_content(connection):
    """
    Add some content to the database.
    """
    # We only define tables and columns we actually need, so this is not a
    # complete mapping of the schema.

    references = sql.table(
        'references',
        sql.column('id', sa.Integer),
        sql.column('accession', sa.String(20)),
        sql.column('checksum', sa.String(32)),
        sql.column('geninfo_identifier', sa.String(13)),
        sql.column('slice_accession', sa.String(20)),
        sql.column('slice_start', sa.Integer),
        sql.column('slice_stop', sa.Integer),
        sql.column('slice_orientation', sa.Enum('forward', 'reverse', name='slice_orentation')),
        sql.column('download_url', sa.String(255)),
        sql.column('added', sa.DateTime))

    assemblies = sql.table(
        'assemblies',
        sql.column('id', sa.Integer),
        sql.column('name', sa.String(30)),
        sql.column('alias', sa.String(10)),
        sql.column('taxonomy_id', sa.Integer),
        sql.column('taxonomy_common_name', sa.String(50)))

    chromosomes = sql.table(
        'chromosomes',
        sql.column('id', sa.Integer),
        sql.column('assembly_id', sa.Integer),
        sql.column('name', sa.String(30)),
        sql.column('accession', sa.String(30)),
        sql.column('organelle', sa.Enum('nucleus', 'mitochondrion',
                                        name='organelle')))

    transcript_mappings = sql.table(
        'transcript_mappings',
        sql.column('chromosome_id', sa.Integer),
        sql.column('reference_type', sa.Enum('refseq', 'lrg',
                                             name='reference_type')),
        sql.column('accession', sa.String(20)),
        sql.column('gene', sa.String(30)),
        sql.column('transcript', sa.Integer),
        sql.column('orientation', sa.Enum('forward', 'reverse',
                                          name='orentation')),
        sql.column('start', sa.Integer),
        sql.column('stop', sa.Integer),
        sql.column('exon_starts', sa.Text),
        sql.column('exon_stops', sa.Text),
        sql.column('select_transcript', sa.Boolean),
        sql.column('source', sa.Enum('ucsc', 'ncbi', 'reference',
                                     name='source')))

    transcript_protein_links = sql.table(
        'transcript_protein_links',
        sql.column('transcript_accession', sa.String(30)),
        sql.column('protein_accession', sa.String(30)),
        sql.column('added', sa.DateTime))

    # Add some common data.
    connection.execute(
        assemblies.insert(),
        name='GRCh37',
        taxonomy_id=9606,
        taxonomy_common_name='Homo sapiens',
        alias='hg19')
    hg19_id = connection.execute(
        assemblies.select(assemblies.c.alias == 'hg19')
        .with_only_columns([assemblies.c.id])
    ).fetchone()[0]

    connection.execute(
        chromosomes.insert(),
        assembly_id=hg19_id,
        name='chr1',
        accession='NC_000001.10',
        organelle='nucleus')
    chr1_id = connection.execute(
        chromosomes.select(chromosomes.c.name == 'chr1')
        .with_only_columns([chromosomes.c.id])
    ).fetchone()[0]

    # Data for migration 402ff01b0d5d:
    # Fix GRCm38 chromosome accession number versions.
    connection.execute(
        chromosomes.insert(),
        assembly_id=hg19_id,
        name='chr11',
        accession='NC_000077.60',
        organelle='nucleus')

    # Data for migration 2e062969eb54:
    # Rename GRCh36 assembly to NCBI36.
    connection.execute(
        assemblies.insert(),
        name='GRCh36',
        taxonomy_id=9606,
        taxonomy_common_name='Homo sapiens',
        alias='hg18')

    # Data for migration 4bafcc5086dd:
    # Fix zero-exon transcript mappings.
    connection.execute(
        transcript_mappings.insert(),
        chromosome_id=chr1_id,
        reference_type='refseq',
        accession='NC_001807',
        gene='ATP6',
        transcript=1,
        orientation='forward',
        start=8528,
        stop=9208,
        exon_starts='8528',
        exon_stops='9208',
        select_transcript=True,
        source='ncbi')
    connection.execute(
        transcript_mappings.insert(),
        chromosome_id=chr1_id,
        reference_type='refseq',
        accession='NC_001807',
        gene='ATP8',
        transcript=1,
        orientation='forward',
        start=8367,
        stop=8573,
        exon_starts='',
        exon_stops='',
        select_transcript=True,
        source='ncbi')

    # Data for migration 3492d2ee8884:
    # Transcript protein links have nullable transcript and unique protein.
    connection.execute(
        transcript_protein_links.insert(),
        transcript_accession='NM_052818',
        protein_accession='NP_438169',
        added=datetime.now())
    connection.execute(
        transcript_protein_links.insert(),
        transcript_accession='NM_001079691',
        protein_accession=None,
        added=datetime.now())

    # Data for migration c92f30c03b62:
    # Add source and source_data columns to references.
    connection.execute(
        references.insert(),
        accession='UD_146488222017',
        checksum='d0a72a13ac975e320f7f0ccee19e89ab',
        geninfo_identifier=None,
        slice_accession=None,
        slice_start=None,
        slice_stop=None,
        slice_orientation=None,
        download_url=None,
        added=datetime.now())
    connection.execute(
        references.insert(),
        accession='UD_146488248189',
        checksum='4742c1ec7cce47fb6ef2cf0475b0bfd1',
        geninfo_identifier=None,
        slice_accession=None,
        slice_start=None,
        slice_stop=None,
        slice_orientation=None,
        download_url='https://barmsijs.lumc.nl/HG/mvermaat/UD_132118785483.gb',
        added=datetime.now())
    connection.execute(
        references.insert(),
        accession='UD_146488309212',
        checksum='d57af526097b5d353792b1ef46b08fa0',
        geninfo_identifier=None,
        slice_accession='NC_000023.11',
        slice_start=31117219,
        slice_stop=33344609,
        slice_orientation='reverse',
        download_url=None,
        added=datetime.now())
    connection.execute(
        references.insert(),
        accession='UD_146488315868',
        checksum='e8a906c451e4228e636c24cdc3140429',
        geninfo_identifier=None,
        slice_accession='NC_000003.11',
        slice_start=3000000,
        slice_stop=3005000,
        slice_orientation='forward',
        download_url=None,
        added=datetime.now())
    connection.execute(
        references.insert(),
        accession='UD_146488317444',
        checksum='4002698fe1ae3f2960cbbe840983b330',
        geninfo_identifier=None,
        slice_accession='NC_000007.13',
        slice_start=5000000,
        slice_stop=5005000,
        slice_orientation='reverse',
        download_url=None,
        added=datetime.now())
    connection.execute(
        references.insert(),
        accession='AB026906.1',
        checksum='29b003d5a71af74dc61a92d2ef5cd5d9',
        geninfo_identifier='5295993',
        slice_accession=None,
        slice_start=None,
        slice_stop=None,
        slice_orientation=None,
        download_url=None,
        added=datetime.now())
    connection.execute(
        references.insert(),
        accession='NG_012337.1',
        checksum='8f3dd98b84999a2449dceb7ee6236313',
        geninfo_identifier='254039638',
        slice_accession=None,
        slice_start=None,
        slice_stop=None,
        slice_orientation=None,
        download_url=None,
        added=datetime.now())
    connection.execute(
        references.insert(),
        accession='LRG_9',
        checksum='30cca9e429227628f365598eab24407e',
        geninfo_identifier=None,
        slice_accession=None,
        slice_start=None,
        slice_stop=None,
        slice_orientation=None,
        download_url='ftp://ftp.ebi.ac.uk/pub/databases/lrgex/SCHEMA_1_7_ARCHIVE/pending/LRG_9.xml',
        added=datetime.now())
