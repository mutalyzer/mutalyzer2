"""
Command line interface to Mutalyzer administrative tools.
"""


import argparse
import json
import os

import alembic.command
import alembic.config
from alembic.migration import MigrationContext
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm.exc import NoResultFound

from .. import announce
from .. import db
from ..db import session
from ..db.models import Assembly, Chromosome
from .. import mapping
from .. import output
from .. import sync


class UserError(Exception):
    pass


def add_assembly(assembly_file):
    """
    Add genome assembly definition from a JSON file.
    """
    try:
        definition = json.load(assembly_file)
    except ValueError:
        raise UserError('Not a valid JSON file: %s' % assembly_file.name)

    try:
        Assembly.by_name_or_alias(definition['name'])
    except NoResultFound:
        pass
    else:
        raise UserError('Assembly with this name or alias already exists: %s'
                        % definition['name'])

    if definition['alias'] is not None:
        try:
            Assembly.by_name_or_alias(definition['alias'])
        except NoResultFound:
            pass
        else:
            raise UserError('Assembly with this name or alias already '
                            'exists: %s' % definition['alias'])

    assembly = Assembly(definition['name'], definition['taxonomy_id'],
                        definition['taxonomy_common_name'],
                        definition['alias'])
    session.add(assembly)

    for chromosome_definition in definition['chromosomes']:
        chromosome = Chromosome(assembly, chromosome_definition['name'],
                                chromosome_definition['accession'],
                                chromosome_definition['organelle'])
        session.add(chromosome)

    session.commit()


def list_assemblies():
    """
    List genome assemblies.
    """
    assemblies = Assembly.query \
        .order_by(Assembly.taxonomy_common_name.asc(),
                  Assembly.name.asc()) \
        .all()

    for assembly in assemblies:
        if assembly.alias:
            name = '%s (%s)' % (assembly.name, assembly.alias)
        else:
            name = assembly.name
        print '%s, %s (%s)' % (name, assembly.taxonomy_common_name,
                               assembly.taxonomy_id)


def import_mapview(assembly_name_or_alias, mapview_file, group_label):
    """
    Import transcript mappings from an NCBI mapview file.
    """
    try:
        assembly = Assembly.by_name_or_alias(assembly_name_or_alias)
    except NoResultFound:
        raise UserError('Not a valid assembly: %s' % assembly_name_or_alias)

    try:
        mapping.import_from_mapview_file(assembly, mapview_file, group_label)
    except mapping.MapviewSortError as e:
        raise UserError(str(e))


def import_gene(assembly_name_or_alias, gene):
    """
    Import transcript mappings for a gene from the UCSC database.
    """
    try:
        assembly = Assembly.by_name_or_alias(assembly_name_or_alias)
    except NoResultFound:
        raise UserError('Not a valid assembly: %s' % assembly_name_or_alias)

    mapping.import_from_ucsc_by_gene(assembly, gene)


def import_reference(assembly_name_or_alias, reference):
    """
    Import transcript mappings from a genomic reference. Note that this is
    currently hard-coded to only work with mtDNA transcripts.
    """
    try:
        assembly = Assembly.by_name_or_alias(assembly_name_or_alias)
    except NoResultFound:
        raise UserError('Not a valid assembly: %s' % assembly_name_or_alias)

    mapping.import_from_reference(assembly, reference)


def sync_cache(wsdl_url, url_template, history=7):
    """
    Synchronize the database cache with another Mutalyzer instance.

    This program is intended to be run daily from cron. Example:

        25 5 * * *  mutalyzer-cache-sync 'http://dom1/?wsdl' 'http://dom1/{file}' -H 7
        55 5 * * *  mutalyzer-cache-sync 'http://dom2/?wsdl' 'http://dom2/{file}' -H 7
    """
    cache_sync = sync.CacheSync(output.Output(__file__))
    cache_sync.sync_with_remote(wsdl_url, url_template, history)


def set_announcement(body, url=None):
    """
    Set announcement to show to the user.
    """
    announce.set_announcement(body, url=url)


def unset_announcement():
    """
    Unset announcement to show to the user.
    """
    announce.unset_announcement()


def setup_database(alembic_config_path=None, destructive=False):
    """
    Setup database tables (if they do not yet exist).
    """
    if alembic_config_path and not os.path.isfile(alembic_config_path):
        raise UserError('Cannot find Alembic configuration: %s'
                        % alembic_config_path)

    if destructive:
        db.Base.metadata.drop_all(db.session.get_bind())

    db.Base.metadata.create_all(db.session.get_bind())
    db.session.commit()

    if alembic_config_path:
        context = MigrationContext.configure(db.session.connection())
        if destructive or context.get_current_revision() is None:
            # We need to close the current session before running Alembic.
            db.session.remove()
            alembic_config = alembic.config.Config(alembic_config_path)
            alembic.command.stamp(alembic_config, 'head')


def main():
    """
    Command-line interface to Mutalyzer administrative tools.
    """
    assembly_parser = argparse.ArgumentParser(add_help=False)
    assembly_parser.add_argument(
        '-a', '--assembly', metavar='ASSEMBLY', dest='assembly_name_or_alias',
        default='hg19', help='assembly to import to (default: hg19)')

    parser = argparse.ArgumentParser(
        description='Mutalyzer administrative tools.')
    subparsers = parser.add_subparsers(
        title='subcommands', dest='subcommand', help='subcommand help')

    # Subparsers for 'assemblies'.
    s = subparsers.add_parser(
        'assemblies', help='manage genome assemblies',
        description='Manage genome assemblies and their transcript mappings.'
        ).add_subparsers()

    # Subparser 'assemblies list'.
    p = s.add_parser(
        'list', help='list assemblies',
        description=list_assemblies.__doc__.split('\n\n')[0])
    p.set_defaults(func=list_assemblies)

    # Subparser 'assemblies add'.
    p = s.add_parser(
        'add', help='add assembly definition from JSON file',
        description=add_assembly.__doc__.split('\n\n')[0])
    p.set_defaults(func=add_assembly)
    p.add_argument(
        'assembly_file', metavar='FILE', type=argparse.FileType('r'),
        help='genome assembly definition JSON file (example: '
        'extras/assemblies/GRCh37.json)')

    # Subparser 'assemblies import-mapview'.
    p = s.add_parser(
        'import-mapview', help='import mappings from NCBI mapview file',
        parents=[assembly_parser],
        description=import_mapview.__doc__.split('\n\n')[0],
        epilog='Note: We require that FILE is sorted on the `feature_id` '
        '(#11) and `chromosome` (#2) columns. This can be done with a '
        '`sort -k 11,11 -k 2,2` command.')
    p.set_defaults(func=import_mapview)
    p.add_argument(
        'mapview_file', metavar='FILE', type=argparse.FileType('r'),
        help='file from NCBI mapview (example: seq_gene.md), see note below')
    p.add_argument(
        'group_label', metavar='GROUP_LABEL',
        help='use only entries with this group label (example: '
        'GRCh37.p2-Primary Assembly)')

    # Subparser 'assemblies import-gene'.
    p = s.add_parser(
        'import-gene', help='import mappings by gene from UCSC database',
        parents=[assembly_parser],
        description=import_gene.__doc__.split('\n\n')[0],
        epilog='Intended use is whenever transcript mappings for a specific '
        'gene are required that are not available from our usual source '
        ' (i.e., NCBI mapview).')
    p.set_defaults(func=import_gene)
    p.add_argument(
        'gene', metavar='GENE_SYMBOL',
        help='gene to import all transcript mappings for from the UCSC '
        'database (example: TTN)')

    # Subparser 'assemblies import-reference'.
    p = s.add_parser(
        'import-reference', help='import mappings from reference',
        parents=[assembly_parser],
        description=import_reference.__doc__.split('\n\n')[0],
        epilog='Intended use is whenever transcript mappings from a specific '
        'genomic reference are required that are not available from our '
        'usual source (i.e., NCBI mapview).')
    p.set_defaults(func=import_reference)
    p.add_argument(
        'reference', metavar='ACCESSION',
        help='genomic reference to import all genes from (example: '
        'NC_012920.1)')

    # Subparsers for 'announcement'.
    s = subparsers.add_parser(
        'announcement', help='manage user announcement',
        description='Manage announcement to show to the user.',
        epilog='The announcement is shown on every page of the website.'
        ).add_subparsers()

    # Subparser 'announcement set'.
    p = s.add_parser(
        'set', help='set user announcement',
        description=set_announcement.__doc__.split('\n\n')[0])
    p.set_defaults(func=set_announcement)
    p.add_argument(
        'body', metavar='ANNOUNCEMENT',
        help='announcement text to show to the user')
    p.add_argument(
        '--url', metavar='URL', dest='url',
        help='URL to more information on the announcement')

    # Subparser 'announcement unset'.
    p = s.add_parser(
        'unset', help='unset user announcement',
        description=unset_announcement.__doc__.split('\n\n')[0])
    p.set_defaults(func=unset_announcement)

    # Subparser 'sync-cache'.
    p = subparsers.add_parser(
        'sync-cache', help='synchronize cache with remote Mutalyzer',
        description=sync_cache.__doc__.split('\n\n')[0],
        epilog='Intended use is to run daily from cron.')
    p.add_argument(
        'wsdl_url', metavar='WSDL_URL',
        help='location of the remote WSDL description')
    p.add_argument(
        'url_template', metavar='URL_TEMPLATE',
        help='URL for remote downloads, in which the filename is to be '
        'substituted for {file}')
    p.add_argument(
        '-H', '--history', metavar='DAYS', dest='history', type=int,
        default=7, help='number of days to go back in the remote cache '
        '(default: 7)')

    # Subparser 'setup-database'.
    p = subparsers.add_parser(
        'setup-database', help='setup database',
        description=setup_database.__doc__.split('\n\n')[0],
        epilog='If Alembic config is given (--alembic-config), this also '
        'prepares the database for future migrations with Alembic '
        '(recommended).')
    p.add_argument(
        '--destructive', dest='destructive', action='store_true',
        help='delete any existing tables and data')
    p.add_argument(
        '-c', '--alembic-config', metavar='ALEMBIC_CONFIG',
        dest='alembic_config_path', help='path to Alembic configuration file')
    p.set_defaults(func=setup_database)

    args = parser.parse_args()

    try:
        args.func(**{k: v for k, v in vars(args).items()
                     if k not in ('func', 'subcommand')})
    except UserError as e:
        parser.error(str(e))


if __name__ == '__main__':
    main()
