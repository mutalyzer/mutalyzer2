"""
Command line interface to Mutalyzer administrative tools.
"""


from __future__ import unicode_literals

import argparse
import codecs
import json
import locale
import os

import alembic.command
import alembic.config
from alembic.migration import MigrationContext
import sqlalchemy
from sqlalchemy.orm.exc import NoResultFound

from . import _cli_string
from .. import announce
from ..redisclient import client
from .. import db
from ..db import session
from ..db.models import Assembly, BatchJob, BatchQueueItem, Chromosome
from .. import mapping
from .. import output
from .. import sync
from .. import util


class UserError(Exception):
    pass


def add_assembly(assembly_file, encoding):
    """
    Add genome assembly definition from a JSON file.
    """
    assembly_file = codecs.getreader(encoding)(assembly_file)

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
        .order_by(*Assembly.order_by_criteria) \
        .all()

    for assembly in assemblies:
        if assembly.alias:
            name = '%s (%s)' % (assembly.name, assembly.alias)
        else:
            name = assembly.name
        print '%s, %s (%s)' % (name, assembly.taxonomy_common_name,
                               assembly.taxonomy_id)


def import_mapview(assembly_name_or_alias, mapview_file, encoding,
                   group_label):
    """
    Import transcript mappings from an NCBI mapview file.
    """
    # For long-running processes it can be convenient to have a short and
    # human-readable process name.
    util.set_process_name('mutalyzer: mapview-import')

    mapview_file = codecs.getreader(encoding)(mapview_file)

    try:
        assembly = Assembly.by_name_or_alias(assembly_name_or_alias)
    except NoResultFound:
        raise UserError('Not a valid assembly: %s' % assembly_name_or_alias)

    try:
        mapping.import_from_mapview_file(assembly, mapview_file, group_label)
    except mapping.MapviewSortError as e:
        raise UserError(unicode(e))


def import_lrgmap(assembly_name_or_alias, lrgmap_file, encoding):
    """
    Import transcript mappings from an EBI LRG transcripts map file.
    """
    lrgmap_file = codecs.getreader(encoding)(lrgmap_file)

    try:
        assembly = Assembly.by_name_or_alias(assembly_name_or_alias)
    except NoResultFound:
        raise UserError('Not a valid assembly: %s' % assembly_name_or_alias)

    mapping.import_from_lrgmap_file(assembly, lrgmap_file)


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

    .. code-block:: none

        25 5 * * *  mutalyzer-cache-sync 'http://dom1/?wsdl' 'http://dom1/{file}' -H 7
        55 5 * * *  mutalyzer-cache-sync 'http://dom2/?wsdl' 'http://dom2/{file}' -H 7
    """
    # For long-running processes it can be convenient to have a short and
    # human-readable process name.
    util.set_process_name('mutalyzer: cache-sync')

    cache_sync = sync.CacheSync(output.Output(__file__))
    inserted, downloaded = cache_sync.sync_with_remote(wsdl_url, url_template,
                                                       history)

    print ('Added %d entries to the local cache and downloaded %d cache files.'
           % (inserted, downloaded))


def list_batch_jobs():
    """
    List batch jobs.
    """
    # For getting all batch jobs and their item counts, the following query
    # might be more obvious at first thought. However, our current query below
    # turns out to be more than twice as fast (and shorter).
    #
    #     sq = session.query(
    #         BatchQueueItem.batch_job_id,
    #         sqlalchemy.func.count(BatchQueueItem.id).label('count')
    #     ).group_by(BatchQueueItem.batch_job_id).subquery()
    #     session.query(
    #         BatchJob,
    #         sq.c.count
    #     ).join(sq, BatchJob.id == sq.c.batch_job_id)
    #
    batch_jobs_with_counts = session.query(
        BatchJob,
        session.query(sqlalchemy.func.count('*')).filter(
            BatchQueueItem.batch_job_id == BatchJob.id
        ).label('count')
    ).order_by(BatchJob.added.asc()).all()

    if len(batch_jobs_with_counts) < 1:
        return

    lengths = {
        'id_len': max(len(str(j.id)) for j, _ in batch_jobs_with_counts),
        'type_len': max(len(j.job_type) for j, _ in batch_jobs_with_counts),
        'count_len': max(len(str(c)) for _, c in batch_jobs_with_counts),
        'email_len': max(len(j.email) for j, _ in batch_jobs_with_counts)
    }

    template = ('{id:>{id_len}}  {type:{type_len}}  {added:%Y-%m-%d %H:%M:%S}'
                '  {count:>{count_len}}  {email:{email_len}}')

    for batch_job, count in batch_jobs_with_counts:
        print template.format(
            id=batch_job.id,
            type=batch_job.job_type,
            added=batch_job.added,
            count=count,
            email=batch_job.email,
            **lengths)


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

    bind = db.session.get_bind()

    if destructive:
        db.Base.metadata.drop_all(bind)

    if destructive or not bind.has_table(Assembly.__tablename__):
        # We assume our migrations will take care of everything if at least
        # the Assembly table exists.
        db.Base.metadata.create_all(bind)

    if alembic_config_path:
        context = MigrationContext.configure(db.session.connection())
        if destructive or context.get_current_revision() is None:
            # We need to close the current session before running Alembic.
            db.session.remove()
            alembic_config = alembic.config.Config(alembic_config_path)
            alembic.command.stamp(alembic_config, 'head')


def add_email_to_blacklist(email):
    """
    Add an e-mail address to the blacklist.
    """
    client.sadd('set:batch-job/website:email_blacklist', email)


def remove_email_from_blacklist(email):
    """
    Remove an e-mail address from the blacklist.
    """
    client.srem('set:batch-job/website:email_blacklist', email)


def list_emails_in_blacklist():
    """
    List all e-mail addresses in the blacklist.
    """
    for email in client.smembers('set:batch-job/website:email_blacklist'):
        print email


def main():
    """
    Command-line interface to Mutalyzer administrative tools.
    """
    default_encoding = locale.getpreferredencoding()

    assembly_parser = argparse.ArgumentParser(add_help=False)
    assembly_parser.add_argument(
        '-a', '--assembly', metavar='ASSEMBLY', type=_cli_string,
        dest='assembly_name_or_alias', default='hg19',
        help='assembly to import to (default: hg19)')

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
        'assembly_file', metavar='FILE', type=argparse.FileType('rb'),
        help='genome assembly definition JSON file (example: '
        'extras/assemblies/GRCh37.json)')
    p.add_argument(
        '--encoding', metavar='ENCODING', type=_cli_string,
        default=default_encoding,
        help='input file encoding (default: %s)' % default_encoding)

    # Subparser 'assemblies import-mapview'.
    p = s.add_parser(
        'import-mapview', help='import mappings from NCBI mapview file',
        parents=[assembly_parser],
        description=import_mapview.__doc__.split('\n\n')[0],
        epilog='Note: We require that FILE is sorted on the `feature_id` '
        '(#11) and `chromosome` (#2) columns. This can be done with a '
        '`sort -t $\'\\t\' -k 11,11 -k 2,2` command.')
    p.set_defaults(func=import_mapview)
    p.add_argument(
        'mapview_file', metavar='FILE', type=argparse.FileType('rb'),
        help='file from NCBI mapview (example: seq_gene.md), see note below')
    p.add_argument(
        '--encoding', metavar='ENCODING', type=_cli_string,
        default=default_encoding,
        help='input file encoding (default: %s)' % default_encoding)
    p.add_argument(
        'group_label', metavar='GROUP_LABEL', type=_cli_string,
        help='use only entries with this group label (example: '
        'GRCh37.p2-Primary Assembly)')

    # Subparser 'assemblies import-lrgmap'.
    p = s.add_parser(
        'import-lrgmap',
        help='import mappings from EBI LRG transcripts map file',
        parents=[assembly_parser],
        description=import_lrgmap.__doc__.split('\n\n')[0])
    p.set_defaults(func=import_lrgmap)
    p.add_argument(
        'lrgmap_file', metavar='FILE', type=argparse.FileType('rb'),
        help='EBI LRG transcript map file (example: '
        'list_LRGs_transcripts_GRCh37.txt)')
    p.add_argument(
        '--encoding', metavar='ENCODING', type=_cli_string,
        default=default_encoding,
        help='input file encoding (default: %s)' % default_encoding)

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
        'gene', metavar='GENE_SYMBOL', type=_cli_string,
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
        'reference', metavar='ACCESSION', type=_cli_string,
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
        'body', metavar='ANNOUNCEMENT', type=_cli_string,
        help='announcement text to show to the user')
    p.add_argument(
        '--url', metavar='URL', dest='url', type=_cli_string,
        help='URL to more information on the announcement')

    # Subparser 'announcement unset'.
    p = s.add_parser(
        'unset', help='unset user announcement',
        description=unset_announcement.__doc__.split('\n\n')[0])
    p.set_defaults(func=unset_announcement)

    # Subparser 'batch-jobs'.
    p = subparsers.add_parser(
        'batch-jobs', help='list batch jobs',
        description=list_batch_jobs.__doc__.split('\n\n')[0])
    p.set_defaults(func=list_batch_jobs)

    # Subparser 'sync-cache'.
    p = subparsers.add_parser(
        'sync-cache', help='synchronize cache with remote Mutalyzer',
        description=sync_cache.__doc__.split('\n\n')[0],
        epilog='Intended use is to run daily from cron.')
    p.add_argument(
        'wsdl_url', metavar='WSDL_URL', type=_cli_string,
        help='location of the remote WSDL description')
    p.add_argument(
        'url_template', metavar='URL_TEMPLATE', type=_cli_string,
        help='URL for remote downloads, in which the filename is to be '
        'substituted for {file}')
    p.add_argument(
        '-H', '--history', metavar='DAYS', dest='history', type=int,
        default=7, help='number of days to go back in the remote cache '
        '(default: 7)')
    p.set_defaults(func=sync_cache)

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
        '-c', '--alembic-config', metavar='ALEMBIC_CONFIG', type=_cli_string,
        dest='alembic_config_path', help='path to Alembic configuration file')
    p.set_defaults(func=setup_database)

    # Subparsers for 'email_blackist'.
    s = subparsers.add_parser(
        'email-blackist', help='manage e-mail blacklist',
        description='Manage blacklisted e-mails.',
        epilog='Blacklisting will notify the user that an e-mail address is '
            'invalid. Using an invalid e-mail address causes mail delivery '
            'problems (of which mutalyzer@humgen.nl gets notifications).'
        ).add_subparsers()

    # Subparser 'email_blackist add'.
    p = s.add_parser(
        'add', help='add an e-mail address to the blacklist',
        description=add_email_to_blacklist.__doc__.split('\n\n')[0])
    p.add_argument(
        'email', metavar='EMAIL', type=_cli_string, help='e-mail address')
    p.set_defaults(func=add_email_to_blacklist)

    # Subparser 'email_blackist del'.
    p = s.add_parser(
        'del', help='remove an e-mail address from the blacklist',
        description=remove_email_from_blacklist.__doc__.split('\n\n')[0])
    p.add_argument(
        'email', metavar='EMAIL', type=_cli_string, help='e-mail address')
    p.set_defaults(func=remove_email_from_blacklist)

    # Subparser 'email_blackist list'.
    p = s.add_parser(
        'list', help='list e-mail addresses in the blacklist',
        description=list_emails_in_blacklist.__doc__.split('\n\n')[0])
    p.set_defaults(func=list_emails_in_blacklist)

    args = parser.parse_args()

    try:
        args.func(**{k: v for k, v in vars(args).items()
                     if k not in ('func', 'subcommand')})
    except UserError as e:
        parser.error(unicode(e))


if __name__ == '__main__':
    main()
