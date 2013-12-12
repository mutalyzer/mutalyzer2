"""
Update the database with mapping information for the given gene or genomic
reference.
"""


import argparse

from .. import output
from .. import mapping


def import_gene(database, gene):
    """
    Update the database with information from the UCSC.

    .. todo: Also report how much was added/updated.
    """
    o = output.Output(__file__)
    o.addMessage(__file__, -1, 'INFO',
                 'Starting UCSC mapping data update (gene: %s)' % gene)

    updater = mapping.UCSCUpdater(database)
    updater.load(gene)
    updater.merge()

    o.addMessage(__file__, -1, 'INFO',
                 'UCSC mapping data update end (gene: %s)' % gene)


def import_reference(database, reference):
    """
    Update the database with information from the given reference.

    .. todo: Also report how much was added/updated.

    .. note: Currently no exon locations are supported, this has only been
       tested on mtDNA.
    """
    o = output.Output(__file__)
    o.addMessage(__file__, -1, 'INFO',
                 'Starting reference mapping data update (reference: %s)' % reference)

    updater = mapping.ReferenceUpdater(database)
    updater.load(reference, o)
    updater.merge()

    o.addMessage(__file__, -1, 'INFO',
                 'Reference mapping data update end (reference: %s)' % reference)


def main():
    """
    Command-line interface to the mapping importer.
    """
    database_parser = argparse.ArgumentParser(add_help=False)
    database_parser.add_argument(
        '-d', '--database', metavar='DATABASE', dest='database',
        default='hg19', help='database to import to (default: hg19)')

    parser = argparse.ArgumentParser(
        description='Mutalyzer mapping importer.',
        epilog='This program is intended to be run manually whenever '
        'transcript mappings for specific genes are required that are not '
        'yet in our database (i.e., they are not yet available from the '
        'NCBI, or they are mtDNA genes). It will not overwrite '
        'transcript/version entries that are already in our database.',
        parents=[database_parser])

    subparsers = parser.add_subparsers(
        title='subcommands', dest='subcommand', help='subcommand help')

    p = subparsers.add_parser(
        'gene', help='import gene', parents=[database_parser],
        description='Import gene mapping from the UCSC.')
    p.add_argument(
        'gene', metavar='GENE_SYMBOL',
        help='gene to import all transcript mappings for from the UCSC '
        'database (example: TTN)')

    p = subparsers.add_parser(
        'reference', help='import reference', parents=[database_parser],
        description='Import genomic reference file')
    p.add_argument(
        'reference', metavar='FILE',
        help='genomic reference to import all genes from (example: '
        'NC_012920.1)')

    args = parser.parse_args()
    if args.subcommand == 'gene':
        import_gene(args.database, args.gene)
    if args.subcommand == 'reference':
        import_reference(args.database, args.reference)


if __name__ == '__main__':
    main()
