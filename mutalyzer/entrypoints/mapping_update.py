"""
Update the database with mapping information from the NCBI.

This program is intended to be run daily from cron. Example:

    25 6 * * *  mutalyzer-mapping-update hg19 /tmp/seq_gene.md reference
"""


import argparse

from .. import output
from .. import mapping


def update_mapping(database, mapping_file, assembly):
    """
    Update the database with information from the NCBI.

    .. todo: Also report how much was added/updated.
    """
    o = output.Output(__file__)
    o.addMessage(__file__, -1, 'INFO',
                 'Starting NCBI mapping data update')

    updater = mapping.NCBIUpdater(database)
    updater.load(mapping_file, assembly)
    updater.merge()

    o.addMessage(__file__, -1, 'INFO', 'NCBI mapping data update end')


def main():
    """
    Command-line interface to the mapping updater.
    """
    parser = argparse.ArgumentParser(
        description='Mutalyzer mapping updater.')

    parser.add_argument(
        'mapping', metavar='FILE',
        help='Path to the NCBI mapping information (example: seq_gene.md)')
    parser.add_argument(
        'assembly', metavar='ASSEMBLY',
        help='use only entries from this assembly (this is the group_name '
        'column in the NCBI mapping file)')
    parser.add_argument(
        '-d', '--database', metavar='DATABASE', dest='database',
        default='hg19', help='database to update (default: hg19)')

    args = parser.parse_args()
    update_mapping(args.database, args.mapping, args.assembly)


if __name__ == '__main__':
    main()
