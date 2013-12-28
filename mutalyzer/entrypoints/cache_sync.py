"""
Synchronize the database cache with other Mutalyzer instances.

This program is intended to be run daily from cron. Example:

    25 5 * * *  mutalyzer-cache-sync 'http://dom1/?wsdl' 'http://dom1/{file}' -H 7
    55 5 * * *  mutalyzer-cache-sync 'http://dom2/?wsdl' 'http://dom2/{file}' -H 7
"""


import argparse

from .. import output
from .. import sync


def sync_cache(remote_wsdl, url_template, history=7):
    """
    Synchronize the database cache with other Mutalyzer instances.
    """
    output = output.Output(__file__)

    cache_sync = sync.CacheSync(output)
    cache_sync.sync_with_remote(remote_wsdl, url_template, history)


def main():
    """
    Command-line interface to the cache synchronizer.
    """
    parser = argparse.ArgumentParser(
        description='Mutalyzer cache synchronizer.')
    parser.add_argument(
        'wsdl', metavar='WSDL',
        help='location of the remote WSDL description')
    parser.add_argument(
        'url_template', metavar='URL_TEMPLATE',
        help='URL for remote downloads, in which the filename is to be '
        'substituted for {{file}}')
    parser.add_argument(
        '-H', '--history', metavar='DAYS', dest='history', type=int,
        default=7, help='number of days to go back in the remote cache '
        '(default: 7)')

    args = parser.parse_args()
    sync_cache(args.wsdl, args.url_template, history=args.history)


if __name__ == '__main__':
    main()
