"""
Some utility functions for our simple migrations.

@todo: Perhaps this should be moved to the mutalyzer package.
"""


import sys
import MySQLdb


COLOR_INFO = '\033[32m'
COLOR_WARNING = '\033[33m'
COLOR_ERROR = '\033[31m'
COLOR_END = '\033[0m'


def print_color(message, color=None):
    if color is None:
        print message
    else:
        print color + message + COLOR_END


def info(message):
    print_color(message, COLOR_INFO)


def warning(message):
    print_color(message, COLOR_WARNING)


def error(message):
    print_color(message, COLOR_ERROR)


def fatal(message):
    error(message)
    sys.exit(1)


def db_connect(database):
    try:
        connection = MySQLdb.connect(host='localhost',
                                     user='mutalyzer',
                                     passwd='',
                                     db=database)
    except MySQLdb.Error as e:
        fatal('Error %d: %s' % (e.args[0], e.args[1]))
    return connection


def main(check, migrate):
    needed = check()
    if needed:
        warning('This migration is needed.')
        if len(sys.argv) > 1 and sys.argv[1] == 'migrate':
            print 'Performing migration.'
            migrate()
            print 'Performed migration.'
    else:
        info('This migration is not needed.')
