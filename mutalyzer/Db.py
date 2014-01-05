#!/usr/bin/python

"""
Module for database access.
The Db class is a superclass of the rest of the classes and should not be
used as such. The superclass mainly consists of a wrapper for SQL
statements.

@requires: MySQLdb
@requires: types
@requires: time
@requires: os
"""

#Public classes:
#    - Db      ; Log in to a database and keep it open for queries.
#    - Mapping ; Mapping of transcripts and genes.
#    - Cache   ; Cache administration.
#    - Batch   ; Batch checker.


import types
import warnings

import MySQLdb

from mutalyzer import util
from mutalyzer.config import settings


#
# Note that compound queries are split into single queries because of a bug
# in MySQLdb. The functions load_Update(),  merge_cdsUpdates() and
# merge_Update (search for MYSQL_BUG in this file) are affected and may be
# rewritten when this bug is fixed.
#

class Db():
    """
    Query the database server (and lazily keep a connection open to it).

    This class is subclassed below to create specific interfaces to the
    database.
    """
    def __init__(self, database, user, host):
        """
        Create an interface to the database.

        @arg database: Name of the database to use.
        @type database: str
        @arg user: User name for the database.
        @type user: str
        @arg host: Host name for the database.
        @type host: str
        """
        self._database = database
        self._user = user
        self._host = host

        # The connection to the database server is created lazily in the query
        # method.
        self._connection = None
    #__init__

    def _connect(self):
        """
        Connect to the database server.

        Note: We would like to automatically reconnect to the database server.
            This is especially useful for long-running processes such as the
            batch deamon, which would otherwise loose their connection on an
            event such as restarting the database server.
            The MySQL client libraries provide a reconnect option, but this
            is unfortunately not implemented in (most versions of) the Python
            MySQLdb module.
            Therefore we manually implement automatic reconnects in the query
            method.
            Also see Trac ticket #91.
        """
        self._connection = MySQLdb.connect(
                user=self._user, db=self._database, host=self._host)
    #_connect

    def query(self, statement):
        """
        Query the database.

        @arg statement: The statement that is to be queried
        @type statement: tuple (string, (args))

        @return: The result of the query
        @rtype: list
        """

        # Convert the arguments to a tuple.
        if type(statement[1]) != types.TupleType :
            args = statement[1],
        else :
            args = statement[1]

        # Escape the input to prevent SQL injections.
        escaped_args = []
        if args != (None,) : # Don't escape the empty string.
            for i in args :
                if i :
                    if type(i) in [types.StringType, types.UnicodeType]:
                        escaped_args.append(MySQLdb.escape_string(str(i)))
                    else :
                        escaped_args.append(i)
                else :
                    escaped_args.append(None)
        #if

        # Do the query, but first connect to the database server if needed.
        # This makes sure lost connections are re-created automatically (e.g.
        # in case the server was restarted for maintenance).
        try:
            cursor = self._connection.cursor()
            cursor.execute(statement[0], tuple(escaped_args))
        except (AttributeError, MySQLdb.OperationalError):
            self._connect()
            cursor = self._connection.cursor()
            cursor.execute(statement[0], tuple(escaped_args))

        result = cursor.fetchall()
        cursor.close()

        return result
    #query
#Db


class Counter(Db):
    """
    Database functions for the service counters.

    Special methods:
        - __init__() ; Initialise the class.

    Inherited methods from Db:
        - query(statement) ; General query function.

    SQL tables from internalDb:
        - Counter   ; Service count information.
    """
    def __init__(self):
        """
        Initialise the Db parent class. Use the internalDb.
        """
        Db.__init__(self, settings.MYSQL_DATABASE,
                    settings.MYSQL_USER,
                    settings.MYSQL_HOST)

    def increment(self, service, interface):
        """
        Increment the counter for given service and interface.

        SQL tables from internalDb:
            - Counter ; Service count information.
        """
        statement = """
            UPDATE `Counter` SET
              `count` = `count` + 1
            WHERE `service` = %s
            AND `interface` = %s;
        """, (service, interface)

        self.query(statement)
#Counter
