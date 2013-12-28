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


class Mapping(Db) :
    """
    Database functions for mapping of transcripts and genes.

    Special methods:
        - __init__(build) ; Initialise the class.

    Public methods:
        - get_protAcc(mrnaAcc)      ; Query the database for a protein ID.
        - get_NM_info(mrnaAcc)      ; Retrieve various data for an NM number.
        - get_NM_version(mrnaAcc)   ; Get the version number of an accession
                                      number.
        - get_Transcripts(chrom, p1, p2, overlap) ; Get a list of transcripts,
                                      given a chromosome and a range.
        - get_GeneName(mrnaAcc)     ; Get the gene name, given an NM number.
        - isChrom(name)             ; Check whether we know this name to be
                                      a chromosome name.

    Inherited methods from Db:
        - query(statement) ; General query function.

    SQL tables from dbNames:
        - Mapping; Accumulated mapping info.
    """

    def __init__(self, build) :
        """
        Initialise the Db parent class. Use the local database for a certain
        build.

        @arg build: The version of the mapping database
        @type build: string
        """
        Db.__init__(self, build, settings.MYSQL_USER, settings.MYSQL_HOST)
    #__init__

    def get_NM_version(self, mrnaAcc) :
        """
        Get the version number of an accession number.

        SQL tables from dbNames:
            - Mapping ; Accumulated mapping info.

        @arg mrnaAcc: The ID of an mRNA
        @type mrnaAcc: string

        @return: The version number
        @rtype: integer
        """

        statement = """
            SELECT version
              FROM Mapping
              WHERE transcript = %s;
        """, mrnaAcc

        return [int(i[0]) for i in self.query(statement)]
    #get_NM_version

    def getAllFields(self, mrnaAcc, version=None, selector=None, selector_version=None):
        """
        Get all Fields of an accession number and version number.
        If the version number is None, use the "newest" version number.

        Optionally also provide a gene and transcript version to select.

        SQL tables from dbNames:
            - Mapping ; Accumulated mapping info.

        @arg mrnaAcc: The ID of an mRNA
        @type mrnaAcc: string
        @arg version: The version number
        @type version: integer

        @return: The version number
        @rtype: integer

        @todo: The 'order by chrom asc' is a quick hack to make sure we first
            get a primary assembly mapping instead of some haplotype mapping
            for genes in the HLA cluster.
            A better fix is to return the entire list of mappings, and/or
            remove all secondary mappings for the HLA cluster.
            See also test_converter.test_hla_cluster and bug #58.
        """
        q = """
                select  transcript,
                        selector, selector_version,
                        start, stop,
                        cds_start, cds_stop,
                        exon_starts, exon_stops,
                        gene, chromosome,
                        orientation, protein
                from Mapping
        """

        where = []
        order = []
        args = []

        if version is None:
            where.append('transcript = %s')
            order.append('version desc')
            order.append('chromosome asc')
            args.append(mrnaAcc)
        else:
            where.append('transcript = %s')
            where.append('version = %s')
            order.append('chromosome asc')
            args.append(mrnaAcc)
            args.append(version)

        # The fallback to NULL selector (and selector_version) is necessary
        # to accept transcript selection in NM references. To be safe, we also
        # prefer entries with a match over entries with NULL.
        # Example: NM_017780.2(CHD7_v001):c.109A>T
        if selector is not None:
            where.append('(selector = %s or (selector IS NULL and selector_version IS NULL and gene = %s and %s = 1))')
            args.append(selector)
            args.append(selector)
            args.append(selector_version)
            order.append('(selector is null) asc')

        if selector_version is not None:
            where.append('(selector_version = %s or (selector_version IS NULL and selector_version IS NULL and gene = %s and %s = 1))')
            args.append(selector_version)
            args.append(selector)
            args.append(selector_version)
            order.append('(selector_version is null) asc')

        q += """
             where """ + ' AND '.join(where) + """
             order by """ + ', '.join(order) + ';'

        statement = q, tuple(args)
        try:
            return self.query(statement)[0]
        except IndexError:
            return None

    def get_Transcripts(self, chrom, p1, p2, overlap) :
        """
        Get a list of transcripts, given a chromosome and a range. If
        all transcripts that are hit should be returned, set overlap to 1,
        if only the transcripts that completely reside within a range
        should be returned, set overlap to 0.

        SQL tables from dbNames:
            - Mapping ; Accumulated mapping info.

        @arg chrom: The chromosome (coded as "chr1", ..., "chrY")
        @type chrom: string
        @arg p1: The position relative to the start of the chromosome
        @type p1: integer
        @arg p2: The position relative to the start of the chromosome
        @type p2: integer
        @arg overlap: Specify the behaviour of the selection:
          - 0 ; Return only the transcripts that completely fall in the
                range [p1, p2]
          - 1 ; Return all hit transcripts
        @type overlap: boolean

        @return: All accession numbers that are hit according to the overlap
        criterium
        @rtype: list
        """
        q = """
                select  transcript,
                        selector, selector_version,
                        start, stop,
                        cds_start, cds_stop,
                        exon_starts, exon_stops,
                        gene, chromosome,
                        orientation, protein,
                        version
                from Mapping
        """
        if overlap:
            q += """
              WHERE chromosome = %s AND
                    start <= "%s" AND
                    stop >= "%s";
            """
            statement = q, (chrom, p2, p1)

        else:
            q += """
              WHERE chromosome = %s AND
                    start >= "%s" AND
                    stop <= "%s";
            """
            statement = q, (chrom, p1, p2)

        return self.query(statement)
    #get_Transcripts

    def get_TranscriptsByGeneName(self, gene):
        """
            Get a list of transcripts, given a gene name.

            Arguments:
                geneName ; Name of a gene.

            SQL tables from dbNames:
                Mapping ; Accumulated mapping info.

            Returns:
                list ; A list of transcripts.
        """
        statement = """
                SELECT  transcript,
                        selector, selector_version,
                        start, stop,
                        cds_start, cds_stop,
                        exon_starts, exon_stops,
                        gene, chromosome,
                        orientation, protein,
                        version
                FROM Mapping
                WHERE gene = %s;
        """, gene

        return self.query(statement)
    #get_TranscriptsByGeneName

    def get_GeneName(self, mrnaAcc) :
        """
        Get the name of a gene, given a transcript identifier (NM number).

        SQL tables from dbNames:
            - Mapping ; Accumulated mapping info.

        @arg mrnaAcc: The ID of an mRNA
        @type mrnaAcc: string

        @return: The gene name
        @rtype: string
        """

        statement = """
            SELECT gene
              FROM Mapping
              WHERE transcript = %s;
        """, mrnaAcc

        ret = self.query(statement)
        if ret :
            return ret[0][0]
        return None
    #get_GeneName

    def isChrom(self, name) :
        """
        Check if the given name is a valid chromosome name.

        SQL tables from dbNames:
            - Mapping ; Accumulated mapping info.

        @arg name: The name to be tested
        @type name: string

        @return: True if the name is found to be a chromosome name, False
        otherwise
        @rtype: boolean
        """

        statement = """
            SELECT COUNT(*)
              FROM Mapping
              WHERE chromosome = %s;
        """, name

        if int(self.query(statement)[0][0]) > 0 :
            return True
        return False
    #isChrom

    def chromName(self, accNo) :
        """
        Get the name of a chromosome, given an accession number.

        SQL tables from dbNames:
            - ChrName ; Assembly release notes.

        @arg accNo: The accession number of a chromosome
        @type accNo: string

        @return: The name of a chromosome
        @rtype: string
        """

        statement = """
            SELECT name
              FROM ChrName
              WHERE AccNo = %s;
        """, accNo

        ret = self.query(statement)
        if ret :
            return ret[0][0]
        return None
    #chromName

    def chromAcc(self, name) :
        """
        Get the accession number of a chromosome, given a name.

        SQL tables from dbNames:
            - ChrName ; Assembly release notes.

        @arg name: The name of a chromosome
        @type name: string

        @return: The accession number of a chromosome
        @rtype: string
        """

        statement = """
            SELECT AccNo, organelle_type
              FROM ChrName
              WHERE name = %s;
        """, name

        ret = self.query(statement)
        if ret :
            return ret[0]
        return None
    #chromAcc

    def get_chromName(self, acc) :
        """
        Get the chromosome name, given a transcript identifier (NM number).

        SQL tables from dbNames:
            - Mapping ; Accumulated mapping info.

        @arg acc: The NM accession number (version NOT included)
        @type acc: string

        @return: The chromosome name (e.g. chr1)
        @rtype: string
        """

        statement = """
            SELECT chromosome
              FROM Mapping
              WHERE transcript = %s;
        """, acc

        ret = self.query(statement)
        if ret :
            return ret[0][0]
        return None
    #get_chromName

    def merge_update(self):
        """
        Merge existing mapping information with new mapping information, which
        should be in table 'MappingTemp'.

        The strategy is as follows. Existing mappings (accumulated by
        Mutalyzer in the past) that are not in the new mapping information are
        added to the new mapping information. The resulting set is used as the
        mapping information from now on.

        This way, we get the latest updates for existing mappings and keep old
        mappings not in the updated information.

        SQL tables from dbNames:
            - Mapping ; Accumulated mapping info.
            - MappingTemp ; New mapping info.
            - MappingBackup ; Backup of accumulated mapping info.

        @note: We temporarily suppress warnings during some queries, since
            they are expected and clutter the console output (e.g. warnings
            for existing tables).
        @todo: Return number of entries added/updated.
        """
        statement = """
            CREATE TABLE IF NOT EXISTS MappingTemp LIKE Mapping;
        """, None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.query(statement)

        statement = """
            INSERT INTO MappingTemp
            SELECT * FROM Mapping AS OldM
            WHERE NOT EXISTS (
                SELECT * FROM MappingTemp AS NewM
                WHERE OldM.transcript = NewM.transcript
                AND OldM.version = NewM.version
            );
        """, None
        self.query(statement)

        statement = """
            DROP TABLE IF EXISTS MappingBackup;
        """, None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.query(statement)

        statement = """
            RENAME TABLE Mapping TO MappingBackup, MappingTemp TO Mapping;
        """, None
        self.query(statement)
    #merge_update

    def ncbi_create_temporary_tables(self):
        """
        Create temporary tables to import NCBI mapping into.

        SQL tables from dbNames:
            - Genes ; Gene names from NCBI.
            - Transcripts ; Transcript mappings from NCBI.
            - Exons ; Exon mappings from NCBI.
        """
        self.ncbi_drop_temporary_tables()

        statement = """
            CREATE TABLE Genes (
                id varchar(20) NOT NULL DEFAULT '',
                name varchar(255) DEFAULT NULL,
                PRIMARY KEY (id)
            );
        """, None
        self.query(statement)

        statement = """
            CREATE TABLE Transcripts (
                name varchar(20) NOT NULL DEFAULT '',
                gene_id varchar(20) DEFAULT NULL,
                chromosome char(2) DEFAULT NULL,
                start int(11) DEFAULT NULL,
                stop int(11) DEFAULT NULL,
                orientation char(1) DEFAULT NULL,
                PRIMARY KEY (name,start)
            );
        """, None
        self.query(statement)

        statement = """
            CREATE TABLE Exons (
                transcript varchar(20) NOT NULL DEFAULT '',
                chromosome char(2) DEFAULT NULL,
                start int(11) DEFAULT NULL,
                stop int(11) DEFAULT NULL,
                cds_start int(11) DEFAULT NULL,
                cds_stop int(11) DEFAULT NULL,
                protein varchar(20) DEFAULT NULL,
                PRIMARY KEY (transcript,start)
            );
        """, None
        self.query(statement)
    #ncbi_create_temporary_table

    def ncbi_drop_temporary_tables(self):
        """
        Drop temporary tables used for importing NCBI mapping information.

        SQL tables from dbNames:
            - Genes ; Gene names from NCBI.
            - Transcripts ; Transcript mappings from NCBI.
            - Exons ; Exon mappings from NCBI.
        """
        statement = """
            DROP TABLE IF EXISTS Genes, Transcripts, Exons;
        """, None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.query(statement)
    #ncbi_drop_temporary_tables

    def ncbi_import_gene(self, id, name):
        """
        Import a (gene id, gene name) pair in a temporary table.

        SQL tables from dbNames:
            - Genes ; Gene names from NCBI.
        """
        statement = """
            INSERT IGNORE INTO Genes (id, name) VALUES (%s, %s);
        """, (id, name)

        self.query(statement)
    #ncbi_import_gene

    def ncbi_import_transcript(self, name, gene, chromosome, start, stop,
                               orientation):
        """
        Import a transcript mapping in a temporary table.

        SQL tables from dbNames:
            - Transcripts ; Transcript mappings from NCBI.
        """
        statement = """
            INSERT IGNORE INTO Transcripts
                (name, gene_id, chromosome, start, stop, orientation)
            VALUES
                (%s, %s, %s, %s, %s, %s);
        """, (name, gene, chromosome, start, stop, orientation)

        self.query(statement)
    #ncbi_import_transcript

    def ncbi_import_exon(self, transcript, chromosome, start, stop, cds_start,
                         cds_stop, protein):
        """
        Import an exon mapping in a temporary table.

        SQL tables from dbNames:
            - Exons ; Exon mappings from NCBI.
        """
        statement = """
            INSERT IGNORE INTO Exons
                (transcript, chromosome, start, stop, cds_start, cds_stop, protein)
            VALUES
                (%s, %s, %s, %s, %s, %s, %s);
        """, (transcript, chromosome, start, stop, cds_start, cds_stop, protein)

        self.query(statement)
    #ncbi_import_exon

    def ncbi_aggregate_mapping(self):
        """
        Aggregate gene, transcript and exon mapping information from the NCBI
        into one table.

        @note: Default MySQL value for group_concat_max_len is 1024, meaning
            that the GROUP_CONCAT aggregate function returns values of at most
            1024 bytes long. This is not enough (currently we need around 3000
            bytes), so we explicitely set this to a higher value.
        @note: We use MAX(E.protein) since MySQL does not have an ANY()
            aggregator.
        @note: Some genes (e.g. in the PAR) are mapped on both the X and Y
            chromosomes. Therefore, we group not only by transcript name, but
            also by chromosome, and add conditions on exon positions. The flaw
            here is that we miss genes that are mapped to two locations on one
            chromosome, but I don't think we have any of those.
        """
        statement = """
            SET group_concat_max_len = 32768;
        """, None
        self.query(statement)

        statement = """
            DROP TABLE IF EXISTS MappingTemp;
        """, None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.query(statement)

        statement = """
            CREATE TABLE MappingTemp LIKE Mapping;
        """, None
        self.query(statement)

        statement = """
            INSERT INTO MappingTemp
            SELECT
                G.name as gene,
                SUBSTRING(T.name FROM 1 FOR LOCATE('.', T.name) - 1) as transcript,
                SUBSTRING(T.name FROM LOCATE('.', T.name) + 1) as version,
                NULL as selector,
                NULL as selector_version,
                CONCAT('chr', T.chromosome) as chromosome,
                T.orientation as orientation,
                MIN(T.start) as start,
                MAX(T.stop) as stop,
                REPLACE(MIN(COALESCE(E.cds_start, 1000000000)), 1000000000, NULL) as cds_start,
                MAX(E.cds_stop) as cds_stop,
                GROUP_CONCAT(DISTINCT E.start ORDER BY E.start ASC) as exon_starts,
                GROUP_CONCAT(DISTINCT E.stop ORDER BY E.stop ASC) as exon_stops,
                MAX(E.protein) as protein,
                'NCBI' as source
            FROM Transcripts as T, Genes as G, Exons as E
            WHERE T.gene_id = G.id AND T.name = E.transcript
            AND E.chromosome = T.chromosome AND E.start >= T.start AND E.stop <= T.stop
            GROUP BY T.name, T.chromosome;
        """, None
        self.query(statement)
    #ncbi_aggregate_mapping

    def ucsc_load_mapping(self, transcripts, overwrite=False):
        """
        Load given transcripts into the 'MappingTemp' table.

        Todo: Don't overwrite NCBI/reference entries, but *do* overwrite
            existing UCSC entries.
        """
        statement = """
            DROP TABLE IF EXISTS MappingTemp;
        """, None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.query(statement)

        statement = """
            CREATE TABLE MappingTemp LIKE Mapping;
        """, None
        self.query(statement)

        for (gene, transcript, version, chromosome, orientation, start,
             stop, cds_start, cds_stop, exon_starts, exon_stops, protein) in transcripts:
            exon_starts = ','.join(str(i) for i in exon_starts)
            exon_stops = ','.join(str(i) for i in exon_stops)

            statement = """
                INSERT INTO `MappingTemp`
                  (`gene`, `transcript`, `version`, `chromosome`, `orientation`,
                   `start`, `stop`, `cds_start`, `cds_stop`, `exon_starts`, `exon_stops`,
                   `protein`, `source`)
                SELECT
                  %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s
                FROM DUAL WHERE NOT EXISTS
                  (SELECT * FROM `Mapping`
                   WHERE `transcript` = %s AND `version` = %s AND 1 = %s);
                """, (gene, transcript, version, chromosome, orientation, start,
                      stop, cds_start, cds_stop, exon_starts, exon_stops, protein,
                      'UCSC', transcript, version, int(overwrite) + 1)
            self.query(statement)
    #ucsc_load_mapping

    def reference_load_mapping(self, transcripts, overwrite=False):
        """
        Load given transcripts into the 'MappingTemp' table.

        Todo: Don't overwrite NCBI/UCSC entries, but *do* overwrite existing
            reference entries.
        """
        statement = """
            DROP TABLE IF EXISTS MappingTemp;
        """, None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.query(statement)

        statement = """
            CREATE TABLE MappingTemp LIKE Mapping;
        """, None
        self.query(statement)

        for t in transcripts:
            exon_starts = ','.join(str(i) for i in t['exon_starts'])
            exon_stops = ','.join(str(i) for i in t['exon_stops'])
            statement = """
                INSERT INTO `MappingTemp`
                  (`gene`, `transcript`, `version`, `selector`, `selector_version`,
                   `chromosome`, `orientation`, `start`, `stop`, `cds_start`, `cds_stop`,
                   `exon_starts`, `exon_stops`, `protein`, `source`)
                SELECT
                  %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 'reference'
                FROM DUAL WHERE NOT EXISTS
                  (SELECT * FROM `Mapping`
                   WHERE `transcript` = %s AND `version` = %s AND 1 = %s);
                """, (t['gene'], t['transcript'], t['version'], t['selector'],
                      t['selector_version'], t['chromosome'], t['orientation'],
                      t['start'], t['stop'], t['cds_start'], t['cds_stop'],
                      exon_starts, exon_stops,
                      t['protein'], t['transcript'], t['version'], int(overwrite) + 1)
            self.query(statement)
    #reference_load_mapping
#Mapping


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
