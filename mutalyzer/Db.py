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


import MySQLdb # connect(), escape_string()
import types   # TupleType
from mutalyzer import util

#
# Note that compound queries are split into single queries because of a bug
# in MySQLdb. The functions load_Update(),  merge_cdsUpdates() and
# merge_Update (search for MYSQL_BUG in this file) are affected and may be
# rewritten when this bug is fixed.
#

class Db() :
    """
    Log in to a database and keep it open for queries.

    Private variables:
        - __db ; Interface to the database.

    Special methods:
        - __init__(dbName, mySqlUser, mySqlHost) ; Do the login.

    Public methods:
        - query(statement) ; General query function.
    """

    def __init__(self, dbName, mySqlUser, mySqlHost) :
        """
        Log in to the database.

        Private variables (altered):
            - __db ; The interface to the database.

        @arg dbName: The name of the database to use
        @type dbName: string
        @arg mySqlUser: User name for the database
        @type mySqlUser: string
        @arg mySqlHost: Host name for the database
        @type mySqlHost: string
        """

        self.__db = MySQLdb.connect(user = mySqlUser, db = dbName,
                                    host = mySqlHost)
    #__init__

    def query(self, statement) :
        """
        Query the database.

        Private variables:
            - __db ; Interface to the database.

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

        # And do the query.
        cursor = self.__db.cursor()
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
        - __init__(build, config) ; Initialise the class.

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

    def __init__(self, build, config) :
        """
        Initialise the Db parent class. Use the local database for a certain
        build.

        @arg build: The version of the mapping database
        @type build: string
        @arg config: Configuration variables
        @type config: class instance
        """

        Db.__init__(self, build, config.LocalMySQLuser, config.LocalMySQLhost)
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

    def getAllFields(self, mrnaAcc, version):
        """
        Get all Fields of an accession number and version number.
        If the version number is None, use the "newest" version number.

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
                        start, stop,
                        cds_start, cds_stop,
                        exon_starts, exon_stops,
                        gene, chromosome,
                        orientation, protein
                from Mapping
        """
        if version is None:
            q += """
                where transcript = %s
                order by version desc, chromosome asc;
                """
            statement = (q, mrnaAcc)
        else:
            q += """
                where transcript = %s and
                      version = %s
                order by chromosome asc;
                """
            statement = q, (mrnaAcc, version)

        return self.query(statement)[0]

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

    def get_TranscriptsByGeneName(self, geneName) :
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
            SELECT transcript, version
                FROM Mapping
                WHERE gene = %s;
        """, geneName

        ret = self.query(statement)
        if ret :
            l = []
            for i in ret :
                l.append(i[0] + '.' + str(i[1]))
            return l
        #if
        return []
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
            SELECT AccNo
              FROM ChrName
              WHERE name = %s;
        """, name

        ret = self.query(statement)
        if ret :
            return ret[0][0]
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

        @todo: Return number of entries added/updated.
        """
        statement = """
            CREATE TABLE IF NOT EXISTS MappingTemp LIKE Mapping;
        """, None
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

    def ncbi_import_exon(self, transcript, start, stop, cds_start, cds_stop,
                         protein):
        """
        Import an exon mapping in a temporary table.

        SQL tables from dbNames:
            - Exons ; Exon mappings from NCBI.
        """
        statement = """
            INSERT IGNORE INTO Exons
                (transcript, start, stop, cds_start, cds_stop, protein)
            VALUES
                (%s, %s, %s, %s, %s, %s);
        """, (transcript, start, stop, cds_start, cds_stop, protein)

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
        """
        statement = """
            SET group_concat_max_len = 32768;
        """, None
        self.query(statement)

        statement = """
            DROP TABLE IF EXISTS MappingTemp;
        """, None
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
                CONCAT('chr', T.chromosome) as chromosome,
                T.orientation as orientation,
                MIN(T.start) as start,
                MAX(T.stop) as stop,
                MAX(E.cds_start) as cds_start,
                MAX(E.cds_stop) as cds_stop,
                GROUP_CONCAT(DISTINCT E.start ORDER BY E.start ASC) as exon_starts,
                GROUP_CONCAT(DISTINCT E.stop ORDER BY E.stop ASC) as exon_stops,
                MAX(E.protein) as protein,
                'NCBI' as source
            FROM Transcripts as T, Genes as G, Exons as E
            WHERE T.gene_id = G.id AND T.name = E.transcript
            GROUP BY T.name;
        """, None
        self.query(statement)
    #ncbi_aggregate_mapping
#Mapping


class Cache(Db) :
    """
    Database functions for cache administration.

    Special methods:
        - __init__(config) ; Initialise the class.

    Public methods:
        - insertGB(accNo, GI, fileHash, ChrAccVer, ChrStart, ChrStop,
          orientation, url) ; Insert info about a GenBank record.
        - updateHash(accNo, fileHash) ; Update the hash of an accession number.
        - getGBFromLoc(ChrAccVer, ChrStart, ChrStop, orientation) ; Get the
                                    accession number from slicing information.
        - getGBFromHash(fileHash) ; Get the accession number from its hash.
        - getGBFromGI(GI)         ; Get the accession number from its GI
                                    number.
        - getLoc(accNo)           ; Get the slicing information of an
                                    accession number.
        - getHash(accNo)          ; Get the hash of a GenBank record.
        - getUrl(accNo)           ; Get the URL of an accession number.

    Inherited methods from Db:
        - query(statement) ; General query function.

    SQL tables from internalDb:
        - GBInfo ; Information about cached and uploaded GenBank files.
    """

    def __init__(self, config) :
        """
        Initialise the Db parent class. Use the internalDb.

        @arg config: Configuration variables
        @type config: class instance
        """

        Db.__init__(self, config.internalDb, config.LocalMySQLuser,
                    config.LocalMySQLhost)
    #__init__

    def insertGB(self, accNo, GI, fileHash, ChrAccVer, ChrStart,
                 ChrStop, orientation, url) :
        """
        Insert information about a GenBank record in the internal database.

        The accNo and fileHash arguments are mandatory.
            - If the record is a normal RefSeq, then the GI number should be
            provided.
            - If the record is a chromosome slice, then the ChrAccVer, ChrStart,
            ChrStop and orientation variables should be specified.
            - If the record is downloaded from the internet, the url should be
            provided.
            - If all fields except the mandatory ones are empty, the record is
            assumed to be uploaded.

        SQL tables from internalDb (altered):
            - GBInfo ; Information about cached and uploaded GenBank files.

        @arg accNo: The name associated with this record
        @type accNo: string
        @arg GI: The GI number (if available)
        @type GI: string
        @arg fileHash: The hash of the content of the record
        @type fileHash: string
        @arg ChrAccVer: The accession number of the chromosome (if available)
        @type ChrAccVer: string
        @arg ChrStart: The start of the record in chromosomal
                          coordinates (if available)
        @type ChrStart: integer
        @arg ChrStop: The end of the record in chromosomal coordinates
                      (if available)
        @type ChrStop: integer
        @arg orientation: The orientation of the record relative to the
                          chromosome (if available) (1 = forward,
                          2 = reverse complement)
        @type orientation: integer
        @arg url: The originating URL (if available)
        @type url: string
        """

        statement = """
            INSERT INTO GBInfo
             (AccNo, GI, hash, ChrAccVer, ChrStart, ChrStop, orientation, url)
              VALUES (%s, %s, %s, %s, %s, %s, %s, %s);
        """, (accNo, GI, fileHash, ChrAccVer, ChrStart, ChrStop, orientation,
              url)

        self.query(statement)
    #insertGB

    def insertLRG(self, accNo, fileHash, url):
        """
        Insert information about a LRG record in the internal database.

        See insertGB() for more information.

        @arg accNo: The name associated with this record
        @type accNo: string
        @arg fileHash: The hash of the content of the record
        @type fileHash: string
        @arg url:  The originating URL (if available)
        @type url: string
        """

        statement = """
            INSERT INTO GBInfo
             (AccNo, GI, hash, ChrAccVer, ChrStart, ChrStop, orientation, url)
              VALUES (%s, %s, %s, %s, %s, %s, %s, %s);
        """, (accNo, None, fileHash, None, None, None, None, url)

        self.query(statement)
    #insertLRG


    def updateHash(self, accNo, fileHash) :
        """
        Update the hash of an accession number.

        SQL tables from internalDb (altered):
            - GBInfo ; Information about cached and uploaded GenBank files.

        @arg accNo: The accession number of a GenBank record
        @type accNo: string
        @arg fileHash: The hash of a GenBank record
        @type fileHash: string
        """

        statement = """
            UPDATE GBInfo
              SET hash = %s
              WHERE AccNo = %s;
        """, (fileHash, accNo)

        self.query(statement)
    #updateHash

    def getGBFromLoc(self, ChrAccVer, ChrStart, ChrStop, orientation) :
        """
        Get the accession number from a chromosomic location, used
        to make a slice, typically this only affects UD-numbers.

        SQL tables from internalDb:
            - GBInfo ; Information about cached and uploaded GenBank files.

        @arg ChrAccVer: The accession number of the chromosome
        @type ChrAccVer: string
        @arg ChrStart: Start position of the slice
        @type ChrStart: integer
        @arg ChrStop: End position of the slice
        @type ChrStop: integer
        @arg orientation: Orientation of the slice:
                          1. Forward
                          2. Reverse complement
        @type orientation: integer

        @return: The accession number
        @rtype: string
        """

        statement = """
            SELECT AccNo
              FROM GBInfo
              WHERE ChrAccVer = %s
              AND ChrStart = %s
              AND ChrStop = %s
              AND orientation = %s;
        """, (ChrAccVer, ChrStart, ChrStop, orientation)

        ret = self.query(statement)
        if ret :
            return ret[0][0]
        return None
    #getGBFromLoc

    def getGBFromHash(self, fileHash) :
        """
        Get the accession number from its hash.

        SQL tables from internalDb:
            - GBInfo ; Information about cached and uploaded GenBank files.

        @arg fileHash: The hash of a GenBank record
        @type fileHash: string

        @return: The accession number
        @rtype: string
        """

        statement = """
            SELECT AccNo
              FROM GBInfo
              WHERE hash = %s;
        """, fileHash

        ret = self.query(statement)
        if ret :
            return ret[0][0]
        return None
    #getGBFromHash

    def getGBFromGI(self, GI) :
        """
        Get the accession number from its GI number, typically this only
        affects RefSeq sequences.

        SQL tables from internalDb:
            - GBInfo ; Information about cached and uploaded GenBank files.

        @arg GI: The GI number of a GenBank record
        @type GI: string

        @return: The accession number
        @rtype: string
        """

        statement = """
            SELECT AccNo
              FROM GBInfo
              WHERE GI = %s;
        """, GI

        ret = self.query(statement)
        if ret :
            return ret[0][0]
        return None
    #getGBFromGI

    def getGBSince(self, created_since):
        """
        Get all accession number entries with creation date {created_since}
        or later.

        SQL tables from internalDb:
            - GBInfo ; Information about cached and uploaded GenBank files.

        @arg created_since: Only entries with later creation dates are returned.
        @type created_since: datatime.datetime

        @return: The accession number
        @rtype: string
        """
        statement = """
            SELECT AccNo, GI, hash, ChrAccVer, ChrStart,
                   ChrStop, orientation, url, created
            FROM GBInfo
            WHERE created >= %s;
        """, created_since

        return self.query(statement)
    #getGBSince

    def getLoc(self, accNo) :
        """
        Get the slicing information of an accession number, typically this
        only affects UD numbers.

        SQL tables from internalDb:
            - GBInfo ; Information about cached and uploaded GenBank files.

        @arg accNo: The accession number of a genbank record
        @type accNo: string

        @return: The slicing information:
                   - ChrAccVer   ; Accession number of the chromosome
                   - ChrStart    ; Start position of the slice
                   - ChrStop     ; End position of the slice
                   - orientation ; Orientation of the slice (1 = forward,
                                   2 = reverse complement)
        @rtype: list
        """

        statement = """
            SELECT ChrAccVer, ChrStart, ChrStop, orientation
              FROM GBInfo
              WHERE AccNo = %s;
        """, accNo

        ret = self.query(statement)
        if ret :
            return list(ret[0])
        return None
    #getLoc

    def getHash(self, accNo) :
        """
        Get the hash of a GenBank record identified by an accession number.

        SQL tables from internalDb:
            - GBInfo ; Information about cached and uploaded GenBank files.

        @arg accNo: The accession number of a genbank record
        @type accNo: string

        @return: The hash of the GenBank record
        @rtype: string
        """

        statement = """
            SELECT hash
              FROM GBInfo
              WHERE AccNo = %s;
        """, accNo

        ret = self.query(statement)
        if ret :
            return ret[0][0]
        return None
    #getHash

    def getUrl(self, accNo) :
        """
        Get the URL of an accession number, typically this only affects
        uploaded UD numbers.

        SQL tables from internalDb:
            - GBInfo ; Information about cached and uploaded GenBank files.

        @arg accNo: The accession number of a genbank record
        @type accNo: string

        @return: The URL of the GenBank record
        @rtype: string
        """

        statement = """
            SELECT url
              FROM GBInfo
              WHERE AccNo = %s;
        """, accNo

        ret = self.query(statement)
        if ret :
            return ret[0][0]
        return None
    #getHash

    def getGI(self, accNo) :
        """
        Get the GI number that is connected to the accession number.

        SQL tables from internalDb:
            - GBInfo ; Information about cached and uploaded GenBank files.

        @arg accNo: The accession number
        @type accNo: string

        @return: GI number
        @rtype: string
        """

        statement = """
            SELECT GI
              FROM GBInfo
              WHERE AccNo = %s;
        """, accNo

        ret = self.query(statement)
        if ret :
            return ret[0][0]
        return None
    #getGI

    def getProtAcc(self, mrnaAcc) :
        """
        Gets the protein accession number for the given mRNA accession
        number.

        SQL tables from internalDb:
            - Link ; mRNA and associated protein IDs.

        @arg mrnaAcc: The ID of an mRNA
        @type mrnaAcc: string

        @return: The protein accession number
        @rtype: string
        """

        statement = """
            SELECT protAcc
              FROM Link
              WHERE mrnaAcc = %s;
        """, mrnaAcc

        ret = self.query(statement)
        if ret :
            return ret[0][0]
        return None
    #getProtAcc

    def getmrnaAcc(self, protAcc) :
        """
        Gets the mRNA accession number for a given protein accession number.

        SQL tables from internalDb:
            - Link ; mRNA and associated protein IDs.

        @arg protAcc: The protein ID
        @type protAcc: string

        @return: The mRNA accession number
        @rtype: string
        """

        statement = """
            SELECT mrnaAcc
              FROM Link
              WHERE protAcc = %s;
        """, protAcc

        ret = self.query(statement)
        if ret :
            return ret[0][0]
        return None

    #getProtAcc

    def insertLink(self, mrnaAcc, protAcc) :
        """
        Inserts the given mRNA and protein accession numbers into the Link
        table.

        SQL tables from internalDb:
            - Link ; mRNA and associated protein IDs.

        @arg protAcc: The protein ID
        @type protAcc: string
        @arg mrnaAcc: The ID of an mRNA
        @type mrnaAcc: string
        """

        statement = """
            INSERT INTO Link
              VALUES (%s, %s);
        """, (mrnaAcc, protAcc)

        self.query(statement)
    #insertLink
#Cache

class Batch(Db) :
    """
    Database functions for the batch checker.

    Special methods:
        - __init__(config) ; Initialise the class.

    Public methods:
        - isJobListEmpty()     ; See if there are active jobs.
        - addJob(outputFilter, email, fromHost); Add a job and give it a unique
                                                 ID.
        - getJobs()            ; Get a list of active jobs.
        - removeJob(jobID)     ; Remove a job and return information about
                                 the job submitter.
        - addToQueue(jobID, accNo, gene, variant) ; Add a request belonging to a
                                                    certain job to the queue.
        - getFromQueue(jobID)  ; Get a request belonging to a certain job
                                 from the queue.

    Inherited methods from Db:
        - query(statement) ; General query function.

    SQL tables from internalDb:
        - BatchJob   ; Job information.
        - BatchQueue ; Requests.
    """

    def __init__(self, config) :
        """
        Initialise the Db parent class. Use the internalDb.

        @arg config: Configuration variables
        @type config: class instance
        """

        Db.__init__(self, config.internalDb, config.LocalMySQLuser,
                    config.LocalMySQLhost)
    #__init__

    def isJobListEmpty(self) :
        """
        See if there are active jobs.

        SQL tables from internalDb:
            - BatchJob ; Job information.

        @return: False if there are active jobs, True otherwise
        @rtype: boolean
        """

        statement = """
            SELECT COUNT(*)
              FROM BatchJob;
        """, None

        if int(self.query(statement)[0][0]) :
            return False
        return True
    #isJobListEmpty

    def entriesLeftForJob(self, jobID):
        """
        Count the number of entries of a job that are still to be processed

        SQL tables from internalDB:
            - BatchQueue  ; Queue information

        @arg jobID: The JobID of interest
        @type jobID: string

        @return: The number of entries
        @rtype: integer
        """
        statement = """
            SELECT COUNT(*)
              FROM BatchQueue
              WHERE JobID = %s;
        """, jobID

        return int(self.query(statement)[0][0])
    #entriesLeftForJob


    def addJob(self, outputFilter, email, fromHost, jobType, Arg1):
        """
        Add a job and give it a unique ID.

        SQL tables from internalDb (altered):
            - BatchJob ; Job information.

        @arg outputFilter: Output settings for all requests in this job
        @type outputFilter: string
        @arg email: Contact information of the submitter
        @type email: string
        @arg jobType: The type of batch job
        @type jobType: string
        @arg Arg1: Possible argument.
        @type Arg1: string

        @return: A job ID
        @rtype: integer
        """
        jobID = util.generate_id()
        statement = """
            INSERT INTO BatchJob
              VALUES (%s, %s, %s, %s, %s, %s)
        """, (jobID, outputFilter, email, fromHost, jobType, Arg1)

        self.query(statement)
        return jobID
    #addJob

    def getJobs(self) :
        """
        Get a list of active jobs.

        SQL tables from internalDb:
            - BatchJob ; Job information.

        @return: List of tuples (job ID, job Type)
        @rtype: list
        """

        statement = """
            SELECT JobID, JobType, Arg1
              FROM BatchJob;
        """, None

        return self.query(statement)
    #getJobs

    def removeJob(self, jobID) :
        """
            Remove a job (because the queue for this job is empty) and return
            information needed to alert the job submitter.

            SQL tables from internalDb (altered):
                - BatchJob ; Job information.

            @arg jobID: Identifier of a job
            @type jobID: string

            @return: Data for the job submitter
            @rtype: triple
        """

        # First retrieve all information about this job.
        statement = """
            SELECT EMail, Filter, FromHost
              FROM BatchJob
              WHERE JobID = %s;
        """, jobID
        data = self.query(statement)[0]

        # Remove the job.
        statement = """
            DELETE
              FROM BatchJob
              WHERE JobID = %s;
        """, jobID

        self.query(statement)
        return data
    #removeJob

    def addToQueue(self, jobID, inputl, flag):
        """
        Add a request belonging to a certain job to the queue.

        SQL tables from internalDb (altered):
            - BatchQueue ; Requests.

        @arg jobID: Identifier of a job
        @type jobID: string
        @arg inputl: The input line of an entry
        @type inputl: string
        """

        # The first value (QueueID) will be auto increased by MySQL.
        statement = """
            INSERT INTO BatchQueue
              VALUES (%s, %s, %s, %s);
        """, (None, jobID, inputl, flag)

        self.query(statement)
    #addToQueue

    def updateBatchDb(self, jobID, old, new, flag, whereNot):
        """
        Update the Entries of a BatchJob. This is used to alter
        batch entries that would otherwise take a long time to process.
        e.g. a batch job with a lot of the same accession numbers without
        version numbers would take a long time because mutalyzer would
        fetch the file from the NCBI for each entry. A database update
        over all entries with the same accession number speeds up the
        job considerably.

        SQL tables from internalDb (altered):
            - BatchQueue ; Requests.

        @arg jobID: Identifier of a job
        @type jobID: string
        @arg old: String to be replaced
        @type old: string
        @arg new: String to replace old with
        @type new: string
        @arg flag: The reason of subsitution
        @type flag: string
        @arg whereNot: A negative selector to prevent false positives
        @type whereNot: string
        """
        #update whereNot to escape parenthesis
        whereNot = whereNot.replace("(","[(]").replace(")","[)]")
        statement = """
            update `BatchQueue` set
               `Input` = REPLACE(`Input`, %s, %s),
               `Flags` = CONCAT(IFNULL(`Flags`,""), %s)
               WHERE `JobID` = %s and NOT
               `Input` RLIKE %s;
        """, (old, new, flag, jobID, whereNot)

        self.query(statement)

    def skipBatchDb(self, jobID, where, flag):
        """
        Flag batch entries to be skipped. This is used if it is certain
        that an entry will cause an error, or that its output is ambiguous.

        SQL tables from internalDB (alterd):
            - BatchQueue  ; Requests

        @arg jobID: Identifier of a job
        @type jobID: string
        @arg where: Look for occurencus of this string
        @type where: string
        @arg flag: The reason of skipping
        @type flag: string
        """
        #update where to escape parenthesis
        where = where.replace("(","[(]").replace(")","[)]")

        statement = """
            update `BatchQueue` set
                `Flags` = CONCAT(IFNULL(`Flags`, ""), %s)
                where `JobID` = %s AND
                `Input` RLIKE %s;
        """, (flag, jobID, where)

        self.query(statement)


    def getFromQueue(self, jobID) :
        """
        Get a request belonging to a certain job from the queue. If a
        request is found, remove it from the queue and return it. Otherwise
        return nothing.

        SQL tables from internalDb (altered):
            - BatchQueue ; Requests.

        @arg jobID: Identifier of a job
        @type jobID: string

        @return:
            - accNo   ; The accession number of a request
            - gene    ; The gene and transcript variant information
            - variant ; The variant
        @rtype: triple
        """

        # To optimize this query, make sure to have two indices on the
        # table:
        # - UNIQUE INDEX (QueueID)
        # - INDEX (JobID, QueueID)
        statement = """
            SELECT QueueID, Input, Flags
            FROM BatchQueue
            WHERE QueueID = (
                SELECT MIN(QueueID)
                FROM BatchQueue
                GROUP BY JobID
                HAVING JobID = %s
            );
        """, jobID

        results = self.query(statement)
        if results :
            queueID, inputl, flags = results[0]
        else :
            return None, None

        # We have found a request, so remove it from the queue.
        statement = """
            DELETE
              FROM BatchQueue
              WHERE QueueID = %s;
        """, queueID

        self.query(statement)
        return inputl, flags
    #getFromQueue
#Batch
