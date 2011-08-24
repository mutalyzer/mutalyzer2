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
#    - Remote  ; Retrieving updates for the mapping databases.
#    - Update  ; Updating the mapping databases.
#    - Cache   ; Cache administration.
#    - Batch   ; Batch checker.


import MySQLdb # connect(), escape_string()
import types   # TupleType
import time    # strftime()
import os      # os.remove()
import tempfile
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
        - map ; Accumulated mapping info.
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

    def get_protAcc(self, mrnaAcc) :
        """
        Query the database for a protein ID given an mRNA ID.

        SQL tables from dbNames:
            - map ; Accumulated mapping info.

        @arg mrnaAcc: The ID of an mRNA
        @type mrnaAcc: string

        @return: The protein ID
        @rtype: string
        """

        statement = """
            SELECT protAcc
              FROM map
              WHERE acc = %s;
        """, mrnaAcc

        return self.query(statement)[0][0]
    #get_protAcc

    def get_NM_info(self, mrnaAcc, version = None) :
        """
        Retrieve various data for an NM number.

        SQL tables from dbNames:
            - map ; Accumulated mapping info.

        @arg mrnaAcc: The ID of an mRNA
        @type mrnaAcc: string
        @arg version: version number of the accession number (not used)
        @type version: integer

        @return: 
                 - exonStarts ; List of exon start sites.
                 - exonEnds   ; List of exon end sites.
                 - cdsStart   ; Position of the start codon.
                 - cdsEnd     ; Position of the end codon.
                 - strand     ; Orientation of the gene (+ = forward,
                                                         - = reverse)
        @rtype: list
        """

        statement = """
            SELECT exonStarts, exonEnds, cdsStart, cdsEnd, strand
              FROM map
              WHERE acc = %s;
        """, mrnaAcc

        return self.query(statement)[0]
    #get_NM_info

    def get_NM_version(self, mrnaAcc) :
        """
        Get the version number of an accession number.

        SQL tables from dbNames:
            - map ; Accumulated mapping info.

        @arg mrnaAcc: The ID of an mRNA
        @type mrnaAcc: string

        @return: The version number
        @rtype: integer
        """

        statement = """
            SELECT version
              FROM map
              WHERE acc = %s;
        """, mrnaAcc

        return [int(i[0]) for i in self.query(statement)]
    #get_NM_version

    def getAllFields(self, mrnaAcc, version):
        """
        Get all Fields of an accession number and version number.
        If the version number is None, use the "newest" version number.

        SQL tables from dbNames:
            - map ; Accumulated mapping info.

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
                select  acc,
                        txStart, txEnd,
                        cdsStart, cdsEnd,
                        exonStarts, exonEnds,
                        geneName, chrom,
                        strand, protAcc
                from map
        """
        if version is None:
            q += """
                where acc = %s
                order by version desc, chrom asc;
                """
            statement = (q, mrnaAcc)
        else:
            q += """
                where acc = %s and
                      version = %s
                order by chrom asc;
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
            - map ; Accumulated mapping info.

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
                select  acc,
                        txStart, txEnd,
                        cdsStart, cdsEnd,
                        exonStarts, exonEnds,
                        geneName, chrom,
                        strand, protAcc,
                        version
                from map
        """
        if overlap:
            q += """
              WHERE chrom = %s AND
                    txStart <= "%s" AND
                    txEnd >= "%s";
            """
            statement = q, (chrom, p2, p1)

        else:
            q += """
              WHERE chrom = %s AND
                    txStart >= "%s" AND
                    txEnd <= "%s";
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
                map ; Accumulated mapping info.

            Returns:
                list ; A list of transcripts.
        """

        statement = """
            SELECT acc, version
                FROM map
                WHERE geneName = %s;
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
            - map ; Accumulated mapping info.

        @arg mrnaAcc: The ID of an mRNA
        @type mrnaAcc: string

        @return: The gene name
        @rtype: string
        """

        statement = """
            SELECT geneName
              FROM map
              WHERE acc = %s;
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
            - map ; Accumulated mapping info.

        @arg name: The name to be tested
        @type name: string

        @return: True if the name is found to be a chromosome name, False
        otherwise
        @rtype: boolean
        """

        statement = """
            SELECT COUNT(*)
              FROM map
              WHERE chrom = %s;
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
            - map ; Accumulated mapping info.

        @arg acc: The NM accession number (version NOT included)
        @type acc: string

        @return: The chromosome name (e.g. chr1)
        @rtype: string
        """

        statement = """
            SELECT chrom
              FROM map
              WHERE acc = %s;
        """, acc
        print acc
        ret = self.query(statement)
        if ret :
            return ret[0][0]
        return None
    #get_chromName
#Mapper

class Remote(Db) :
    """
    Database functions for retrieving updates for the mapping databases.

    Special methods:
        - __init__(config) ; Initialise the class.

    Public methods:
        - get_Update()        ; Retrieve new mapping info from the UCSC.

    Inherited methods from Db:
        - query(statement) ; General query function.

    SQL tables from dbNames:
        - gbStatus ; acc -> version mapping (NM to NM + version),
                     type, modDate
        - refGene  ; name -> geneName mapping (NM to gene name),
                     txStart, txEnd, cdsStart, cdsEnd, exonStarts,
                     exonEnds, chrom, strand.
        - refLink  ; mrnaAcc -> protAcc mapping (NM to NP).
    """

    def __init__(self, build, config) :
        """
        Initialise the Db parent class. Use the remote database for a
        certain build.

        Private variables (altered):
            - __config ; Configuration variables.

        @arg build: The version of the mapping database
        @type build: string
        @arg config: Configuration variables
        @type config: class instance
        """

        self.__config = config
        Db.__init__(self, build, config.RemoteMySQLuser, config.RemoteMySQLhost)
    #__init__

    def get_Update(self) :
        """
        Retrieve all mapping updates from the UCSC within a certain time
        window (defined in the configuration file) and gather the results
        into one mapping table.

        The results will be written to a temporary file to be imported in
        the local database with the load_Update() function.

        Return temporary filename used to store the results.

        @return: Filename used to store the results.
        @rtype: string

        SQL tables from dbNames:
            - gbStatus ; acc -> version mapping (NM to NM + version),
                         type, modDate
            - refGene  ; name -> geneName mapping (NM to gene name),
                         txStart, txEnd, cdsStart, cdsEnd, exonStarts,
                         exonEnds, chrom, strand.
            - refLink  ; mrnaAcc -> protAcc mapping (NM to NP).
        """

        statement = """
            SELECT DISTINCT acc, version, txStart, txEnd, cdsStart, cdsEnd,
                            exonStarts, exonEnds, name2 AS geneName, chrom,
                            strand, protAcc
              FROM gbStatus, refGene, refLink
              WHERE type = "mRNA"
              AND refGene.name = acc
              AND acc = mrnaAcc
              AND time >= DATE_SUB(CURDATE(), INTERVAL %s DAY);
        """, self.__config.UpdateInterval

        handle, filename = tempfile.mkstemp(text=True)

        # Convert the results to a tab delimited file.
        for i in self.query(statement) :
            for j in i :
                os.write(handle, str(j) + chr(0x09))  # 0x09 is a TAB.
            os.write(handle, '\n')
        #for

        os.close(handle)
        return filename
    #get_Update
#Remote

class Update(Db) :
    """
    Database functions for updating the mapping databases.

    Public methods:
        - load_Update()       ; Load new mapping info into the local database.
        - count_Updates()     ; Count the number of entries in the new
                                mapping info table.
        - backup_cdsUpdates() ; Make a backup of updates that overwrite the
                                old mapping info.
        - count_cdsUpdates()  ; Count the number of updates that overwrite
                                the old mapping info.
        - merge_cdsUpdates()  ; Merge the backup of old mapping info with the
                                other old info.
        - merge_Update()      ; Merge the new mapping info from the UCSC with
                                what we already have.

    Inherited methods from Db:
        - query(statement) ; General query function.

    SQL tables from dbNames:
        - map                ; Accumulated mapping info.
        - map_temp           ; Newly found data.
        - map_new            ; Merge of map_temp and map.
        - map_cdsBackup_temp ; Entries that were updated without an increment
                               of the version number.
        - map_cdsBackup      ; Merge of map_cdsBackup_temp and itself.
    """

    def __init__(self, build, config) :
        """
        Initialise the Db parent class. Use the remote database for a
        certain build.

        Private variables (altered):
            - __config ; Configuration variables.

        @arg build: The version of the mapping database
        @type build: string
        @arg config: Configuration variables
        @type config: class instance
        """

        self.__config = config
        Db.__init__(self, build, config.LocalMySQLuser, config.LocalMySQLhost)
    #__init__

    def load_Update(self, filename) :
        """
        Load the updates from the temporary file created by the get_Update()
        function and import it in the local database.

        @arg filename: Filename to read the update from.
        @type filename: string

        SQL tables from dbNames (altered):
            - map_temp ; Created and loaded with data from TempFile.

        SQL tables from dbNames:
            - map ; Accumulated mapping info.
        """

        # The statements in this function may be combined when MYSQL_BUG is
        # solved.

        statement = """
            CREATE TABLE map_temp LIKE map;
        """, None
        self.query(statement)
        statement = """
            LOAD DATA LOCAL INFILE %s
              INTO TABLE map_temp;
        """, filename

        self.query(statement)

        os.remove(filename)
    #load_Update

    def count_Updates(self) :
        """
        Count the number of updates. This function will only work if it
        is preceeded by the load_Update() function. Otherwise the map_temp
        table may not exist. This function can not be used after the
        merge_Update() function has been executed, since it drops the
        map_temp table.

        @return: The number of entries in the table of updated mapping info
        @rtype: integer
        """

        statement = """
            SELECT COUNT(*)
              FROM map_temp;
        """, None

        return int(self.query(statement)[0][0])
    #count_Updates

    def backup_cdsUpdates(self) :
        """
        Copy all mapping entries where there was an update, but no
        increment in the version number, to a backup table. Note that
        we use acc, version, txStart as the primary key because members
        of a gene family are mapped multiple times.

        SQL tables from dbNames (altered):
            - map_cdsBackup_temp ; Created and filled with entries that
                                   were updated without an increment of the
                                   version number.

        SQL tables from dbNames:
            - map      ; Accumulated mapping info.
            - map_temp ; Freshly downloaded mapping info.
        """

        statement = """
            CREATE TABLE map_cdsBackup_temp
              SELECT map.*
                FROM map, map_temp
                WHERE map.acc = map_temp.acc
                AND map.version = map_temp.version
                AND map.txStart = map_temp.txStart
                AND (
                  map.cdsStart != map_temp.cdsStart
                  OR map.cdsEnd != map_temp.cdsEnd
                );
        """, None

        self.query(statement)
    #backup_cdsUpdates

    def count_cdsUpdates(self) :
        """
        Count the number of mapping entries that have changed without an
        increment in the version number. This function can only be called
        after backup_cdsUpdates() has been executed and before
        merge_cdsUpdates has been executed.

        SQL tables from dbNames:
            - map_cdsBackup_temp ; Entries that wre updated without an
                                   increment of the version number.

        @return: The number of mapping entries that have changed without an
        increment in the version number
        @rtype: integer
        """

        statement = """
            SELECT COUNT(*)
              FROM map_cdsBackup_temp;
        """, None

        return int(self.query(statement)[0][0])
    #count_cdsUpdates

    def merge_cdsUpdates(self) :
        """
        Merge the mapping entries that have changed without an increment in
        the version number with a table that contains backups of these
        entries.

        SQL tables from dbNames (altered):
            - map_cdsBackup      ; Extended with the entries in
                                   map_cdsBackup_temp.
            - map_cdsBackup_temp ; Dropped.
        """

        # The statements in this function may be combined when MYSQL_BUG is
        # solved.

        statement = """
            INSERT INTO map_cdsBackup
              SELECT *
                FROM map_cdsBackup_temp;
        """, None
        self.query(statement)
        statement = """
            DROP TABLE map_cdsBackup_temp;
        """, None

        self.query(statement)
    #merge_cdsUpdates

    def merge_Update(self) :
        """
        Merge the new mapping data with the old ones.

        SQL tables from dbNames (altered):
            - map_new  ; Created and filled with the merge of map_temp and map.
                         Dropped after use.
            - map_temp ; Merged with map to form map_new. Dropped after use.
            - map      ; Overwritten with the merged info in map_new.
        """

        # The statements in this function may be combined when MYSQL_BUG is
        # solved.

        statement = """
            CREATE TABLE map_new
              SELECT *
                FROM map_temp
              UNION
              SELECT *
                FROM map
                WHERE NOT EXISTS (
                  SELECT *
                    FROM map_temp
                    WHERE map.acc = map_temp.acc
                    AND map.version = map_temp.version
                    AND map.txStart = map_temp.txStart
                );
        """, None
        self.query(statement)
        statement = """
            DROP TABLE map;
        """, None
        self.query(statement)
        statement = """
            CREATE TABLE map
              SELECT *
                FROM map_new;
        """, None
        self.query(statement)
        statement = """
            DROP TABLE map_new;
        """, None
        self.query(statement)
        statement = """
            DROP TABLE map_temp;
        """, None

        self.query(statement)
    #merge_Update
#Update

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
