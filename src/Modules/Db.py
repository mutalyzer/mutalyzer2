#!/usr/bin/python

"""
    Module for database access.
    The Db class is a superclass of the rest of the classes and should not be
    used as such. The superclass mainly consists of a wrapper for SQL
    statements.


    Public classes:
        Db      ; Log in to a database and keep it open for queries.
        Mapping ; Mapping of transcripts and genes.
        Remote  ; Retrieving updates for the mapping databases.
        Update  ; Updating the mapping databases.
        Cache   ; Cache administration.
        Batch   ; Batch checker.
"""

import MySQLdb # connect(), escape_string()
import types   # TupleType
import time    # strftime()
import os      # os.remove()

from Modules import Misc # ID()

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
            __db ; Interface to the database.

        Special methods:
            __init__(dbName, mySqlUser, mySqlHost) ; Do the login.

        Public methods:
            query(statement) ; General query function.
    """

    def __init__(self, dbName, mySqlUser, mySqlHost) :
        """
            Log in to the database. 

            Arguments:
                dbName    ; The name of the database to use.
                mySqlUser ; User name for the database.
                mySqlHost ; Host name for the database.

            Private variables (altered):
                __db       ; The interface to the database.
        """

        self.__db = MySQLdb.connect(user = mySqlUser, db = dbName,
                                    host = mySqlHost)
    #__init__

    def query(self, statement) :
        """
            Query the database.

            Arguments:
                statement ; The statement that is to be queried, consists of
                            a tuple: (string, (args)).

            Returns:
                list ; The result of the query.

            Private variables:
                __db ; Interface to the database.
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
                    if type(i) == types.StringType :
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
            __init__(build, config) ; Initialise the class.

        Public methods:
            get_protAcc(mrnaAcc)      ; Query the database for a protein ID.
            get_NM_info(mrnaAcc)      ; Retrieve various data for an NM number.
            get_NM_version(mrnaAcc)   ; Get the version number of an accession 
                                        number.
            get_Transcripts(chrom,    ; Get a list of transcripts, given a
                            position,   chromosome and a range. 
                            overlap)
            get_GeneName(mrnaAcc)     ; Get the gene name, given an NM number.
            isChrom(name)             ; Check whether we know this name to be
                                        a chromosome name.

        Inherited methods from Db:
            query(statement) ; General query function.

        SQL tables from dbNames:
            map ; Accumulated mapping info.
    """

    def __init__(self, build, config) :
        """
            Initialise the Db parent class. Use the local database for a
            certain build.

            Arguments:
                build  ; The version of the mapping database.
                config ; Configuration variables.
        """

        Db.__init__(self, build, config.LocalMySQLuser, config.LocalMySQLhost)
    #__init__

    def get_protAcc(self, mrnaAcc) :
        """
            Query the database for a protein ID given an mRNA ID.

            Arguments:
                mrnaAcc ; The ID of an mRNA.

            SQL tables from dbNames:
                map ; Accumulated mapping info.

            Returns:
                string ; The protein ID .
        """

        statement = """
            SELECT protAcc
              FROM map
              WHERE acc = %s;
        """, mrnaAcc

        return self.query(statement)[0][0]
    #get_protAcc

    def get_NM_info(self, mrnaAcc) :
        """
            Retrieve various data for an NM number.

            Arguments:
                mrnaAcc ; The ID of an mRNA.

            SQL tables from dbNames:
                map ; Accumulated mapping info.

            Returns:
                list:
                    exonStarts ; List of exon start sites.
                    exonEnds   ; List of exon end sites.
                    cdsStart   ; Position of the start codon.
                    cdsEnd     ; Position of the end codon.
                    strand     ; The orientation of the gene (+ = forward,
                                                              - = reverse).
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

            Arguments:
                mrnaAcc ; The ID of an mRNA.

            SQL tables from dbNames:
                map ; Accumulated mapping info.

            Returns:
                integer ; The version number.
        """

        statement = """
            SELECT version
              FROM map
              WHERE acc = %s;
        """, mrnaAcc

        ret = self.query(statement)
        if ret :
            return int(ret[0][0])
        return 0
    #get_NM_version

    def get_Transcripts(self, chrom, p1, p2, overlap) :
        """
            Get a list of transcripts, given a chromosome and a range. If 
            all transcripts that are hit should be returned, set overlap to 1,
            if only the transcripts that completely reside within a range
            should be returned, set overlap to 0.

            Arguments:
                chrom   ; The chromosome (coded as "chr1", ..., "chrY").
                p1      ; The position relative to the start of the chromosome.
                p2      ; The position relative to the start of the chromosome.
                overlap ; Specify the behaviour of the selection:
                          0 ; Return only the transcripts that completely fall
                              in the range [p1, p2].
                          1 ; Return all hit transcripts.

            SQL tables from dbNames:
                map ; Accumulated mapping info.

            Returns:
                list ; All accession numbers that are hit according to the
                       overlap criterium.
        """

        if overlap :
            statement = """
                SELECT acc
                  FROM map
                  WHERE chrom = %s AND
                        txStart <= "%s" AND
                        txEnd >= "%s";
            """, (chrom, p2, p1)
        #if
        else :
            statement = """
                SELECT acc
                  FROM map
                  WHERE chrom = %s AND
                        txStart >= "%s" AND
                        txEnd <= "%s";
            """, (chrom, p1, p2)
        #else

        ret = [] # Convert the results to a normal list.
        for i in self.query(statement) :
            ret.append(i[0] + '.' + str(self.get_NM_version(i[0])))
        return ret
    #get_Transcripts

    def get_GeneName(self, mrnaAcc) :
        """
            Get the name of a gene, given a transcript identifier (NM number).

            Arguments:
                mrnaAcc ; The ID of an mRNA.

            SQL tables from dbNames:
                map ; Accumulated mapping info.

            Returns:
                string ; The gene name.
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

            Arguments:
                name ; The name to be tested.

            SQL tables from dbNames:
                map ; Accumulated mapping info.

            Returns:
                boolean ; True if the name is found to be a chromosome name,
                          False otherwise.
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

            Arguments:
                accNo ; The accession number of a chromosome.

            SQL tables from dbNames:
                ChrName ; Assembly release notes.

            Returns:
                string ; The name of a chromosome.
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

            Arguments:
                name ; The name of a chromosome.

            SQL tables from dbNames:
                ChrName ; Assembly release notes.

            Returns:
                string ; The accession number of a chromosome.
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
            
            Arguments:
                acc ; The NM accession number (version NOT included)

            SQL tables from dbNames:
                map ; .

            Returns:
                string  ; The chromosome name (e.g. chr1)
            
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
            __init__(config) ; Initialise the class.
        
        Public methods:
            get_Update()        ; Retrieve new mapping info from the UCSC.

        Inherited methods from Db:
            query(statement) ; General query function.

        SQL tables from dbNames:
            gbStatus ; acc -> version mapping (NM to NM + version), 
                       type, modDate
            refGene  ; name -> geneName mapping (NM to gene name), 
                       txStart, txEnd, cdsStart, cdsEnd, exonStarts, 
                       exonEnds, chrom, strand.
            refLink  ; mrnaAcc -> protAcc mapping (NM to NP).
    """

    def __init__(self, build, config) :
        """
            Initialise the Db parent class. Use the remote database for a 
            certain build.

            Arguments:
                build  ; The version of the mapping database.
                config ; Configuration variables.

            Private variables (altered):
                __config ; Configuration variables.
        """

        self.__config = config
        Db.__init__(self, build, config.RemoteMySQLuser, config.RemoteMySQLhost)
    #__init__

    def get_Update(self) :
        """
            Retrieve all mapping updates from the UCSC within a certain time
            window (defined in the configuration file) and gather the results
            into one mapping table.

            The results will be written to a temporary file (also defined in
            the configuration file) to be imported in the local database with
            the load_Update() function.

            
            SQL tables from dbNames:
                gbStatus ; acc -> version mapping (NM to NM + version), 
                           type, modDate
                refGene  ; name -> geneName mapping (NM to gene name), 
                           txStart, txEnd, cdsStart, cdsEnd, exonStarts, 
                           exonEnds, chrom, strand.
                refLink  ; mrnaAcc -> protAcc mapping (NM to NP).
        """

        statement = """
            SELECT DISTINCT acc, version, txStart, txEnd, cdsStart, cdsEnd, 
                            exonStarts, exonEnds, name2 AS geneName, chrom,
                            strand, protAcc
              FROM gbStatus, refGene, refLink
              WHERE type = "mRNA"
              AND refGene.name = acc
              AND acc = mrnaAcc
              AND modDate >= DATE_SUB(CURDATE(), INTERVAL %s DAY);
        """, self.__config.UpdateInterval

        handle = open(self.__config.TempFile, "w")

        # Convert the results to a tab delimited file.
        for i in self.query(statement) :
            for j in i :
                handle.write(str(j) + chr(0x09))  # 0x09 is a TAB.
            handle.write('\n')
        #for

        handle.close()
    #get_Update
#Remote    

class Update(Db) :
    """
        Database functions for updating the mapping databases.

        Public methods:
            load_Update()       ; Load new mapping info into the local database.
            count_Updates()     ; Count the number of entries in the new
                                  mapping info table.
            backup_cdsUpdates() ; Make a backup of updates that overwrite the
                                  old mapping info.
            count_cdsUpdates()  ; Count the number of updates that overwrite
                                  the old mapping info.
            merge_cdsUpdates()  ; Merge the backup of old mapping info with the
                                  other old info.
            merge_Update()      ; Merge the new mapping info from the UCSC with
                                  what we already have.

        Inherited methods from Db:
            query(statement) ; General query function.

        SQL tables from dbNames:
            map                ; Accumulated mapping info.
            map_temp           ; Newly found data.
            map_new            ; Merge of map_temp and map.
            map_cdsBackup_temp ; Entries that were updated without an increment
                                 of the version number.
            map_cdsBackup      ; Merge of map_cdsBackup_temp and itself.
    """

    def __init__(self, build, config) :
        """
            Initialise the Db parent class. Use the remote database for a 
            certain build.

            Arguments:
                build  ; The version of the mapping database.
                config ; Configuration variables.

            Private variables (altered):
                __config ; Configuration variables.
        """

        self.__config = config
        Db.__init__(self, build, config.LocalMySQLuser, config.LocalMySQLhost)
    #__init__

    def load_Update(self) :
        """
            Load the updates from the temporary file (defined in the
            configuration file) created by the get_Update() function and import
            it in the local database.

            SQL tables from dbNames (altered):
                map_temp ; Created and loaded with data from TempFile.

            SQL tables from dbNames:
                map ; Accumulated mapping info.
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
        """, self.__config.TempFile

        self.query(statement)

        os.remove(self.__config.TempFile)
    #load_Update

    def count_Updates(self) :
        """
            Count the number of updates. This function will only work if it
            is preceeded by the load_Update() function. Otherwise the map_temp
            table may not exist. This function can not be used after the 
            merge_Update() function has been executed, since it drops the
            map_temp table.

            Returns:
                int ; The number of entries in the table of updated mapping 
                      info.
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
                map_cdsBackup_temp ; Created and filled with entries that
                                     were updated without an increment of the
                                     version number.

            SQL tables from dbNames:
                map      ; Accumulated mapping info.
                map_temp ; Freshly downloaded mapping info.
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
                map_cdsBackup_temp ; Entries that wre updated without an 
                                     increment of the version number.

            Returns:
                int ; The number of mapping entries that have changed without
                      an increment in the version number.
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
                map_cdsBackup      ; Extended with the entries in 
                                     map_cdsBackup_temp.
                map_cdsBackup_temp ; Dropped.
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
                map_new  ; Created and filled with the merge of map_temp and
                           map. Dropped after use.
                map_temp ; Merged with map to form map_new. Dropped after use.
                map      ; Overwritten with the merged info in map_new.
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
            __init__(config) ; Initialise the class.
        
        Public methods:
            insertGB(accNo, GI,       ; Insert info about a GenBank record.
                     fileHash,  
                     ChrAccVer,      
                     ChrStart, 
                     ChrStop, 
                     orientation, 
                     url)
            updateHash(accNo,         ; Update the hash of an accession number.
                       fileHash)
            getGBFromLoc(ChrAccVer,   ; Get the accession number from slicing 
                         ChrStart,      information.
                         ChrStop, 
                         orientation)
            getGBFromHash(fileHash)   ; Get the accession number from its hash.
            getGBFromGI(GI)           ; Get the accession number from its GI 
                                        number.
            getLoc(accNo)             ; Get the slicing information of an 
                                        accession number.
            getHash(accNo)            ; Get the hash of a GenBank record.
            getUrl(accNo)             ; Get the URL of an accession number.

        Inherited methods from Db:
            query(statement) ; General query function.

        SQL tables from internalDb:
            GBInfo ; Information about cached and uploaded GenBank files.
    """

    def __init__(self, config) :
        """
            Initialise the Db parent class. Use the internalDb.

            Arguments:
                config ; Configuration variables.
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
            - If the record is a chromosome slice, then the ChrAccVer, 
              ChrStart, ChrStop and orientation variables should be specified.
            - If the record is downloaded from the internet, the url should
              be provided.
            - If all fields except the mandatory ones are empty, the record
              is assumed to be uploaded.

            Arguments:
                accNo       ; The name associated with this record.
                GI          ; The GI number (if available).
                fileHash    ; The hash of the content of the record.
                ChrAccVer   ; The accession number of the chromosome (if 
                              available).
                ChrStart    ; The start of the record in chromosomal 
                              coordinates (if available).
                ChrStop     ; The end of the record in chromosomal coordinates 
                              (if available).
                orientation ; The orientation of the record relative to the 
                              chromosome (if available) (1 = forward, 
                              2 = reverse complement).
                url         ; The originating URL (if available).

            SQL tables from internalDb (altered):
                GBInfo ; Information about cached and uploaded GenBank files.
        """                 

        statement = """
            INSERT INTO GBInfo 
              VALUES (%s, %s, %s, %s, %s, %s, %s, %s);
        """, (accNo, GI, fileHash, ChrAccVer, ChrStart, ChrStop, orientation, 
              url)

        self.query(statement)
    #insertGB

    def updateHash(self, accNo, fileHash) :
        """
            Update the hash of an accession number.

            Arguments:
                accNo     ; The accession number of a GenBank record.
                fileHash  ; The hash of a GenBank record.

            SQL tables from internalDb (altered):
                GBInfo ; Information about cached and uploaded GenBank files.
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

            Arguments:
                ChrAccVer   ; The accession number of the chromosome.
                ChrStart    ; Start position of the slice.
                ChrStop     ; End position of the slice.
                orientation ; Orientation of the slice:
                              1  ; Forward.
                              2 ; Reverse complement.

            SQL tables from internalDb:
                GBInfo ; Information about cached and uploaded GenBank files.

            Returns:
                string ; The accession number.
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

            Arguments:
                fileHash ; The hash of a GenBank record.

            SQL tables from internalDb:
                GBInfo ; Information about cached and uploaded GenBank files.

            Returns:
                string ; The accession number.
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

            Arguments:
                GI ; The GI number of a GenBank record.

            SQL tables from internalDb:
                GBInfo ; Information about cached and uploaded GenBank files.

            Returns:
                string ; The accession number.
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

    def getLoc(self, accNo) :
        """
            Get the slicing information of an accession number, typically this
            only affects UD numbers.

            Arguments:
                accNo ; The accession number of a genbank record.

            SQL tables from internalDb:
                GBInfo ; Information about cached and uploaded GenBank files.

            Returns:
                list ; The slicing information:
                       ChrAccVer   ; Accession number of the chromosome.
                       ChrStart    ; Start position of the slice.
                       ChrStop     ; End position of the slice.
                       orientation ; Orientation of the slice (1 = forward, 
                                     2 = reverse complement).
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

            Arguments:
                accNo ; The accession number of a genbank record.

            SQL tables from internalDb:
                GBInfo ; Information about cached and uploaded GenBank files.

            Returns:
                string ; The hash of the GenBank record.
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

            Arguments:
                accNo ; The accession number of a genbank record.

            SQL tables from internalDb:
                GBInfo ; Information about cached and uploaded GenBank files.

            Returns:
                string ; The URL of the GenBank record.
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

            Arguments:
                accNo ; The accession number.        

            SQL tables from internalDb:
                GBInfo ; Information about cached and uploaded GenBank files.
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
#Cache    

class Batch(Db) :
    """
        Database functions for the batch checker.

        Special methods:
            __init__(config) ; Initialise the class.

        Public methods:
            isJobListEmpty()     ; See if there are active jobs.
            addJob(outputFilter, ; Add a job and give it a unique ID.
                   email, 
                   fromHost) 
            getJobs()            ; Get a list of active jobs.
            removeJob(jobID)     ; Remove a job and return information about 
                                   the job submitter.
            addToQueue(jobID,    ; Add a request belonging to a certain job to
                       accNo,      the queue.
                       gene, 
                       variant)
            getFromQueue(jobID)  ; Get a request belonging to a certain job 
                                   from the queue.
        
        Inherited methods from Db:
            query(statement) ; General query function.

        SQL tables from internalDb:
            BatchJob   ; Job information.
            BatchQueue ; Requests.
    """

    def __init__(self, config) :
        """
            Initialise the Db parent class. Use the internalDb.

            Arguments:
                config ; Configuration variables.
        """

        Db.__init__(self, config.internalDb, config.LocalMySQLuser, 
                    config.LocalMySQLhost)
    #__init__

    def isJobListEmpty(self) :
        """
            See if there are active jobs.

            SQL tables from internalDb:
                BatchJob ; Job information.

            Returns:
                boolean ; False if there are active jobs, True otherwise.
        """

        statement = """
            SELECT COUNT(*)
              FROM BatchJob;
        """, None

        if int(self.query(statement)[0][0]) :
            return False
        return True
    #isJobListEmpty

    def addJob(self, outputFilter, email, fromHost) :
        """
            Add a job and give it a unique ID.

            Arguments:
                outputFilter ; Output settings for all requests in this job.
                email        ; Contact information of the submitter.

            SQL tables from internalDb (altered):
                BatchJob ; Job information.

            Returns:
                int ; A job ID.
        """

        M = Misc.Misc()
        jobID = M.ID()
        del M
        statement = """
            INSERT INTO BatchJob
              VALUES (%s, %s, %s, %s);
        """, (jobID, outputFilter, email, fromHost)

        self.query(statement)
        return jobID
    #addJob

    def getJobs(self) :
        """
            Get a list of active jobs.

            SQL tables from internalDb:
                BatchJob ; Job information.

            Returns:
                list ; List of job IDs.
        """
        
        statement = """
            SELECT JobID
              FROM BatchJob;
        """, None

        ret = []
        for i in self.query(statement) :
            ret.append(i[0])
        return ret
    #getJobs

    def removeJob(self, jobID) :
        """
            Remove a job (because the queue for this job is empty) and return
            information needed to alert the job submitter.
            
            Arguments:
                jobID   ; Identifier of a job.

            SQL tables from internalDb (altered):
                BatchJob ; Job information.

            Returns:
                triple ; Data for the job submitter.
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

    def addToQueue(self, jobID, accNo, gene, variant) :
        """
            Add a request belonging to a certain job to the queue.

            Arguments:
                jobID   ; Identifier of a job.
                accNo   ; The accession number of a request.
                gene    ; The gene and transcript variant information.
                variant ; The variant.

            SQL tables from internalDb (altered):
                BatchQueue ; Requests.
        """

        # The first value (QueueID) will be auto increased by MySQL.
        statement = """
            INSERT INTO BatchQueue
              VALUES (%s, %s, %s, %s, %s);
        """, (None, jobID, accNo, gene, variant)

        self.query(statement)
    #addToQueue

    def getFromQueue(self, jobID) :
        """
            Get a request belonging to a certain job from the queue. If a
            request is found, remove it from the queue and return it. Otherwise
            return nothing.

            Arguments:
                jobID ; Identifier of a job.

            SQL tables from internalDb (altered):
                BatchQueue ; Requests.

            Returns:
                triple:
                    accNo   ; The accession number of a request.
                    gene    ; The gene and transcript variant information.
                    variant ; The variant.
        """

        statement = """
            SELECT QueueID, AccNo, Gene, Variant
              FROM BatchQueue
              WHERE JobID = %s
              ORDER BY QueueID
              LIMIT 1;
        """, jobID

        results = self.query(statement)
        if results :
            jobID, accNo, gene, variant = results[0]
        else :
            return None

        # We have found a request, so remove it from the queue.
        statement = """
            DELETE 
              FROM BatchQueue
              WHERE QueueID = %s;
        """, jobID

        self.query(statement)
        return accNo, gene, variant
    #getFromQueue
#Batch    

#
# Unit test.
#
if __name__ == "__main__" :
    pass
#if
