#!/usr/bin/python

import MySQLdb # connect(), escape_string()
import types   # TupleType
import time

#from Output import Output
import os # os.remove()

from Modules import Misc

class Db() :
    """
        Log in to a database and keep it open for queries.

        Private variables:
            __db ; Interface to the database.

        Special methods:
            __init__(config, where) ; Do the login.

        Private methods:
            __query(statement) ; General query function.

        Public methods:
            # For mapping.
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

            # For updating mapping information
            get_Update()              ; Retrieve new mapping info from the UCSC.
            load_Update()             ; Load new mapping info into the local 
                                        database.
            count_Updates()           ; Count the number of entries in the new 
                                        mapping info table.
            backup_cdsUpdates()       ; Make a backup of updates that overwrite
                                        the old mapping info.
            count_cdsUpdates()        ; Count the number of updates that
                                        overwrite the old mapping info.
            merge_cdsUpdates()        ; Merge the backup of old mapping info 
                                        with the other old info.
            merge_Update()            ; Merge the new mapping info from the 
                                        UCSC with what we already have.

            # For cache administration.
            insertGB(AccNo, GI, md5,  ; Insert info about a GenBank record.
                     ChrAccVer,      
                     ChrStart, 
                     ChrStop, 
                     orientation, 
                     url)
            updateHash(AccNo, md5)    ; Update the hash of an accession number.
            getGBFromLoc(ChrAccVer,   ; Get the accession number from slicing 
                         ChrStart,      information.
                         ChrStop, 
                         orientation)
            getGBFromHash(md5)        ; Get the accession number from its hash.
            getGBFromGI(GI)           ; Get the accession number from its GI 
                                        number.
            getLoc(AccNo)             ; Get the slicing information of an 
                                        accession number.
            getHash(AccNo)            ; Get the hash of a GenBank record.
            getUrl(AccNo)             ; Get the URL of an accession number.

        Inherited from Output.Config:
            internalDb      ; Name of the internal database.
            RemoteMySQLuser ; MySQL username for the UCSC database.
            RemoteMySQLhost ; Host name for the UCSC database.
            LocalMySQLuser  ; MySQL username for the local databases.
            UpdateInterval  ; The size of the time window.
            TempFile        ; The name and location of the temporary file. This
                              file is created if it doesn't exist and is
                              overwritten if it does exist. The function
                              load_Update() will remove this file.
    """

    #
    # Note that compound queries are split into single queries because of a bug
    # in MySQLdb. The functions load_Update(),  merge_cdsUpdates() and
    # merge_Update (search for MYSQL_BUG in this file) are affected and may be
    # rewritten when this bug is fixed.
    #

    def __init__(self, where, dbName, config) :
        """
            Log in to the database. The username and the name of the 
            database are given in the configuration file.

            Arguments:
                where  ; A switch to see which database to use:
                         local  ; Use the database on localhost.
                         remote ; Use the UCSC database.
                dbName ; The name of the database to use (hg18 or hg19).

            Private variables (altered):
                __db       ; The interface to the database.

            Inherited variables from Output.Config:
                internalDb      ; Name of the internal database.
                RemoteMySQLuser ; MySQL username for the UCSC database.
                RemoteMySQLhost ; Host name for the UCSC database.
                LocalMySQLuser  ; MySQL username for the local databases.
        """

        #Output.__init__(self, __file__)
        self.__config = config

        self.opened = False
        if dbName in self.__config.dbNames or dbName == self.__config.internalDb :
            if where == "remote" :
                self.__db = MySQLdb.connect(user = self.__config.RemoteMySQLuser, 
                                            db = dbName,
                                            host = self.__config.RemoteMySQLhost)
            else :                                 
                self.__db = MySQLdb.connect(user = self.__config.LocalMySQLuser, 
                                            db = dbName)
            self.opened = True                                                
        #if                                            
    #__init__

    def __query(self, statement) :
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
                    escaped_args.append(MySQLdb.escape_string(str(i)))
                else :
                    escaped_args.append(None)
        #if                    
            
        # And do the query.
        cursor = self.__db.cursor()
        cursor.execute(statement[0], tuple(escaped_args))
        result = cursor.fetchall()
        cursor.close()

        return result
    #__query

    #
    # These methods are used for mapping.
    #

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

        return self.__query(statement)[0][0]
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

        return self.__query(statement)[0]
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

        ret = self.__query(statement)
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
        for i in self.__query(statement) :
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

        return self.__query(statement)[0][0]
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

        if int(self.__query(statement)[0][0]) > 0 :
            return True
        return False
    #isChrom

    #
    # These methods are used for updating the mapping information.
    #

    def get_Update(self) :
        """
            Retrieve all mapping updates from the UCSC within a certain time
            window (defined in the configuration file) and gather the results
            into one mapping table.

            The results will be written to a temporary file (also defined in
            the configuration file) to be imported in the local database with
            the load_Update() function.

            
            Inherited variables from Output.Config:
                UpdateInterval ; The size of the time window.
                TempFile       ; The name and location of the temporary file. 
                                 This file is created if it doesn't exist and
                                 is overwritten if it does exist. The function
                                 load_Update() will remove this file.

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
        for i in self.__query(statement) :
            for j in i :
                handle.write(str(j) + chr(0x09))  # 0x09 is a TAB.
            handle.write('\n')
        #for

        handle.close()
    #get_Update

    def load_Update(self) :
        """
            Load the updates from the temporary file (defined in the
            configuration file) created by the get_Update() function and import
            it in the local database.

            Inherited variables from Config:
                TempFile ; The name and location of the temporary file. This 
                           file is created by the get_Update() function. After
                           the local import is complete, this file will be
                           removed.
            
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
        self.__query(statement)
        statement = """
            LOAD DATA LOCAL INFILE %s 
              INTO TABLE map_temp;
        """, self.__config.TempFile

        self.__query(statement)

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

        return int(self.__query(statement)[0][0])
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

        self.__query(statement)
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

        return int(self.__query(statement)[0][0])
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
        self.__query(statement)
        statement = """
            DROP TABLE map_cdsBackup_temp;
        """, None

        self.__query(statement)
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
        self.__query(statement)
        statement = """
            DROP TABLE map;
        """, None
        self.__query(statement)
        statement = """
            CREATE TABLE map
              SELECT * 
                FROM map_new;
        """, None
        self.__query(statement)
        statement = """
            DROP TABLE map_new;
        """, None
        self.__query(statement)
        statement = """
            DROP TABLE map_temp;
        """, None

        self.__query(statement)
    #merge_Update

    #
    # These methods are used for cache administration.
    #

    def insertGB(self, AccNo, GI, md5, ChrAccVer, ChrStart, 
                 ChrStop, orientation, url) :
        """                 
            Insert information about a GenBank record in the internal database.

            The AccNo and md5 arguments are mandatory.
            - If the record is a normal RefSeq, then the GI number should be
              provided.
            - If the record is a chromosome slice, then the ChrAccVer, 
              ChrStart, ChrStop and orientation variables should be specified.
            - If the record is downloaded from the internet, the url should
              be provided.
            - If all fields except the mandatory ones are empty, the record
              is assumed to be uploaded.

            Arguments:
                AccNo       ; The name associated with this record.
                GI          ; The GI number (if available).
                md5         ; The md5sum of the content of the record.
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
        """, (AccNo, GI, md5, ChrAccVer, ChrStart, ChrStop, orientation, url)

        self.__query(statement)
    #insertGB

    def updateHash(self, AccNo, md5) :
        """
            Update the hash of an accession number.

            Arguments:
                AccNo ; The accession number of a GenBank record.
                hash  ; The hash of a GenBank record.

            SQL tables from internalDb (altered):
                GBInfo ; Information about cached and uploaded GenBank files.
        """

        statement = """
            UPDATE GBInfo
              SET hash = %s
              WHERE AccNo = %s;
        """, (md5, AccNo)

        self.__query(statement)
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

        ret = self.__query(statement)
        if ret :
            return ret[0][0]
        return None
    #getGBFromLoc

    def getGBFromHash(self, md5) :
        """
            Get the accession number from its hash.

            Arguments:
                hash ; The hash of a GenBank record.

            SQL tables from internalDb:
                GBInfo ; Information about cached and uploaded GenBank files.

            Returns:
                string ; The accession number.
        """

        statement = """
            SELECT AccNo 
              FROM GBInfo
              WHERE hash = %s;
        """, md5

        ret = self.__query(statement)
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

        ret = self.__query(statement)
        if ret :
            return ret[0][0]
        return None
    #getGBFromGI

    def getLoc(self, AccNo) :
        """
            Get the slicing information of an accession number, typically this
            only affects UD numbers.

            Arguments:
                AccNo ; The accession number of a genbank record.

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
        """, AccNo

        ret = self.__query(statement)
        if ret :
            return list(ret[0])
        return None
    #getLoc

    def getHash(self, AccNo) :
        """
            Get the hash of a GenBank record identified by an accession number.

            Arguments:
                AccNo ; The accession number of a genbank record.

            SQL tables from internalDb:
                GBInfo ; Information about cached and uploaded GenBank files.

            Returns:
                string ; The hash of the GenBank record.
        """

        statement = """
            SELECT hash 
              FROM GBInfo
              WHERE AccNo = %s;
        """, AccNo

        ret = self.__query(statement)
        if ret :
            return ret[0][0]
        return None
    #getHash

    def getUrl(self, AccNo) :
        """
            Get the URL of an accession number, typically this only affects
            uploaded UD numbers.

            Arguments:
                AccNo ; The accession number of a genbank record.

            SQL tables from internalDb:
                GBInfo ; Information about cached and uploaded GenBank files.

            Returns:
                string ; The URL of the GenBank record.
        """
        
        statement = """
            SELECT url 
              FROM GBInfo
              WHERE AccNo = %s;
        """, AccNo

        ret = self.__query(statement)
        if ret :
            return ret[0][0]
        return None
    #getHash

    def getGI(self, AccNo) :
        """
        """
        statement = """
            SELECT GI 
              FROM GBInfo
              WHERE AccNo = %s;
        """, AccNo

        ret = self.__query(statement)
        if ret :
            return ret[0][0]
        return None
    #getGI

    #
    # These methods are for the batch checker.
    #

    def getWatchDogTimer(self) :
        """
        """

        statement = """
            SELECT Value 
              FROM Var
              WHERE Name = "WatchDog";
        """, None

        return int(self.__query(statement)[0][0])
    #getWatchDogTimer

    def setWatchDogTimer(self) :
        """
        """

        statement = """
            UPDATE Var
              SET Value = %s
              WHERE Name = "WatchDog";
        """, time.strftime("%s")

        self.__query(statement)
    #setWatchDogTimer

    def isJobListEmpty(self) :
        """
        """

        statement = """
            SELECT COUNT(*)
              FROM BatchJob;
        """, None

        if int(self.__query(statement)[0][0]) :
            return False
        return True
    #isJobListEmpty

    def addJob(self, output_filter, email) :
        """
        """

        M = Misc.Misc()
        jobID = M.ID()
        del M
        statement = """
            INSERT INTO BatchJob
              VALUES (%s, %s, %s);
        """, (jobID, output_filter, email)

        self.__query(statement)
        return jobID
    #addJob

    def getJobs(self) :
        """
        """
        
        statement = """
            SELECT JobID
              FROM BatchJob;
        """, None

        ret = []
        for i in self.__query(statement) :
            ret.append(i[0])
        return ret
    #getJobs

    def removeJob(self, jobID) :
        """
        """

        statement = """
            SELECT EMail
              FROM BatchJob
              WHERE JobID = %s;
        """, jobID
        eMail = self.__query(statement)[0][0]

        statement = """
            DELETE 
              FROM BatchJob
              WHERE JobID = %s;
        """, jobID

        self.__query(statement)
        return eMail
    #removeJob

    def addToQueue(self, jobID, AccNo, Gene, Variant) :
        """
        """

        statement = """
            INSERT INTO BatchQueue
              VALUES (%s, %s, %s, %s, %s);
        """, (None, jobID, AccNo, Gene, Variant)

        self.__query(statement)
    #addToQueue

    def getFromQueue(self, jobID) :
        """
        """

        statement = """
            SELECT QueueID, AccNo, Gene, Variant
              FROM BatchQueue
              WHERE JobID = %s
              ORDER BY QueueID
              LIMIT 1;
        """, jobID

        results = self.__query(statement)
        if results :
            queueID, accNo, gene, variant = results[0]
        else :
            return None

        statement = """
            DELETE 
              FROM BatchQueue
              WHERE QueueID = %s;
        """, queueID

        self.__query(statement)
        return accNo, gene, variant
    #getFromQueue
#Db

#
# Unit test.
#
if __name__ == "__main__" :
    # Get the username / db from the config file.
    D = Db("local")

    # Do some basic testing (will crash if MySQL is not set up properly.
    D.get_protAcc("NM_002001")
    D.get_NM_info("NM_002001")
    D.get_NM_version("NM_002001")
    D.get_Transcripts("chr1", 159272155, 159272155, 0)
    D.get_GeneName("NM_002001")
    del D
#if
