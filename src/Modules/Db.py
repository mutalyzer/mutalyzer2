#!/usr/bin/python

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
            get_protAcc(mrnaAcc)      ; Query the database for a protein ID.
            get_NM_info(mrnaAcc)      ; Retrieve various data for an NM number.
            get_NM_version(mrnaAcc)   ; Get the version number of an accession 
                                        number.
            get_Transcripts(chrom,    ; Get a list of transcripts, given a
                            position,   chromosome and a range. 
                            overlap)
            get_GeneName(mrnaAcc)     ; Get the gene name, given an NM number.
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
    """

    #
    # Note that compound queries are split into single queries because of a bug
    # in MySQLdb. The functions load_Update(),  merge_cdsUpdates() and
    # merge_Update (search for MYSQL_BUG in this file) are affected and may be
    # rewritten when this bug is fixed.
    #

    def __init__(self, config, where) :
        """
            Log in to the database. The username and the name of the 
            database are given in the configuration file.

            Arguments:
                config ; The configuration object.
                where  ; A switch to see which database to use:
                         local  ; Use the database on localhost.
                         remote ; Use the UCSC database.

            Private variables (altered):
                __db       ; The interface to the database.
                __tempfile ; A temporary file used for data exchange between
                             the remote and the local databases.
                __interval ; The time window (in days) in which we retrieve
                             updates from the UCSC.
        """

        import MySQLdb

        if where == "remote" :
            self.__db = MySQLdb.connect(user = config.RemoteMySQLuser, 
                                        db = config.dbName,
                                        host = config.RemoteMySQLhost)
        else :                                 
            self.__db = MySQLdb.connect(user = config.LocalMySQLuser, 
                                        db = config.dbName)
                                 
        self.__tempfile = config.TempFile
        self.__interval = config.UpdateInterval
    #__init__

    def __query(self, statement) :
        """
            Query the database.

            Arguments:
                statement ; The statement that is to be queried.

            Returns:
                list ; The result of the query.

            Private variables:
                __db ; Interface to the database.
        """

        cursor = self.__db.cursor()
        cursor.execute(statement)
        result = cursor.fetchall()
        cursor.close()

        return result
    #__query

    def get_protAcc(self, mrnaAcc) :
        """
            Query the database for a protein ID given an mRNA ID.

            Arguments:
                mrnaAcc ; The ID of an mRNA.

            SQL tables:
                map ; Accumulated mapping info.

            Returns:
                string ; The protein ID .
        """

        statement = """
            SELECT protAcc
              FROM map
              WHERE acc = "%s";
        """ % mrnaAcc

        return self.__query(statement)[0][0]
    #get_protAcc

    def get_NM_info(self, mrnaAcc) :
        """
            Retrieve various data for an NM number.

            Arguments:
                mrnaAcc ; The ID of an mRNA.

            SQL tables:
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
              WHERE acc = "%s";
        """ % mrnaAcc

        return self.__query(statement)[0]
    #get_NM_info

    def get_NM_version(self, mrnaAcc) :
        """
            Get the version number of an accession number.

            Arguments:
                mrnaAcc ; The ID of an mRNA.

            SQL tables:
                map ; Accumulated mapping info.

            Returns:
                integer ; The version number.
        """

        statement = """
            SELECT version
              FROM map
              WHERE acc = "%s";
        """ % mrnaAcc

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
                          0 ; Return all hit transcripts.
                          1 ; Return only the transcripts that completely fall
                              in the range [p1, p2].

            SQL tables:
                map ; Accumulated mapping info.

            Returns:
                list ; All accession numbers that are hit according to the
                       overlap criterium.
        """

        if overlap :
            statement = """
                SELECT acc
                  FROM map
                  WHERE chrom = "%s" AND
                        txStart <= "%i" AND
                        txEnd >= "%i";
            """ % (chrom, p2, p1)
        #if
        else :
            statement = """
                SELECT acc
                  FROM map
                  WHERE chrom = "%s" AND
                        txStart >= "%i" AND
                        txEnd <= "%i";
            """ % (chrom, p1, p2)
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

            SQL tables:
                map ; Accumulated mapping info.

            Returns:
                string ; The gene name.
        """

        statement = """
            SELECT geneName
              FROM map
              WHERE acc = "%s";
        """ % mrnaAcc

        return self.__query(statement)[0][0]
    #get_GeneName

    def get_Update(self) :
        """
            Retrieve all mapping updates from the UCSC within a certain time
            window (defined in the configuration file) and gather the results
            into one mapping table.

            The results will be written to a temporary file (also defined in
            the configuration file) to be imported in the local database with
            the load_Update() function.

            
            Private variables:
                __interval ; The size of the time window (from the config file).
                __tempfile ; The name and location of the temporary file (from
                             the config file). This file is created if it
                             doesn't exist and is overwritten if it does exist.
                             The function load_Update() will remove this file.

            SQL tables:
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
              AND modDate >= DATE_SUB(CURDATE(), INTERVAL %i DAY);
        """ % self.__interval

        handle = open(self.__tempfile, "w")

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

            Private variables:
                __tempfile ; The name and location of the temporary file (from
                             the config file). This file is created by the 
                             get_Update() function. After the local import is 
                             complete, this file will be removed.
            
            SQL tables (altered):
                map_temp ; Created and loaded with data from __tempfile.

            SQL tables:
                map ; Accumulated mapping info.
        """

        import os # os.remove()

        # The statements in this function may be combined when MYSQL_BUG is
        # solved.

        statement = """
            CREATE TABLE map_temp LIKE map;
        """
        self.__query(statement)
        statement = """
            LOAD DATA LOCAL INFILE "%s" 
              INTO TABLE map_temp;
        """ % self.__tempfile

        self.__query(statement)

        os.remove(self.__tempfile)
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
        """

        return int(self.__query(statement)[0][0])
    #count_Updates

    def backup_cdsUpdates(self) :
        """
            Copy all mapping entries where there was an update, but no
            increment in the version number, to a backup table. Note that
            we use acc, version, txStart as the primary key because members
            of a gene family are mapped multiple times.

            SQL tables (altered):
                map_cdsBackup_temp ; Created and filled with entries that
                                     were updated without an increment of the
                                     version number.

            SQL tables:
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
        """

        self.__query(statement)
    #backup_cdsUpdates

    def count_cdsUpdates(self) :
        """
            Count the number of mapping entries that have changed without an
            increment in the version number. This function can only be called
            after backup_cdsUpdates() has been executed and before
            merge_cdsUpdates has been executed.

            SQL tables:
                map_cdsBackup_temp ; Entries that wre updated without an 
                                     increment of the version number.

            Returns:
                int ; The number of mapping entries that have changed without
                      an increment in the version number.
        """

        statement = """
            SELECT COUNT(*)
              FROM map_cdsBackup_temp;
        """

        return int(self.__query(statement)[0][0])
    #count_cdsUpdates

    def merge_cdsUpdates(self) :
        """
            Merge the mapping entries that have changed without an increment in
            the version number with a table that contains backups of these 
            entries.

            SQL tables (altered):
                map_cdsBackup      ; Extended with the entries in 
                                     map_cdsBackup_temp.
                map_cdsBackup_temp ; Dropped.
        """

        # The statements in this function may be combined when MYSQL_BUG is
        # solved.

        statement = """
            INSERT INTO map_cdsBackup
              SELECT * FROM map_cdsBackup_temp;
        """
        self.__query(statement)
        statement = """
            DROP TABLE map_cdsBackup_temp;
        """

        self.__query(statement)
    #merge_cdsUpdates

    def merge_Update(self) :
        """
            Merge the new mapping data with the old ones.

            SQL tables (altered):
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
        """
        self.__query(statement)
        statement = """
            DROP TABLE map;
        """
        self.__query(statement)
        statement = """
            CREATE TABLE map
              SELECT * FROM map_new;
        """
        self.__query(statement)
        statement = """
            DROP TABLE map_new;
        """
        self.__query(statement)
        statement = """
            DROP TABLE map_temp;
        """

        self.__query(statement)
    #merge_Update
#Db

#
# Unit test.
#
if __name__ == "__main__" :
    import Config

    # Get the username / db from the config file.
    C = Config.Config()
    D = Db(C)
    del C

    # Do some basic testing (will crash if MySQL is not set up properly.
    D.get_protAcc("NM_002001")
    D.get_NM_info("NM_002001")
    D.get_NM_version("NM_002001")
    D.get_Transcripts("chr1", 159272155)
    D.get_GeneName("NM_002001")
    del D
#if
