#!/usr/bin/python

class Db() :
    """
        Log in to a database and keep it open for queries.

        Private variables:
            __cursor ; Interface to the database.

        Special methods:
            __init__(config) ; Do the login.
            __del__()        ; Close the handle to the database.

        Private methods:
            __query(statement) ; General query function.

        Public methods:
            get_protAcc(mrnaAcc)             ; Query the database for a protein 
                                               ID.
            get_NM_info(mrnaAcc)             ; Retrieve various data for an NM 
                                               number.
            get_NM_version(mrnaAcc)          ; Get the version number of an 
                                               accession number.
            get_Transcripts(chrom, position) ; Get a list of transcripts, given
                                               a chromosome and a position on 
                                               that chromosome.
    """

    def __init__(self, config) :
        """
            Log in to the database. The username and the name of the 
            database are given in the configuration file.

            Arguments:
                config ; The configuration object.

            Private variables (altered):
                __cursor ; The interface to the database.
        """

        import MySQLdb

        db = MySQLdb.connect(user = config.MySQLuser, db = config.dbName)
        self.__cursor = db.cursor();
    #__init__

    def __del__(self) :
        """
            Close the handle to the database.
        """

        self.__cursor.close()
    #__del__

    def __query(self, statement) :
        """
            Query the database.

            Arguments:
                statement ; The statement that is to be queried.

            Returns:
                list ; The result of the query.

            Private variables:
                __cursor ; Interface to the database.
        """

        self.__cursor.execute(statement)
        result = self.__cursor.fetchall()

        return result
    #__query


    def get_protAcc(self, mrnaAcc) :
        """
            Query the database for a protein ID given an mRNA ID.

            Arguments:
                mrnaAcc ; The ID of an mRNA.

            Returns:
                string ; The protein ID .
        """

        statement = """
            SELECT protAcc
            FROM refLink
            WHERE mrnaAcc = "%s";
        """ % mrnaAcc

        return self.__query(statement)[0][0]
    #get_protAcc

    def get_NM_info(self, mrnaAcc) :
        """
            Retrieve various data for an NM number.

            Arguments:
                mrnaAcc ; The ID of an mRNA.

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
            FROM refGene
            WHERE name = "%s";
        """ % mrnaAcc

        return self.__query(statement)[0]
    #get_NM_info

    def get_NM_version(self, mrnaAcc) :
        """
            Get the version number of an accession number.

            Arguments:
                mrnaAcc ; The ID of an mRNA.

            Returns:
                integer ; The version number.
        """

        statement = """
            SELECT version
            FROM gbStatus
            WHERE acc = "%s";
        """ % mrnaAcc

        ret = self.__query(statement)
        if ret :
            return int(ret[0][0])
        return 0
    #get_NM_version

    def get_Transcripts(self, chrom, position) :
        """
            Get a list of transcripts, given a chromosome and a position on
            that chromosome.

            Arguments:
                chrom    ; The chromosome (coded as "chr1", ..., "chrY").
                position ; The position relative to the start of the chromosome.

            Returns:
                list ; All accession numbers that are hit.
        """

        statement = """
            SELECT name
            FROM refGene
            WHERE chrom = "%s" AND
                  txStart <= "%i" AND
                  txEnd >= "%i";
        """ % (chrom, position, position)

        ret = [] # Convert the results to a normal list.
        for i in self.__query(statement) :
            ret.append(i[0] + '.' + str(self.get_NM_version(i[0])))
        return ret
    #get_Transcripts

    def get_GeneName(self, mrnaAcc) :
        """
        """

        statement = """
            SELECT name2
            FROM refGene
            WHERE name = "%s";
        """ % mrnaAcc

        return self.__query(statement)[0][0]
    #get_GeneName
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
