#!/usr/bin/python

class Db() :
    """
        Log in to a database and keep it open for queries.

        Private variables:
            __cursor ; Interface to the database.

        Special methods:
            __init__(config) ; Do the login.

        Private methods:
            __query(statement) ; General query function.

        Public methods:
            get_protAcc(mrnaAcc)    ; Query the database for a protein ID.
            get_NM_info(mrnaAcc)    ; Retrieve various data for an NM number.
            get_NM_version(mrnaAcc) ; Get the version number of an accession 
                                      number.
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

        return self.__query(statement)[0][0]
    #get_NM_version
#Db
