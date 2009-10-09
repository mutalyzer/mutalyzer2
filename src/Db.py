#!/usr/bin/python

class Db() :
    """
        Log in to a database and keep it open for queries.

        Private variables:
            __cursor ; Interface to the database.

        Special methods:
            __init__(config) ; Do the login.

        Public methods:
            get_protAcc(mrnaAcc) ; Query the database for a protein ID.
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

    def get_protAcc(self, mrnaAcc) :
        """
            Query the database for a protein ID given an mRNA ID.

            Arguments:
                mrnaAcc ; The ID of an mRNA.

            Returns:
                string ; The protein ID .

            Private variables:
                __cursor ; The interface to the database.
        """

        statement = """
            SELECT protAcc
            FROM refLink
            WHERE mrnaAcc = "%s";
        """ % mrnaAcc

        self.__cursor.execute(statement)
        result = self.__cursor.fetchall()

        return result[0][0]
    #get_NP
#Db
