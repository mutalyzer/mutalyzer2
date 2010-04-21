#!/usr/bin/python

class Config() :
    """
        Read the configuration file and store the data.

        Public variables:
                ; Used by the Retriever module:
            email      ; Email address used for Entrez.
            cache      ; Location of the cache directory.
            cachesize  ; Maximum size of the cache directory.

                ; Used by the Db module:
            dbName          ; The name of the database.
            LocalMySQLuser  ; The user that has read/write access to the local
                              database.
            RemoteMySQLuser ; A user that has read access to the remote UCSC 
                              database.
            RemoteMySQLhost ; The name of the remote UCSC database.
            UpdateInterval  ; The number of days to search back for updates.
            TempFile        ; The temporary file where the downloaded updates
                              are placed.

                ; Used by the Output module:
            log        ; The name and location of the log file.
            datestring ; Format of the prefix for log messages.


        Special Methods:
            __init__ ; Read the configuration file.
    """
    def __init__(self) :
        """
            Initialise the class with variables read from the configuration 
            file. In principle, this is the only place in the code where a
            hard coded constant is used (the name and path to the configuration
            file).

            Public variables (altered):
                          ; Used by the Retriever module:
                email     ; Email address used for Entrez.
                cache     ; Location of the cache directory.
                cachesize ; Maximum size of the cache directory.

                          ; Used by the Db module:
                MySQLuser ; The user that has access to the database.
                dbName    ; The name of the database.
        """
#        import sys                   # argv
#        sys.path.append("/home/gerard/Projects/web_dev/src/")
#        sys.path.append("/home/gerard/Projects/web_dev/")
#        import pdb
#        print sys.path
#        pdb.set_trace()    
        from configobj import ConfigObj # ConfigObj()
        
        # Figure out where this program is located and go two levels up.
        import os
        myPath = os.path.dirname(__file__) + "/../.."
        os.chdir(myPath)
#        print "De huidige directory is: %s" % os.path.dirname(__file__)
#        config = ConfigObj("./mutalyzer.conf")
        # Assumes the configuration file to be two levels up
        config = ConfigObj("./mutalyzer.conf")

        # Set the variables needed by the Retriever module.
        self.email = config["email"]
        self.cache = config["cache"]
        self.cachesize = int(config["cachesize"])

        # Set the variables needed by the Db module.
        self.dbName = config["dbName"]
        self.LocalMySQLuser = config["LocalMySQLuser"]
        self.RemoteMySQLuser = config["RemoteMySQLuser"]
        self.RemoteMySQLhost = config["RemoteMySQLhost"]
        self.UpdateInterval = int(config["UpdateInterval"])
        self.TempFile = config["TempFile"]

        # Set the variables needed by the Output module.
        self.log = config["log"]
        self.datestring = config["datestring"]
    #__init__
#Config

#
# Unit test.
#
if __name__ == "__main__" :
    C = Config() # Will crash if the config file is not found.
    del C
#if
