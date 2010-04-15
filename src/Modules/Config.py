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
            dbNames         ; The names of the databases.
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
                dbNames    ; The name of the databases.
        """
        from configobj import ConfigObj # ConfigObj()

        config = ConfigObj("./mutalyzer.conf")

        # Set the variables needed by the Retriever module.
        self.email = config["email"]
        self.cache = config["cache"]
        self.cachesize = int(config["cachesize"])

        # Set the variables needed by the Db module.
        self.dbNames = config["dbNames"]
        self.LocalMySQLuser = config["LocalMySQLuser"]
        self.RemoteMySQLuser = config["RemoteMySQLuser"]
        self.RemoteMySQLhost = config["RemoteMySQLhost"]
        self.UpdateInterval = int(config["UpdateInterval"])
        self.TempFile = config["TempFile"]

        # Set the variables needed by the Output module.
        self.log = config["log"]
        self.datestring = config["datestring"]

        # Set the variables needed by the Mutator module.
        self.flanksize = int(config["flanksize"])
        self.maxvissize = int(config["maxvissize"])
        self.flankclipsize = int(config["flankclipsize"])
    #__init__
#Config

#
# Unit test.
#
if __name__ == "__main__" :
    C = Config() # Will crash if the config file is not found.
    del C
#if
