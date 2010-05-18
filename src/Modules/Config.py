#!/usr/bin/python

class Config() :
    """
        Read the configuration file and store the data.

        Public variables:
            # Used by the Retriever module:
            email      ; Email address used for Entrez.
            cache      ; Location of the cache directory.
            cachesize  ; Maximum size of the cache directory in bytes.
            maxDldSize ; Maximum size of a GenBank record in bytes.
            minDldSize ; Minimum size of a GenBank record in bytes.

            # Used by the Db module:
            internalDb      ; Name of the internal database.
            dbNames         ; Name of the mapping databases
            LocalMySQLuser  ; Username for the local databases.
            RemoteMySQLuser ; Username for the remote UCSC database.
            RemoteMySQLhost ; Hostname of the UCSC database server.
            UpdateInterval  ; Time window (in days) to search for updates.
            TempFile        ; Location for downloaded updates.

            # Used by the Output module:
            log        ; Name and location of the logfile.
            datestring ; Prefix for log messages.

            # Used by the Mutator module:
            flanksize     ; Length of the flanking sequences in the 
                            visualisation.
            maxvissize    ; Maximum length of the variation in the 
                            visualisation.
            flankclipsize ; Length of the inserted/deleted flanks.


        Special Methods:
            __init__ ; Read the configuration file.
    """

    class Retriever() :
        pass
    #Retriever

    class Db() :
        pass
    #Db

    class Output() :
        pass
    #Output

    class Mutator() :
        pass
    #Mutator

    class Scheduler() :
        pass
    #Scheduler

    def __init__(self) :
        """
            Initialise the class with variables read from the configuration 
            file. In principle, this is the only place in the code where a
            hard coded constant is used (the name and path to the configuration
            file).

            Public variables (altered):
                # Used by the Retriever module:
                email      ; Email address used for Entrez.
                cache      ; Location of the cache directory.
                cachesize  ; Maximum size of the cache directory in bytes.
                maxDldSize ; Maximum size of a GenBank record in bytes.
                minDldSize ; Minimum size of a GenBank record in bytes.

                # Used by the Db module:
                internalDb      ; Name of the internal database.
                dbNames         ; Name of the mapping databases
                LocalMySQLuser  ; Username for the local databases.
                RemoteMySQLuser ; Username for the remote UCSC database.
                RemoteMySQLhost ; Hostname of the UCSC database server.
                UpdateInterval  ; Time window (in days) to search for updates.
                TempFile        ; Location for downloaded updates.

                # Used by the Output module:
                log        ; Name and location of the logfile.
                datestring ; Prefix for log messages.

                # Used by the Mutator module:
                flanksize     ; Length of the flanking sequences in the 
                                visualisation.
                maxvissize    ; Maximum length of the variation in the 
                                visualisation.
                flankclipsize ; Length of the inserted/deleted flanks.
        """
        from configobj import ConfigObj # ConfigObj()

        config = ConfigObj("./mutalyzer.conf")

        # Set the variables needed by the Retriever module.
        self.Retriever.email = config["email"]
        self.Retriever.cache = config["cache"]
        self.Retriever.cachesize = int(config["cachesize"]) * 1048576
        self.Retriever.maxDldSize = int(config["maxDldSize"]) * 1048576
        self.Retriever.minDldSize = int(config["minDldSize"])

        # Set the variables needed by the Db module.
        self.Db.internalDb = config["internalDb"]
        self.Db.dbNames = config["dbNames"]
        self.Db.LocalMySQLuser = config["LocalMySQLuser"]
        self.Db.RemoteMySQLuser = config["RemoteMySQLuser"]
        self.Db.RemoteMySQLhost = config["RemoteMySQLhost"]
        self.Db.UpdateInterval = int(config["UpdateInterval"])
        self.Db.TempFile = config["TempFile"]

        # Set the variables needed by the Output module.
        self.Output.log = config["log"]
        self.Output.datestring = config["datestring"]
        self.Output.loglevel = int(config["loglevel"])
        self.Output.outputlevel = int(config["outputlevel"])

        # Set the variables needed by the Mutator module.
        self.Mutator.flanksize = int(config["flanksize"])
        self.Mutator.maxvissize = int(config["maxvissize"])
        self.Mutator.flankclipsize = int(config["flankclipsize"])

        # Set the variables needed by the Scheduler module.
        self.Scheduler.watchDogTimeOut = int(config["watchDogTimeOut"])
        self.Scheduler.processName = config["processName"]
    #__init__
#Config

#
# Unit test.
#
if __name__ == "__main__" :
    C = Config() # Will crash if the config file is not found.
    del C
#if
