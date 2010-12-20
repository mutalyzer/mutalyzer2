#!/usr/bin/python

"""
Module for reading the config file and splitting up the variables into
subclasses. Each of these subclasses are used to configure a specific
module.
"""

class Config() :
    """
    Read the configuration file and store the data in subclasses.

    Special Methods:
        - __init__ ; Read the configuration file and initialise the
                     subclasses.
    """
    # Public subclasses:
    #     - Retriever ; Container for the Retriever configuration variables.
    #     - Db        ; Container for the Db configuration variables.
    #     - Output    ; Container for the Output configuration variables.
    #     - Mutator   ; Container for the Mutator configuration variables.
    #     - Scheduler ; Container for the Scheduler configuration variables.
    #     - File      ; Container for the File configuration variables.
    #     - GBparser  ; Container for the File configuration variables.


    class Retriever() :
        """
        Container class for the Retriever configuration variables.

        Public variables:
            - email      ; Email address used for Entrez.
            - cache      ; Location of the cache directory.
            - cachesize  ; Maximum size of the cache directory in bytes.
            - maxDldSize ; Maximum size of a GenBank record in bytes.
            - minDldSize ; Minimum size of a GenBank record in bytes.
            - lrgURL     ; base URL of LRG files.
        """

        pass
    #Retriever

    class Db() :
        """
        Container class for the Db configuration variables.

        Public variables:
            - internalDb      ; Name of the internal database.
            - dbNames         ; Name of the mapping databases
            - LocalMySQLuser  ; Username for the local databases.
            - LocalMySQLhost  ; Hostname of the local databases.

            - RemoteMySQLuser ; Username for the remote UCSC database.
            - RemoteMySQLhost ; Hostname of the UCSC database server.
            - UpdateInterval  ; Time window (in days) to search for
                                updates.
            - TempFile        ; Location for downloaded updates.
            """
    #Db

    class Output() :
        """
        Container class for the Output configuration variables.

        Public variables:
            - log         ; Name and location of the logfile.
            - datestring  ; Prefix for log messages.
            - loglevel    ; Default level for logging.
            - outputlevel ; Default level for output.
        """

        pass
    #Output

    class Mutator() :
        """
        Container class for the Mutator configuration variables.

        Public variables:
            - flanksize     ; Length of the flanking sequences in the
                              visualisation.
            - maxvissize    ; Maximum length of the variation in the
                              visualisation.
            - flankclipsize ; Length of the inserted/deleted flanks.
        """

        pass
    #Mutator

    class Scheduler() :
        """
        Container class for the Scheduler configuration variables.

        Public variables:
            - processName ; Name of the scheduler in the process list.
            - mailFrom    ; Return e-mail address.
            - mailMessage ; Template e-mail.
            - mailSubject ; Subject of the  e-mail.
            - resultsDir  ; Location of the results.
        """

        pass
    #Scheduler

    class Batch() :
        """
        Container class for the Scheduler configuration variables.

        Public variables:
            - PIDfile     ; Location of the PID file.
        """

        pass
    #Batch

    class File() :
        """
        Container class for the File configuration variables.

        Public variables:
            - bufSize   ; Amount of bytes to be read for determining the file
                          type.
            - header    ; The obligatory header in batch request files.
            - tempDir   ; Directory for temporary files.
            - threshold ; The threshold under which the percentage of errors
                          is allowed in a batchfile.
        """

        pass
    #File

    class GBparser() :
        """
        Container class for the GBparser configuration variables.

        Public variables:
            - upstream   ; Number of upstream nucleotides when searching for a
                           transcript.
            - downstream ; Number of downstream nucleotides when searching for a
                           transcript.
        """

        pass
    #GBparser

    class GenRecord() :
        pass

    def __init__(self) :
        """
        Initialise the class with variables read from the configuration
        file. In principle, this is the only place in the code where a
        hard coded constant is used (the name and path to the configuration
        file).

        Public subclasses (altered):
            - Retriever ; Initialised with Retriever configuration variables.
            - Db        ; Initialised with Db configuration variables.
            - Output    ; Initialised with Output configuration variables.
            - Mutator   ; Initialised with Mutator configuration variables.
            - Scheduler ; Initialised with Scheduler configuration variables.
        
        @requires: ConfigObj
        """
        from configobj import ConfigObj # ConfigObj()

        config = ConfigObj("./mutalyzer.conf")

        # Set the variables needed by the Retriever module.
        self.Retriever.email = config["email"]
        self.Retriever.cache = config["cache"]
        self.Retriever.cachesize = int(config["cachesize"]) * 1048576
        self.Retriever.maxDldSize = int(config["maxDldSize"]) * 1048576
        self.Retriever.minDldSize = int(config["minDldSize"])
        self.Retriever.lrgURL = config["lrgurl"]

        # Set the variables needed by the Db module.
        self.Db.internalDb = config["internalDb"]
        self.Db.dbNames = config["dbNames"]
        self.Db.LocalMySQLuser = config["LocalMySQLuser"]
        self.Db.LocalMySQLhost = config["LocalMySQLhost"]
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
        self.Scheduler.processName = config["processName"]
        self.Scheduler.mailFrom = config["mailFrom"]
        self.Scheduler.mailMessage = config["mailMessage"]
        self.Scheduler.mailSubject = config["mailSubject"]
        self.Scheduler.resultsDir = config["resultsDir"]
        self.Scheduler.nameCheckOutHeader = config["nameCheckOutHeader"]
        self.Scheduler.syntaxCheckOutHeader= config["syntaxCheckOutHeader"]
        self.Scheduler.positionConverterOutHeader= config["positionConverterOutHeader"]

        # Set thte variables neede for the Batch module.
        self.Batch.PIDfile = config["PIDfile"]

        # Set the variables needed by the File module.
        self.File.bufSize = int(config["bufSize"])
        self.File.header = config["header"]
        self.File.tempDir = config["tempDir"]
        self.File.threshold = float(config["threshold"])

        # Set the variables needed by the GBparser module.
        self.GBparser.email = config["email"]

        ## Set the variables needed by the File module.
        #self.File.upstream = int(config["upstream"])
        #self.File.downstream = int(config["downstream"])
        self.GenRecord.spliceAlarm = int(config["spliceAlarm"])
        self.GenRecord.spliceWarn = int(config["spliceWarn"])
    #__init__
#Config

#
# Unit test.
#
if __name__ == "__main__" :
    C = Config() # Will crash if the config file is not found.
    del C
#if
