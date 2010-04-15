#!/usr/bin/python

from Config import Config

class Output(Config) :
    """
        Provide an output interface for errors, warnings and logging purposes.

        Private variables:
            __instance   ; The name of the module that made this object.
            __loghandle  ; The handle of the log file.
            __datestring ; Format of the prefix for log messages.
            __errors     ; The number of errors that have been processed.
            __warnings   ; The number of warnings that have been processed.

        Special methods:
            __init__(config, instance) ; Initialise the class with variables
                                         from the config file and the calling
                                         module.
            __del__()                  ; Close the logfile.

        Private methods:                                         
            __NiceName(filename) ; Strip the path and the extention from a
                                   filename.

        Public methods:
            ErrorMsg(filename, message)   ; Print an error message to standard
                                            output and log it.
            WarningMsg(filename, message) ; Print an error message to standard
                                            output.
            LogMsg(filename, message)     ; Log a message.
            Summary()                     ; Print a summary of the number of
                                            errors and warnings.
    """

    def __init__(self, instance) :
        """ 
            Initialise the class private variables with variables from the 
            config file and the calling module.
            
            Arguments:
                config   ; The configuration object.
                instance ; The filename of the module that created this object.

            Private variables (altered):
                __instance   ; Initialised with the name of the module that
                               created this object.
                __loghandle  ; Initialised as the handle of the log file 
                               defined in the configuration file.
                __datestring ; Format of the prefix for log messages.
                __errors     ; Initialised to 0.
                __warnings   ; Initialised to 0.

            Inherited variables from Config:
                log ; Location of the log file.
        """

        Config.__init__(self)
        self.__instance = self.__NiceName(instance)
        self.__loghandle = open(self.log, "a")
        self.__errors = 0
        self.__warnings = 0
    #__init__

    def __del__(self) :
        """
            Close the log file.
            
            Private variables(altered):
                __loghandle ; The handle of the log file defined in the 
                             configuration file.
        """

        self.__loghandle.close()
    #__del__

    def __NiceName(self, filename) :
        """
            Strip the path and the extention from a filename.

            Arguments:
                filename ; A complete path plus extention.

            Returns:
                string ; The bare filename without a path and extention.
        """

        return filename.split('/')[-1].split('.')[0]
    #__NiceName

    def ErrorMsg(self, filename, message) :
        """
            Print an error message to standard output and log it.

            Arguments:
                filename ; The file where the error originated.
                message  ; The error message.

            Private variables (altered):
                __errors ; Increased by one.
        """

        print "Error (%s): %s" % (self.__NiceName(filename), message)
        self.LogMsg(filename, "Error: " + message)
        self.__errors += 1
    #ErrorMsg

    def WarningMsg(self, filename, message) :
        """
            Print an error message to standard output.

            Arguments:
                filename ; The file where the warning originated.
                message  ; The warning message.

            Private variables (altered):
                __warnings ; Increased by one.
        """

        print "Warning (%s): %s" % (self.__NiceName(filename), message)
        self.__warnings += 1
    #WarningMsg

    def LogMsg(self, filename, message) :
        """
            Log a message to the log file defined in the configuration file.

            Arguments:
                filename ; The file where the logging request originated.
                message  ; The message to be logged.

            Private variables:
                __loghandle  ; The handle of the log file defined in the 
                               configuration file.
                __instance   ; The name of the module that created this output
                               object.

            Inherited variables from Config:
                datestring ; Format of the prefix for log messages.
        """

        from time import strftime

        self.__loghandle.write(strftime(self.datestring + ' ') + \
            "%s (%s): %s\n" % (self.__instance, self.__NiceName(filename), 
            message))
        self.__loghandle.flush()
    #LogMsg

    def Summary(self) :
        """
            Print a summary of the number of errors and warnings.

            Private variables:
                __errors   ; The number of errors.
                __warnings ; The number of warnings.
        """

        e_s = 's'
        w_s = 's'
        if self.__errors == 1 :
            e_s = ''
        if self.__warnings == 1 :
            w_s = ''
            
        print "%i Error%s, %i Warning%s." % (self.__errors, e_s, 
                                             self.__warnings, w_s)
    #Summary
#Output

#
# Unit test.
#
if __name__ == "__main__" :
    O = Output(__file__)

    O.LogMsg(__file__, "hallo")
    O.ErrorMsg(__file__, "Ja, er ging wat mis.")
    del O
#if
