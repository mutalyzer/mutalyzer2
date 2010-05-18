#!/usr/bin/python

#from Config import Config
from time import strftime

class Node() :
    """
    """

    def __init__(self, level) :
        self.message = []
        self.level = level
    #__init__
#Node

class Message() :
    """
    """

    def __init__(self, origin, level, code, description) :
        self.origin = origin
        self.level = level
        self.code = code
        self.description = description
    #__init__
#Message

#class Empty() :
#    def __len__(self) :
#        return 0

class Output() :
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
            __niceName(filename) ; Strip the path and the extention from a
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

    def __init__(self, instance, config) :
        """ 
            Initialise the class private variables with variables from the 
            config file and the calling module.
            
            Arguments:
                config   ; The configuration object.
                instance ; The filename of the module that created this object.

            Public variables(altered):
                outputdata ; The output list.

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

        #Config.__init__(self)
        self.__config = config
        self.__outputData = {}
        self.__messages = []
        self.__instance = self.__niceName(instance)
        self.__loghandle = open(self.__config.log, "a")
        self.__errors = 0
        self.__warnings = 0


        #self.createOutputNode("debug", 0)
        #self.createOutputNode("info", 1)
        #self.createOutputNode("warnings", 2)
        #self.createOutputNode("errors", 3)
        #self.createOutputNode("fatalerrors", 4)
        #self.createOutputNode("log", 5)
    #__init__

    def __del__(self) :
        """
            Clean up the output list and close the log file.
            
            Public variables(altered):
                outputdata ; The output list.

            Private variables(altered):
                __loghandle ; The handle of the log file defined in the 
                             configuration file.
        """

        self.__loghandle.close()
        for i in self.__outputData :
            del i
        for i in self.__messages :
            del i
    #__del__

    def __niceName(self, filename) :
        """
            Strip the path and the extention from a filename.

            Arguments:
                filename ; A complete path plus extention.

            Returns:
                string ; The bare filename without a path and extention.
        """

        return filename.split('/')[-1].split('.')[0]
    #__niceName

    def __levelToName(self, level) :
        """
        """

        if level == 0 :
            return "Debug"
        if level == 1 :
            return "Info"
        if level == 2 :
            return "Warning"
        if level == 3 :
            return "Error"
        if level == 4 :
            return "Fatal"
        return ""
    #__levelToName

    #def addToOutputNode(self, filename, name, code, message) :
    #    """
    #    """

    #    niceName = self.__niceName(filename)

    #    self.__outputData[name].message.append(Message(niceName, code, message))

    #    level = self.__outputData[name].level
    #    if level >= self.__config.loglevel :
    #        prefix = ""
    #        if level == 2 :
    #            prefix = "Warning: "
    #        if level == 3 :
    #            prefix = "Error: "
    #        if level == 4 :
    #            prefix = "Fatal: "
    #        self.__loghandle.write(strftime(self.__config.datestring + ' ') + \
    #                               "%s (%s) %s: %s%s\n" % (self.__instance, 
    #                               niceName, code, prefix, message))
    #        self.__loghandle.flush()
    #    #if
    ##addToOutputNode

    def addMessage(self, filename, level, code, description) :
        """
        """

        niceName = self.__niceName(filename)

        self.__messages.append(Message(niceName, level, code, description))

        if level == 2 :
            self.__warnings += 1
        if level > 2 :
            self.__errors += 1
        if level > self.__config.loglevel or level == -1 :
            self.__loghandle.write(strftime(self.__config.datestring + ' ') + \
                                   "%s (%s) %s: %s%s\n" % (self.__instance, 
                                   niceName, code, self.__levelToName(level), 
                                   description))
            self.__loghandle.flush()
        #if
    #addMessage

    def getMessages(self) :
        """
        """

        for i in self.__messages :
            if i.level > self.__config.outputlevel :
                print "%s (%s): %s" % (self.__levelToName(i.level), i.origin,
                                       i.description)
    #getMessages

    def addOutput(self, name, data) :
        """
        """

        if self.__outputData.has_key(name) :
            self.__outputData[name].append(data)
        else :
            self.__outputData[name] = [data]
    #addOutput

    def getOutput(self, name) :
        """
        """

        if self.__outputData.has_key(name) :
            return self.__outputData[name]
        return None
    #getOutput


    #def createOutputNode(self, name, level) :
    #    """
    #    """

    #    self.__outputData[name] = Node(level)
    ##createOutputNode

    #def getData(self, name) :
    #    """
    #    """

    #    if self.__outputData.has_key(name) and \
    #       self.__outputData[name].level >= self.__config.outputlevel :
    #        return self.__outputData[name].message
    #    return []
    ##getdata

    #'''
    #def getMsg(self, name) :
    #    """
    #    """

    #    if self.__outputData[name].serverity >= self.__config.outputlevel :
    #        return 
    ##getMsg    

    #def ErrorMsg(self, filename, message) :
    #    """
    #        Print an error message to standard output and log it.

    #        Arguments:
    #            filename ; The file where the error originated.
    #            message  ; The error message.

    #        Private variables (altered):
    #            __errors ; Increased by one.
    #    """

    #    print "Error (%s): %s" % (self.__niceName(filename), message)
    #    self.LogMsg(filename, "Error: " + message)
    #    #self.addData(filename, "test", "error", "5", message)
    #    self.__errors += 1
    ##ErrorMsg

    #def WarningMsg(self, filename, message) :
    #    """
    #        Print an error message to standard output.

    #        Arguments:
    #            filename ; The file where the warning originated.
    #            message  ; The warning message.

    #        Private variables (altered):
    #            __warnings ; Increased by one.
    #    """

    #    print "Warning (%s): %s" % (self.__niceName(filename), message)
    #    #self.addData(filename, "test", "warning", "5", message)
    #    self.__warnings += 1
    ##WarningMsg

    #def LogMsg(self, filename, message) :
    #    """
    #        Log a message to the log file defined in the configuration file.

    #        Arguments:
    #            filename ; The file where the logging request originated.
    #            message  ; The message to be logged.

    #        Private variables:
    #            __loghandle  ; The handle of the log file defined in the 
    #                           configuration file.
    #            __instance   ; The name of the module that created this output
    #                           object.

    #        Inherited variables from Config:
    #            datestring ; Format of the prefix for log messages.
    #    """


    #    self.__loghandle.write(strftime(self.__config.datestring + ' ') + \
    #        "%s (%s): %s\n" % (self.__instance, self.__niceName(filename), 
    #        message))
    #    self.__loghandle.flush()
    ##LogMsg
    #'''

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
    import Config

    C = Config.Config()

    O = Output(__file__, C.Output)

    O.WarningMsg(__file__, "Ja, er ging wat mis.")
    del O
#if
