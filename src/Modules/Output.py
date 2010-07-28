#!/usr/bin/python

"""
    Module for storing output and messages.
    Output is stored as a named list that can be expanded.
    Messages can be retrieved at a later time to provide flexibility. Message
    levels are defined to increase or decrease the amount of logging and ouput.
    The position of the log file, as well as the levels are defined in the
    configuration file.

    Message levels:
        -1 : Log     ; Specifically log a message.
         0 : Debug   ; Debug information.
         1 : Info    ; Info.
         2 : Warning ; Regular warnings.
         3 : Error   ; Serious errors that can be compensated for.
         4 : Fatal   ; Errors that are not recoverable.
         5 : Off     ; Can be used as a log/output level to turn off output.

    Public classes:
        Message ; Container class for message variables.
        Output  ; Output interface for errors, warnings and logging.
"""

import time # strftime()

class Message() :
    """
        Container class for message variables.

        Special methods:
            __init__(origin, level, code, description) ; Make a message object.

        Public variables:
            origin      ; Name of the module creating this object.
            level       ; Importance of the message.
            code        ; The error code of the message.
            description ; A description of the message.
    """

    def __init__(self, origin, level, code, description) :
        """
            Make a new message object.

            Arguments:
                origin      ; Name of the module creating this object.
                level       ; Importance of the message.
                code        ; The error code of the message.
                description ; A description of the message.

            Public variables (altered):
                origin      ; Name of the module creating this object.
                level       ; Importance of the message.
                code        ; The error code of the message.
                description ; A description of the message.
        """

        self.origin = origin
        self.level = level
        self.code = code
        self.description = description
    #__init__
#Message

class Output() :
    """
        Provide an output interface for errors, warnings and logging purposes.

        Private variables:
            __config     ; Configuration variables.
            __outputdata ; The output dictionary.
            __messages   ; The messages list.
            __instance   ; The name of the module that made this object.
            __loghandle  ; The handle of the log file.
            __errors     ; The number of errors that have been processed.
            __warnings   ; The number of warnings that have been processed.

        Special methods:
            __init__(instance, config) ; Initialise the class with variables
                                         from the config file and the calling
                                         module.
            __del__()                  ; Close the logfile and clean up.

        Private methods:
            __niceName(filename) ; Strip the path and the extention from a
                                   filename.
            __levelToName(level) ; Convert a log level to a readable string.

        Public methods:
            addMessage(filename,    ; Add a message to the message list.
                       level,
                       code,
                       description)
            getMessages()           ; Print all messages that exceed the
                                      configured output level.
            addOutput(name, data)   ; Add output to the output dictionary.
            getOutput(name)         ; Retrieve data from the output dictionary.
            Summary()               ; Print a summary of the number of errors
                                      and warnings.
    """

    def __init__(self, instance, config) :
        """
            Initialise the class private variables with variables from the
            config file and the calling module.

            Arguments:
                instance ; The filename of the module that created this object.
                config   ; The configuration object.

            Private variables (altered):
                __config     ; Configuration variables.
                __outputdata ; The output dictionary.
                __messages   ; The messages list.
                __instance   ; Initialised with the name of the module that
                               created this object.
                __loghandle  ; Initialised as the handle of the log file
                               defined in the configuration file.
                __errors     ; Initialised to 0.
                __warnings   ; Initialised to 0.
        """

        self.__config = config
        self.__outputData = {}
        self.__messages = []
        self.__instance = self.__niceName(instance)
        self.__loghandle = open(self.__config.log, "a")
        self.__errors = 0
        self.__warnings = 0
    #__init__

    def __del__(self) :
        """
            Clean up the output dictionary, the messages list and close the log
            file.

            Private variables(altered):
                __loghandle  ; The handle of the log file defined in the
                               configuration file.
                __outputdata ; The output dictionary.
                __messages   ; The messages list.
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
            Convert a log level to a readable string.

            Arguments:
                level ; A log level (an integer between -1 and 5).

            Returns:
                string ; A readable description of the log level.
        """

        if level == 0 :
            return "Debug: "
        if level == 1 :
            return "Info: "
        if level == 2 :
            return "Warning: "
        if level == 3 :
            return "Error: "
        if level == 4 :
            return "Fatal: "
        return ""
    #__levelToName

    def addMessage(self, filename, level, code, description) :
        """
            Add a message to the message list.
            If the level exceeds the configured loglevel or if the level is -1,
            then the message is also logged.
            If the severity equals 2, then the number of warnings is inreased,
            if it exceeds 2, then the number of errors is increased.

            Arguments:
                filename    ; Name of the calling module.
                level       ; Severity of the message.
                code        ; Error code of the message.
                description ; Description of the message.

            Private variables:
                __messages  ; The messages list.
                __instance  ; Module that created the Output object.
                __config    ; The variables loglevel and datestring are used.
                __loghandle ; Handle to the log file.

            Private variables (altered):
                __warnings ; Increased by one if the severity equals 2.
                __errors   ; Increased by one if the severity exceeds 2.
        """

        niceName = self.__niceName(filename)

        # Append a new message object to the messages list.
        self.__messages.append(Message(niceName, level, code, description))

        if level == 2 :
            self.__warnings += 1
        if level > 2 :
            self.__errors += 1

        # Log the message if the message is important enough, or if it is only
        # meant to be logged (level -1).
        if level > self.__config.loglevel or level == -1 :
            self.__loghandle.write(time.strftime(
                self.__config.datestring + ' ') + "%s (%s) %s: %s%s\n" % (
                    self.__instance, niceName, code, self.__levelToName(level),
                    description))
            self.__loghandle.flush()
        #if
    #addMessage

    def getMessages(self) :
        """
            Print all messages that exceed the configured output level.

            Private variables:
                __messages  ; The messages list.
                __config    ; The variable outputlevel is used.
        """

        ret = []
        for i in self.__messages :
            if i.level > self.__config.outputlevel :
                #print "%s(%s): %s" % (self.__levelToName(i.level), i.origin,
                #                      i.description)
                ret.append("%s(%s): %s" % (self.__levelToName(i.level),
                                           i.origin, i.description))
        return ret
    #getMessages

    def getSoapMessages(self):
        """ Returns a list of SoapMessages for over the wire """
        #TODO: MOVE to top if works
        from Modules.Mapper import SoapMessage

        ret = []
        for i in self.__messages:
            if i.level > self.__config.outputlevel:
                mess = SoapMessage()
                mess.errorcode = i.code
                mess.message = i.description
                ret.append(mess)
        return ret

    def getBatchMessages(self, level):
        """ Returns a list of Messages with an errorlevel >= level

            and removes additional lines from a parseerror
        """
        ret = []
        lastorigin = ""
        for i in self.__messages:
            if i.level >= level:
                if lastorigin == "Parser": continue #Only one parse error
                lastorigin = i.origin
                ret.append("(%s): %s" % (i.origin, i.description))
        return ret


    def addOutput(self, name, data) :
        """
            If the output dictionary already has a node with the specified
            name, the list that this name points to is expanded with the data.
            Otherwise create a node and assign a list containing the data.

            Arguments:
                name ; Name of a node in the output dictionary.
                data ; The data to be stored at this node.

            Private variables:
                __outputData ; The output dictionary.
        """

        if self.__outputData.has_key(name) :
            self.__outputData[name].append(data)
        else :
            self.__outputData[name] = [data]
    #addOutput

    def getOutput(self, name) :
        """
            Return a list of data from the output dictionary.

            Arguments:
                name ; Name of a node in the output dictionary.

            Private variables:
                __outputData ; The output dictionary.
        """

        if self.__outputData.has_key(name) :
            return self.__outputData[name]
        return []
    #getOutput

    def getIndexedOutput(self, name, index) :
        """
        """

        if self.__outputData.has_key(name) :
            if 0 <= index < len(self.__outputData[name]) :
                return self.__outputData[name][index]
        return []
    #getFirst

    def getMessagesWithErrorCode(self, errorcode):
        ret = []
        for i in self.__messages:
            if i.code == errorcode:
                ret.append(i)
        return ret


    def Summary(self) :
        """
            Print a summary of the number of errors and warnings.

            Private variables:
                __errors   ; The number of errors.
                __warnings ; The number of warnings.

            Returns:
                triple:
                    integer ; Number of errors.
                    integer ; Number of warnings.
                    string  ; Summary.
        """

        e_s = 's'
        w_s = 's'
        if self.__errors == 1 :
            e_s = ''
        if self.__warnings == 1 :
            w_s = ''

        return self.__errors, self.__warnings, "%i Error%s, %i Warning%s." % (
            self.__errors, e_s, self.__warnings, w_s)
    #Summary
#Output

#
# Unit test.
#
if __name__ == "__main__" :
    pass
#if
