"""
Module for storing output and messages.
Output is stored as a named list that can be expanded.
Messages can be retrieved at a later time to provide flexibility. Message
levels are defined to increase or decrease the amount of logging and ouput.
The position of the log file, as well as the levels are defined in the
configuration file.

Message levels:
    - E{-}1 : Log     ; Specifically log a message.
    - 0 : Debug   ; Debug information.
    - 1 : Info    ; Info.
    - 2 : Warning ; Regular warnings.
    - 3 : Error   ; Serious errors that can be compensated for.
    - 4 : Fatal   ; Errors that are not recoverable.
    - 5 : Off     ; Can be used as a log/output level to turn off output.
"""
# Public classes:
#     - Message ; Container class for message variables.
#     - Output  ; Output interface for errors, warnings and logging.


import time

from mutalyzer import util
from mutalyzer.Serializers import SoapMessage


class Output() :
    """
    Provide an output interface for errors, warnings and logging purposes.

    Private variables:
        - __config     ; Configuration variables.
        - __outputdata ; The output dictionary.
        - __messages   ; The messages list.
        - __instance   ; The name of the module that made this object.
        - __loghandle  ; The handle of the log file.
        - __errors     ; The number of errors that have been processed.
        - __warnings   ; The number of warnings that have been processed.

    Special methods:
        - __init__(instance, config) ; Initialise the class with variables
                                       from the config file and the calling
                                       module.
        - __del__()                  ; Close the logfile and clean up.

    Public methods:
        - addMessage(filename, level, code, description) ; Add a message to
                                                           the message list.
        - getMessages()           ; Print all messages that exceed the
                                    configured output level.
        - addOutput(name, data)   ; Add output to the output dictionary.
        - getOutput(name)         ; Retrieve data from the output dictionary.
        - Summary()               ; Print a summary of the number of errors
                                    and warnings.
    """

    def __init__(self, instance, config) :
        """
        Initialise the class private variables with variables from the
        config file and the calling module.

        Private variables (altered):
            - __config     ; Configuration variables.
            - __outputdata ; The output dictionary.
            - __messages   ; The messages list.
            - __instance   ; Initialised with the name of the module that
                             created this object.
            - __loghandle  ; Initialised as the handle of the log file
                             defined in the configuration file.
            - __errors     ; Initialised to 0.
            - __warnings   ; Initialised to 0.

        @arg instance: The filename of the module that created this object
        @type instance: string
        @arg config: The configuration object
        @type config: object
        """

        self.__config = config
        self.__outputData = {}
        self.__messages = []
        self.__instance = util.nice_filename(instance)
        self.__loghandle = open(self.__config.log, "a+")
        self.__errors = 0
        self.__warnings = 0
    #__init__

    def __del__(self) :
        """
        Clean up the output dictionary, the messages list and close the log
        file.

        Private variables(altered):
            - __loghandle  ; The handle of the log file defined in the
                             configuration file.
            - __outputdata ; The output dictionary.
            - __messages   ; The messages list.
        """

        self.__loghandle.close()
        for i in self.__outputData :
            del i
        for i in self.__messages :
            del i
    #__del__

    def addMessage(self, filename, level, code, description) :
        """
        Add a message to the message list.
        If the level exceeds the configured loglevel or if the level is -1,
        then the message is also logged.
        If the severity equals 2, then the number of warnings is inreased,
        if it exceeds 2, then the number of errors is increased.

        Private variables:
            - __messages  ; The messages list.
            - __instance  ; Module that created the Output object.
            - __config    ; The variables loglevel and datestring are used.
            - __loghandle ; Handle to the log file.

        Private variables (altered):
            - __warnings ; Increased by one if the severity equals 2.
            - __errors   ; Increased by one if the severity exceeds 2.

        @arg filename:    Name of the calling module
        @arg level:       Severity of the message
        @arg code:        Error code of the message
        @arg description: Description of the message
        """
        nice_name = util.nice_filename(filename)
        message = Message(nice_name, level, code, description)

        # Append a new message object to the messages list.
        self.__messages.append(message)

        if level == 2:
            self.__warnings += 1
        if level > 2:
            self.__errors += 1

        # Log the message if the message is important enough, or if it is only
        # meant to be logged (level -1).
        if level >= self.__config.loglevel or level == -1 :
            self.__loghandle.write(time.strftime(
                self.__config.datestring + ' ') + "%s (%s) %s: %s: %s\n" % (
                    self.__instance, nice_name, code, message.named_level(),
                    description))
            self.__loghandle.flush()
        #if
    #addMessage

    def getMessages(self) :
        """
        Print all messages that exceed the configured output level.

        Private variables:
            - __messages  ; The messages list.
            - __config    ; The variable outputlevel is used.

        @return: A list of messages
        @rtype: list
        """
        return filter(lambda m: m.level >= self.__config.outputlevel,
                      self.__messages)
    #getMessages

#        ret = []
#        for i in self.__messages :
#            if i.level >= self.__config.outputlevel :
#                ret.append('%s: (%s): %s' % (i.named_level(), i.origin,
#                                             i.description))
#        return ret

    def getSoapMessages(self):
        """
        Returns a list of SoapMessages for over the wire

        Private variables:
            - __messages  ; The messages list.
            - __config    ; The variable outputlevel is used.

        @requires: Modules.Serializers.SoapMessage

        @return: list of SoapMessages
        @rtype: list
        """
        ret = []
        for i in self.__messages:
            if i.level >= self.__config.outputlevel:
                mess = SoapMessage()
                mess.errorcode = i.code
                mess.message = i.description
                ret.append(mess)
            #if
        #for
        return ret
    #getSoapMessages

    def getBatchMessages(self, level):
        """
        Returns a list of Messages with an errorlevel >= level
        and removes additional lines from a parseerror

        Private variables:
            - __messages   ; The messages list.

        @arg level: error level
        @type level: integer

        @return: list of Messages
        @rtype: list
        """

        ret = []
        lastorigin = ""
        for i in self.__messages:
            if i.level >= level:
                # Todo: We changed this from 'Parser' to 'grammar', does this
                # still work?
                if lastorigin == 'grammar': #Only one parse error
                    continue
                lastorigin = i.origin
                ret.append("(%s): %s" % (i.origin, i.description))
            #if
        #for
        return ret
    #getBatchMessages


    def addOutput(self, name, data) :
        """
        If the output dictionary already has a node with the specified
        name, the list that this name points to is expanded with the data.
        Otherwise create a node and assign a list containing the data.

        Private variables:
            - __outputData ; The output dictionary.

        @arg name: Name of a node in the output dictionary
        @type name: string
        @arg data: The data to be stored at this node
        @type data: object
        """

        if self.__outputData.has_key(name) :
            self.__outputData[name].append(data)
        else :
            self.__outputData[name] = [data]
    #addOutput

    def getOutput(self, name) :
        """
        Return a list of data from the output dictionary.

        Private variables:
            - __outputData ; The output dictionary.

        @arg name: Name of a node in the output dictionary
        @type name: string

        @return: output dictionary
        @rtype: dictionary
        """

        if self.__outputData.has_key(name) :
            return self.__outputData[name]
        return []
    #getOutput

    def getIndexedOutput(self, name, index) :
        """
        Return an element of a list, the list is called 'name' in de
        __outputData dictionary. If either the list or the element does not
        exist, return None.

        @arg name:  Name of the list.
        @arg index: Index of the element to be retuned.

        Private variables:
            - __outputData ; The output dictionary.

        @return: The requested element or None
        @rtype: any type
        """

        if self.__outputData.has_key(name) :
            if 0 <= index < len(self.__outputData[name]) :
                return self.__outputData[name][index]
        return None
    #getFirst

    def getMessagesWithErrorCode(self, errorcode):
        """
        Retrieve all messages that have a specific error code.

        Private variables:
            - __messages   ; The messages list.

        @arg errorcode: The error code to filter on
        @type errorcode: string

        @return: A filtered list
        @rtype: list
        """
        return filter(lambda m: m.code == errorcode, self.__messages)
    #getMessagesWithErrorCode

    def Summary(self) :
        """
        Print a summary of the number of errors and warnings.

        Private variables:
            - __errors   ; The number of errors.
            - __warnings ; The number of warnings.

        @return:
            triple:
                - Number of errors
                - Number of warnings
                - Summary
        @rtype: integer, integer, string
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

class Message() :
    """
    Container class for message variables.

    Special methods:
        - __init__(origin, level, code, description) ; Make a message object.

    Public variables:
        - origin      ; Name of the module creating this object.
        - level       ; Importance of the message.
        - code        ; The error code of the message.
        - description ; A description of the message.
    """

    def __init__(self, origin, level, code, description) :
        """
        Make a new message object.

        Public variables (altered):
            - origin      ; Name of the module creating this object.
            - level       ; Importance of the message.
            - code        ; The error code of the message.
            - description ; A description of the message.

        @arg origin: Name of the module creating this object
        @type origin: string
        @arg level: Importance of the message
        @type level: integer
        @arg code: The error code of the message
        @type code: string
        @arg description: A description of the message
        @type description: string
        """
        self.origin = origin
        self.level = level
        self.code = code
        self.description = description
    #__init__

    def __repr__(self):
        return 'Message("%s", %i, "%s", "%s")' % \
               (self.origin, self.level, self.code, self.description)

    def __str__(self):
        return '%s (%s): %s' % \
               (self.named_level(), self.origin, self.description)

    def named_level(self):
        """
        Get message log level as readable string.

        @return:     A readable description of the log level.
        @rtype:      string
        """
        if self.level == 0:
            return "Debug"
        if self.level == 1:
            return "Info"
        if self.level == 2:
            return "Warning"
        if self.level == 3:
            return "Error"
        if self.level == 4:
            return "Fatal"
        return ''
    #named_level
#Message
