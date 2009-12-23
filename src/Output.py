class Output() :
    def __init__(self, config) :
        self.__log = config.log
        self.__errors = 0
        self.__warnings = 0
    #__init__

    def __NiceName(self, filename) :
        return filename.split('/')[-1].split('.')[0]
    #__NiceName

    def ErrorMsg(self, filename, message) :
        print "Error (%s): %s" % (self.__NiceName(filename), message)
        self.LogMsg(filename, "Error: " + message)
        self.__errors += 1
    #ErrorMsg

    def WarningMsg(self, filename, message) :
        print "Warning (%s): %s" % (self.__NiceName(filename), message)
        self.__warnings += 1
    #WarningMsg

    def LogMsg(self, filename, message) :
        from time import strftime

        handle = open(self.__log, "a")
        handle.write(strftime("%Y-%m-%d %H:%M:%S ") + \
            self.__NiceName(filename) + ": " + message + "\n")
        handle.close()
    #LogMsg

    def Summary(self) :
        print "%i Errors, %i Warnings." % (self.__errors, self.__warnings)
#Output

if __name__ == "__main__" :
    import Config

    C = Config.Config()
    O = Output(C)

    O.LogMsg(__file__, "hallo")
    O.ErrorMsg(__file__, "Ja, er ging wat mis.")
#if
