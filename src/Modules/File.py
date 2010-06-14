#!/usr/bin/python

"""
    Module for parsing CSV files and spreadsheets.

    Public classes:
        File ; Parse CSV files and spreadsheets.
"""

import magic           # open(), MAGIC_MIME, MAGIC_NONE
import csv             # Sniffer(), reader(), Error
import xlrd            # open_workbook()
import zipfile         # ZipFile()
import xml.dom.minidom # parseString()
import os              # remove()
import types           # UnicodeType

from Modules import Misc

class File() :
    """
        Parse CSV files and spreadsheets.

        Private variables:
                __config ; Configuration variables.
                __output ; The Output object.

        Special methods:
            __init__(config, output) ; Initialse the class.

        Private methods:
            __tempFileWrapper(func,   ; Call func() with a filename.
                              handle)
            __getMimeType(handle)     ; Get the mime type of a stream.
            __parseCsvFile(handle)    ; Parse a CSV file.
            __parseXlsFile(handle)    ; Parse an Excel file.
            __parseOdsFile(handle)    ; Parse an OpenDocument Spreadsheet file.
            __checkBatchFormat(job)   ; Check a batch job and sanitize it.

        Public methods:
            parseFileRaw(handle)   ; Parse a stream with the appropriate parser.
            parseBatchFile(handle) ; Parse a stream with the appropriate parser
                                     and sanitize the output.
    """

    def __init__(self, config, output) :
        """
            Initialise the class.

            Private variables (altered):
                __config ; Initialised with configuration variables.
                __output ; Set to the Output object.
        """

        self.__config = config
        self.__output = output
    #__init__

    def __tempFileWrapper(self, func, handle) :
        """
            Make a temporary file, put the content of a stream in it and pass
            the filename to a general function. Return whatever this function
            returns.

            Arguments:
                func   ; A general function that needs a file name as argument.
                handle ; A stream.

            Returns:
                unknown ; The output of func().
        """

        # Generate an unique filename in the tempDir directory.
        MiscInstance = Misc.Misc()
        fileName = self.__config.tempDir + '/' + str(MiscInstance.ID())
        del MiscInstance

        # Dump the content of the stream pointed to by handle into the file.
        handle.seek(0)
        writeHandle = open(fileName, "w")
        writeHandle.write(handle.read())
        writeHandle.close()

        # Open the file with func().
        ret = func(fileName)
        os.remove(fileName)

        return ret
    #__tempFileWrapper

    def __getMimeType(self, handle) :
        """
            Get the mime type of a stream by inspecting a fixed number of bytes.
            The stream is not rewinded after use.

            Arguments:
                handle ; A handle to a stream.

            Private variables:
                __config ; The bufSize configuration variables.

            Returns:
                string ; The mime type of a file.
        """

        handle.seek(0)
        buf = handle.read(self.__config.bufSize)

        MagicInstance = magic.open(magic.MAGIC_MIME)
        MagicInstance.load()
        mimeType = MagicInstance.buffer(buf).split(';')[0]
        MagicInstance.close()
        MagicInstance = magic.open(magic.MAGIC_NONE)
        MagicInstance.load()
        description = MagicInstance.buffer(buf)
        del MagicInstance
        
        return mimeType, description
    #__getMimeType

    def __parseCsvFile(self, handle) :
        """
            Parse a CSV file.
            The stream is not rewinded after use.

            Arguments:
                handle ; A handle to a stream.

            Private variables:
                __config ; The bufSize configuration variables.

            Returns:
                list ; A list of lists.
        """

        handle.seek(0)
        buf = handle.read(self.__config.bufSize)

        try :
            dialect = csv.Sniffer().sniff(buf)
        except csv.Error, e :
            self.__output.addMessage(__file__, 4, "EBPARSE", e)
            return None
        #except

        handle.seek(0)
        reader = csv.reader(handle, dialect)

        ret = []
        for i in reader :
            ret.append(i)

        return ret
    #__parseCsvFile

    def __parseXlsFile(self, handle) :
        """
            Parse an Excel file.
            The stream is not rewinded after use.

            Arguments:
                handle ; A handle to a stream.

            Returns:
                list ; A list of lists.
        """

        workBook = self.__tempFileWrapper(xlrd.open_workbook, handle)
        sheet = workBook.sheet_by_index(0)

        ret = []
        for i in range(sheet.nrows) :
            row = []
            for j in sheet.row_values(i) :
                if type(j) == types.UnicodeType : # Convert the data to strings.
                    row.append(j.encode("utf8"))
                else :
                    row.append(str(j))
            #for
            ret.append(row)
        #for

        del sheet, workBook

        return ret
    #__parseXlsFile

    def __parseOdsFile(self, handle) :
        """
            Parse an OpenDocument Spreadsheet file.
            The stream is not rewinded after use.

            Arguments:
                handle ; A handle to a stream.

            Returns:
                list ; A list of lists.

        """

        zipFile = self.__tempFileWrapper(zipfile.ZipFile, handle)
        doc = xml.dom.minidom.parseString(zipFile.read("content.xml"))
        zipFile.close()

        ret = []
        for i in doc.getElementsByTagName("table:table-row") :
            row = []
            for j in i.getElementsByTagName("table:table-cell") :
                c = j.getElementsByTagName("text:p")
                if c :
                    row.append(c[0].lastChild.data.encode("utf8"))
                #if
            #for
            ret.append(row)
        #for

        return ret
    #__parseOdsFile

    def __checkBatchFormat(self, job) :
        """
            Check if a job is of the correct format.
            - Each row should consist of three elements.
            - The first and the last element should be non-zero.
            - The first line should be the header defined in the config file.
            - Silently ignore all empty lines.

            Arguments:
                job ; list of lists.

            Private variables:
                __config ; The header configuration variable.

            Returns:
                list ; A sanitised list of lists (without a header or empty
                       lines).
        """
        
        if job[0] != self.__config.header :
            self.__output.addMessage(__file__, 4, "EBPARSE", 
                                     "Header not valid.")
            return None
        #if

        for i in range(0, len(job)) :
            if job[i] :                                      # Non empty line.
                if len(job[i]) == 3 :
                    if job[i][0] or job[i][1] or job[i][2] : # Non empty line.
                        if not job[i][0] :
                            self.__output.addMessage(__file__, 4, "EBPARSE",
                                "The first column may not be empty in line " \
                                "%i." % i)
                            return None
                        #if
                        if not job[i][2] :
                            self.__output.addMessage(__file__, 4, "EBPARSE",
                                "The last column may not be empty in line " \
                                "%i." % i)
                            return None
                        #if
                    #if
                #if
                else :
                    self.__output.addMessage(__file__, 4, "EBPARSE",
                        "Wrong amount of columns in line %i.\n" % i)
                    return None
                #else
            #if
        #for

        # All tests are passed, now we do some trimming.
        ret = []
        for i in range(1, len(job)) :
            if job[i] and job[i] != ['', '', ''] :
                ret.append(job[i])

        return ret
    #__checkBatchFormat

    def parseFileRaw(self, handle) :
        """
            Check which format a stream has and parse it with the appropriate
            parser if the stream is recognised.

            Arguments:
                handle ; A handle to a stream.

            Returns:
                list ; A list of lists, None if an error occured.
        """

        mimeType = self.__getMimeType(handle)
        if mimeType[0] == "text/plain" :
            return self.__parseCsvFile(handle)
        if mimeType[0] == "application/vnd.ms-office" :
            return self.__parseXlsFile(handle)
        if mimeType == ("application/octet-stream", 
                        "OpenDocument Spreadsheet") :
            return self.__parseOdsFile(handle)

        return None
    #parseFile

    def parseBatchFile(self, handle) :
        """
            Check which format a stream has and parse it with the appropriate
            parser if the stream is recognised.

            Arguments:
                handle ; A handle to a stream.

            Returns:
                list ; A sanitised list of lists (without a header or empty
                       lines), or None if an error occured.
        """

        job = self.parseFileRaw(handle)
        if job :
            return self.__checkBatchFormat(job)
        return None
    #parseBatchFile
#File
