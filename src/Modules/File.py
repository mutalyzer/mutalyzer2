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
        # Apperantly apache will remove this file even when opened by the 
        # function *func
        os.remove(fileName)

        return ret
    #__tempFileWrapper

    def getMimeType(self, handle) :
        """
            Get the mime type of a stream by inspecting a fixed number of bytes.
            The stream is rewinded after use.

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
        handle.seek(0)
        
        return mimeType, description
    #getMimeType

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

        #Watch out for : delimiter FIXME and for the . delimiter
        if dialect.delimiter == ":":
            dialect.delimiter = "\t"

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

        #zipFile = self.__tempFileWrapper(zipfile.ZipFile, handle)
        zipFile = zipfile.ZipFile(handle)
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
        #remove empty lines (store original line numbers line 1 = job[0])
        jobl = [(l+1, row) for l, row in enumerate(job) if row and any(row)]

        #TODO:  Add more new style old style logic
        #       Should empty lines be stored
        #       Can we concatonate 3 column entries to 1 column

        if jobl[0][1] == self.__config.header : #Old style NameCheckBatch job
            #collect all lines where the amount of arguments != 3
            notthree = filter(lambda i: len(i[1])!=3, jobl)

            [jobl.remove(r) for r in notthree]       # update job

            #collect all lines where the first or third argument is empty
            emptyfield = filter(lambda i: not(i[1][0] and i[1][2]), jobl)

            [jobl.remove(r) for r in emptyfield]     # update job

            #Create output Message for incompatible fields
            if any(notthree):
                lines = ", ".join([str(i[0]) for i in notthree])
                self.__output.addMessage(__file__, 4, "EBPARSE",
                        "Wrong amount of columns in line(s): %s.\n" % lines)

            if any(emptyfield):
                lines = ", ".join([str(i[0]) for i in emptyfield])
                self.__output.addMessage(__file__, 4, "EBPARSE",
                        "The first and last column can't be left empty on "
                        "line(s): %s.\n" % lines)

            if notthree or emptyfield:
                return None

            #Create a Namechecker batch entry
            ret = []
            for line, job in jobl[1:]:
                if job[1]:
                    if job[0].startswith("LRG"):
                        inputl = "%s%s:%s" % tuple(job)
                    else:
                        inputl = "%s(%s):%s" % tuple(job)
                else:
                    inputl = "%s:%s" % (job[0], job[2])
                ret.append(inputl)
            return ret

        else:   #No Header, possibly a new BatchType
            #collect all lines with data in fields other than the first
            lines = ", ".join([str(row[0]) for row in jobl if any(row[1][1:])])
            if any(lines):
                self.__output.addMessage(__file__, 4, "EBPARSE",
                    "New Type Batch jobs (see help) should contain one "
                    "entry per line, please check line(s): %s" % lines)
                self.__output.addMessage(__file__, 2, "DEBUG",
                        `job`)
            else:
                return [job[0] for line, job in jobl]


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

        mimeType = self.getMimeType(handle)
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
