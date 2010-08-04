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
        #jobl = [(l+1, row) for l, row in enumerate(job) if row and any(row)]
        jobl = [(l+1, row) for l, row in enumerate(job)]

        #TODO:  Add more new style old style logic

        if jobl[0][1] == self.__config.header : #Old style NameCheckBatch job
            ret = []
            notthree = []
            emptyfield = []
            for line, job in jobl[1:]:
                if not any(job):    #Empty line
                    ret.append("~!")
                    continue

                inputl = ""
                if len(job)!=3:     #Need three columns
                    notthree.append(line)
                elif (not(job[0] and job[2])):
                    emptyfield.append(line)
                else:
                    if job[1]:
                        if job[0].startswith("LRG"):
                            inputl = "%s%s:%s" % tuple(job)
                        else:
                            inputl = "%s(%s):%s" % tuple(job)
                    else:
                        inputl = "%s:%s" % (job[0], job[2])

                if not inputl:
                    #try to make something out of it
                    inputl = "~!InputFields: " #start with the skip flag
                    inputl+= "|".join(job)

                ret.append(inputl)
            #for

            #Create output Message for incompatible fields
            if any(notthree):
                lines = makeList(notthree, 10)
                self.__output.addMessage(__file__, 3, "EBPARSE",
                        "Wrong amount of columns in %i line(s): %s.\n" %
                        (len(notthree), lines))

            if any(emptyfield):
                lines = makeList(emptyfield, 10)
                self.__output.addMessage(__file__, 3, "EBPARSE",
                        "The first and last column can't be left empty in "
                        "%i line(s): %s.\n" % (len(emptyfield), lines))

            errlist = notthree + emptyfield
        #if

        else:   #No Header, possibly a new BatchType
            if len(jobl) == 0: return
            #collect all lines with data in fields other than the first
            errlist = [line for line, row in jobl if any(row[1:])]
            if any(errlist):
                self.__output.addMessage(__file__, 3, "EBPARSE",
                    "New Type Batch jobs (see help) should contain one "
                    "entry per line, please check %i line(s): %s" %
                    (len(errList), makeList(errlist)))

            ret = []
            for line, job in jobl:
                if not any(job):    #Empty line
                    ret.append("~!")
                    continue
                if line in lines:
                    inputl = "~!InputFields: "   #Dirty Escape BatchEntries
                else:
                    inputl = ""
                ret.append(inputl+"|".join(job))
        #else

        if not ret: return None     #prevent divide by zero

        err = float(len(errlist))/len(ret)
        if err == 0:
            return ret
        elif err < 0.05:
            #allow a 5 percent threshold for errors in batchfiles
            self.__output.addMessage(__file__, 3, "EBPARSE",
                    "There were errors in your batch entry file, they are "
                    "omitted and your batch is started.")
            self.__output.addMessage(__file__, 3, "EBPARSE",
                    "Please check the batch input file help at the top of "
                    "this page for additional information.")
            return ret
        else:
            return None
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

def makeList(l, maxlen=10):
    ret = ", ".join(str(i) for i in l[:maxlen])
    if len(l)>maxlen:
        return ret+", ..."
    else:
        return ret
