#!/usr/bin/python

"""
Module for parsing CSV files and spreadsheets.

@todo: Check ODS, XLS compatibility
@requires: magic
@requires: csv
@requires: xlrd
@requires: zipfile
@requires: xml.dom.minidom
@requires: os
@requires: types
"""
# Public classes:
#     - File ; Parse CSV files and spreadsheets.


import magic           # open(), MAGIC_MIME, MAGIC_NONE
import csv             # Sniffer(), reader(), Error
import xlrd            # open_workbook()
import zipfile         # ZipFile()
import xml.dom.minidom # parseString()
import os              # remove()
import tempfile
import types           # UnicodeType
from cStringIO import StringIO

from mutalyzer import util

class File() :
    """
    Parse CSV files and spreadsheets.

    Private variables:
        - __config ; Configuration variables.
        - __output ; The Output object.

    Special methods:
        - __init__(config, output) ; Initialise the class.

    Private methods:
        - __tempFileWrapper(func, handle) ; Call func() with a filename.
        - __parseCsvFile(handle)    ; Parse a CSV file.
        - __parseXlsFile(handle)    ; Parse an Excel file.
        - __parseOdsFile(handle)    ; Parse an OpenDocument Spreadsheet file.
        - __checkBatchFormat(job)   ; Check a batch job and sanitize it.

    Public methods:
        - getMimeType(handle)    ; Get the mime type of a stream.
        - parseFileRaw(handle)   ; Parse a stream with the appropriate parser.
        - parseBatchFile(handle) ; Parse a stream with the appropriate parser
                                   and sanitize the output.
    """

    def __init__(self, config, output) :
        """
        Initialise the class.

        Private variables (altered):
            - __config ; Initialised with configuration variables.
            - __output ; Set to the Output object.
        
        @arg config: Configuration variables
        @type config: class instance
        @arg output: Output object
        @type output: class instance
        """

        self.__config = config
        self.__output = output #: The Output object
    #__init__

    def __tempFileWrapper(self, func, handle) :
        """
        Make a temporary file, put the content of a stream in it and pass
        the filename to a general function. Return whatever this function
        returns.

        @arg func: general function that needs a file name as argument
        @type func: function
        @arg handle: A stream
        @type handle: stream

        @return: unknown; the output of func().
        @rtype: ?
        """
        write_handle, filename = tempfile.mkstemp(text=True)

        # Dump the content of the stream pointed to by handle into the file.
        handle.seek(0)
        write_handle.write(handle.read())
        write_handle.close()

        # Open the file with func().
        ret = func(filename)
        # Apperantly apache will remove this file even when opened by the
        # function *func
        os.remove(filename)

        return ret
    #__tempFileWrapper

    def __parseCsvFile(self, handle) :
        """
        Parse a CSV file.
        The stream is not rewinded after use.

        Private variables:
            - __config ; The bufSize configuration variables.

        @arg handle: A handle to a stream
        @type handle: stream

        @return: list of lists
        @rtype: list
        """

        # If we naively assume the input file uses \n characters as
        # newlines, the CSV parser can trip over e.g. Windows style
        # newlines. It will probably complain with a message like:
        #
        #   new-line character seen in unquoted field - do you need to open
        #   the file in universal-newline mode?
        #
        # Here we try to support multiplatform newline modes in the input
        # file, so \n, \r, \r\n are all recognized as a newline.
        #
        # This can be done by opening the file in 'U' mode, but in this case
        # we already have an opened file (probably, if the call originated
        # from a web request, opened by the web.py input handler, which uses
        # the Python cgi module for opening uploaded files).
        #
        # A fix is to get the handle's file descriptor and create a new
        # handle for it, using 'U' mode.
        #
        # However, sometimes our handler has no .fileno(), for example when
        # the input file is quite small (< 1kb). In that case, the cgi
        # module seems to optimize things and use a StringIO for the file,
        # which of course has no .fileno().
        #
        # So our solution is:
        # - We have .fileno(): Create a new handle, using 'U' mode.
        # - We have no .fileno(): Replace all newlines by \n (in-memory)
        #   and wrap the result in a new StringIO.

        if hasattr(handle, 'fileno'):
            # Todo: We get the following error in our logs (exact origin
            # unknown):
            #
            #   close failed in file object destructor:
            #   IOError: [Errno 9] Bad file descriptor
            #
            # I am unable to find the reason for this. Everything seems to
            # be working though.
            new_handle = os.fdopen(handle.fileno(), 'rU')
        elif hasattr(handle, 'getvalue'):
            data = handle.getvalue()
            data = data.replace('\r\n', '\n').replace('\r', '\n')
            new_handle = StringIO(data)
        else:
            self.__output.addMessage(__file__, 4, "EBPARSE",
                                     "Fatal error parsing input file, please"
                                     " report this as a bug including the"
                                     " input file.")
            return None

        # I don't think the .seek(0) is needed now we created a new handle
        new_handle.seek(0)
        buf = new_handle.read(self.__config.bufSize)

        # Default dialect
        dialect = 'excel'

        # The idea is that for new-style batch input files we have only
        # one column and the sniffer cannot find a delimiter.

        try :
            # Todo: delimiters in config file
            dialect = csv.Sniffer().sniff(buf, delimiters="\t ;|,")
        except csv.Error, e :
            #self.__output.addMessage(__file__, 4, "EBPARSE", e)
            #return None
            pass
        #except

        #Watch out for : delimiter FIXME and for the . delimiter
#        if dialect.delimiter == ":":
#            dialect.delimiter = "\t"

        new_handle.seek(0)
        reader = csv.reader(new_handle, dialect)

        ret = []
        for i in reader :
            ret.append(i)

        new_handle.close()
        return ret
    #__parseCsvFile

    def __parseXlsFile(self, handle) :
        """
        Parse an Excel file.
        The stream is not rewinded after use.

        @arg handle: A handle to a stream
        @type handle: stream

        @return: A list of lists
        @rtype: list
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

        @arg handle: A handle to a stream
        @type handle: stream

        @return: A list of lists
        @rtype: list
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
           - The first and the last element should be non-empty.
           - The first line should be the header defined in the config file.

        Private variables:
            - __config ; The header configuration variable.
        
        @todo: Add more new style old style logic
        @todo: if not inputl: try to make something out of it

        @arg job: list of lists
        @type job: list

        @return: A sanitised list of lists (without a header or empty lines)
        @rtype: list
        """
        #store original line numbers line 1 = job[0]
        jobl = [(l+1, row) for l, row in enumerate(job)]

        #TODO:  Add more new style old style logic

        if jobl[0][1] == self.__config.header : #Old style NameCheckBatch job
            ret = []
            notthree = []
            emptyfield = []
            for line, job in jobl[1:]:

                #Empty line
                if not any(job):
                    ret.append("~!")
                    continue

                inputl = ""
                if len(job)!=3:     #Need three columns
                    notthree.append(line)
                elif (not(job[0] and job[2])):
                    # First and last column cant be empty
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
                    #TODO: try to make something out of it
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
                if line in errlist:
                    inputl = "~!InputFields: "   #Dirty Escape BatchEntries
                else:
                    inputl = ""
                ret.append(inputl+"|".join(job))
        #else

        if not ret: return None     #prevent divide by zero

        err = float(len(errlist))/len(ret)
        if err == 0:
            return ret
        elif err < self.__config.threshold:
            #allow a 5 (default) percent threshold for errors in batchfiles
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

    def getMimeType(self, handle) :
        """
        Get the mime type of a stream by inspecting a fixed number of bytes.
        The stream is rewinded after use.

        Private variables:
            - __config: The bufSize configuration variables.

        @arg handle: A handle to a stream
        @type handle: stream

        @return: The mime type of a file
        @rtype: string
        """

        handle.seek(0)
        buf = handle.read(self.__config.bufSize) #: The bufSize configuration variables.


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

    def parseFileRaw(self, handle) :
        """
        Check which format a stream has and parse it with the appropriate
        parser if the stream is recognised.

        @arg handle: A handle to a stream
        @type handle: stream

        @return: A list of lists, None if an error occured
        @rtype: list
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

        @arg handle: A handle to a stream
        @type handle: stream

        @return: A sanitised list of lists (without a header or empty lines),
        or None if an error occured
        @rtype: list
        """

        job = self.parseFileRaw(handle)
        if job :
            return self.__checkBatchFormat(job)
        return None
    #parseBatchFile
#File

def makeList(l, maxlen=10):
    """
    Converts a list of lines to a string to be used in output messages for
    incompatible fields.

    @arg l: list of lines
    @type l: list
    @arg maxlen: maximum length of the string you want to return
    @type maxlen: integer
    @return: a list converted to a string with comma's and spaces
    @rtype: string
    """
    ret = ", ".join(str(i) for i in l[:maxlen])
    if len(l)>maxlen:
        return ret+", ..."
    else:
        return ret
#makeList
