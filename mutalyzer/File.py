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


from __future__ import unicode_literals

import codecs
import re
import magic           # open(), MAGIC_MIME, MAGIC_NONE
import csv             # Sniffer(), reader(), Error
import xlrd            # open_workbook()
import zipfile         # ZipFile()
import xml.dom.minidom # parseString()
import chardet

from mutalyzer.config import settings


# Amount of bytes to be read from a file at a time (this is also the amount
# read for determining the file type).
BUFFER_SIZE = 32768


class _UniversalNewlinesByteStreamIter(object):
    """
    The codecs module doesn't provide universal newline support. This class is
    used as a stream wrapper that provides this functionality.

    The wrapped stream must yield byte strings. We decode it using the given
    encoding, normalise newlines, and yield UTF-8 encoded data (read method)
    or lines (as iterator).

    Adaptation from an old Cython version:
    https://github.com/cython/cython/blob/076fac3/Cython/Utils.py
    """
    normalise_newlines = re.compile('\r\n?|\n').sub

    def __init__(self, stream, encoding='utf-8', buffer_size=0x1000):
        # let's assume .read() doesn't change
        self.stream = codecs.getreader(encoding)(stream)
        self._read = self.stream.read
        self.buffer_size = buffer_size

    def _read_normalised(self, count=None):
        count = count or self.buffer_size
        data = self._read(count)
        if '\r' not in data:
            return data
        if data.endswith('\r'):
            # may be missing a '\n'
            data += self._read(1)
        return self.normalise_newlines('\n', data)

    def _readlines(self):
        buffer = []
        data = self._read_normalised()
        while data:
            buffer.append(data)
            lines = ''.join(buffer).splitlines(True)
            for line in lines[:-1]:
                yield line
            buffer = [lines[-1]]
            data = self._read_normalised()

        if buffer[0]:
            yield buffer[0]

    def seek(self, pos):
        if pos == 0:
            self.stream.seek(0)
        else:
            raise NotImplementedError

    def read(self, count=-1):
        return self._read_normalised(count).encode('utf-8')

    def __iter__(self):
        return (line.encode('utf-8') for line in self._readlines())


class File() :
    """
    Parse CSV files and spreadsheets.

    Private variables:
        - __output ; The Output object.

    Special methods:
        - __init__(config, output) ; Initialise the class.

    Private methods:
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

    def __init__(self, output) :
        """
        Initialise the class.

        Private variables (altered):
            - __output ; Set to the Output object.

        @arg output: Output object
        @type output: class instance
        """
        self.__output = output #: The Output object
    #__init__

    def __parseCsvFile(self, handle) :
        """
        Parse a CSV file. Does not reset the file handle to start.

        @arg handle: CSV file. Must be a seekable binary file object.
        @type handle: file object

        @return: list of lists
        @rtype: list
        """
        buf = handle.read(BUFFER_SIZE)
        result = chardet.detect(buf)
        handle.seek(0)

        if result['confidence'] > 0.5:
            encoding = unicode(result['encoding'])
        else:
            encoding = 'utf-8'

        # Python 2.7 makes it extraordinarily hard to do this correctly. We
        # have a binary file object containing lines of text in a certain
        # encoding with unknown style of line-endings.
        #
        # We want to correctly decode the file contents, accept any style of
        # line-endings, parse the lines with the `csv` module, and return
        # unicode strings.
        #
        # 1. `codecs.getreader` does not have a universal newlines mode.
        # 2. `io.TextIOWrapper` cannot be wrapped around our file object,
        #    since it is required to be an `io.BufferedIOBase`, which it
        #    usually will not be.
        # 3. The `csv` module cannot read unicode.
        #
        # Ugh.
        #
        # So, we use a stream wrapper that consumes byte strings, decodes to
        # unicode, normalises newlines, and produces the result UTF-8 encoded.
        # That's what we feed the `csv` module. We decode what it gives back
        # to unicode strings. What a mess.
        handle = _UniversalNewlinesByteStreamIter(handle, encoding=encoding,
                                                  buffer_size=BUFFER_SIZE)

        try:
            buf = handle.read(BUFFER_SIZE)
        except UnicodeDecodeError:
            self.__output.addMessage(__file__, 3, 'EBPARSE',
                                     'Could not decode file (using %s encoding).'
                                     % encoding)
            return None

        # Default dialect
        dialect = 'excel'

        # The idea is that for new-style batch input files we have only
        # one column and the sniffer cannot find a delimiter.

        try :
            # Todo: delimiters in config file
            dialect = csv.Sniffer().sniff(buf, delimiters="\t ;|,")
            dialect.skipinitialspace = True
        except csv.Error:
            #self.__output.addMessage(__file__, 4, "EBPARSE", e)
            #return None
            pass
        #except

        #Watch out for : delimiter FIXME and for the . delimiter
#        if dialect.delimiter == ":":
#            dialect.delimiter = "\t"

        handle.seek(0)
        reader = csv.reader(handle, dialect)

        ret = []
        try:
            for i in reader:
                ret.append([c.decode('utf-8') for c in i])
        except UnicodeDecodeError:
            self.__output.addMessage(__file__, 3, 'EBPARSE',
                                     'Could not decode file (using %s encoding).'
                                     % encoding)
            return None

        return ret
    #__parseCsvFile

    def __parseXlsFile(self, handle) :
        """
        Parse an Excel file. Does not reset the file handle to start.

        @arg handle: Excel file. Must be a binary file object.
        @type handle: file object

        @return: A list of lists
        @rtype: list
        """

        try:
            workBook = xlrd.open_workbook(file_contents=handle.read())
        except xlrd.XLRDError:
            return None

        sheet = workBook.sheet_by_index(0)

        ret = []
        for i in range(sheet.nrows) :
            row = []
            for j in sheet.row_values(i) :
                row.append(j)
            #for
            ret.append(row)
        #for

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
        # Todo: Use a library for this.

        zipFile = zipfile.ZipFile(handle)
        doc = xml.dom.minidom.parseString(zipFile.read("content.xml"))
        zipFile.close()

        ret = []
        for i in doc.getElementsByTagName("table:table-row") :
            row = []
            for j in i.getElementsByTagName("table:table-cell") :
                c = j.getElementsByTagName("text:p")
                if c :
                    row.append(c[0].lastChild.data)
                #if
            #for
            if row:
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

        @todo: Add more new style old style logic
        @todo: if not inputl: try to make something out of it

        @arg job: list of lists
        @type job: list

        @return: A sanitised list of lists (without a header or empty lines)
                 and the number of columns.
        @rtype: tuple(list, int)
        """
        columns = 1
        max_column_length = 200

        #store original line numbers line 1 = job[0]
        jobl = [(l+1, row) for l, row in enumerate(job)]

        #TODO:  Add more new style old style logic
        if jobl[0][1] == ['AccNo', 'Genesymbol', 'Mutation']: #Old style NameCheckBatch job
            ret = []
            notthree = []
            emptyfield = []
            toolong = []
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

                if len(inputl) > max_column_length:
                    toolong.append(line)
                else:
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

            if any(toolong):
                lines = makeList(toolong, 10)
                self.__output.addMessage(__file__, 3, "EBPARSE",
                        "Batch input field exceeds %d characters in %i line(s): %s.\n" %
                        (max_column_length, len(toolong), lines))

            errlist = notthree + emptyfield + toolong
        #if

        else:   #No Header, possibly a new BatchType
            if len(jobl) == 0:
                return (None, columns)
            # Determine number of columns from first line.
            columns = len(jobl[0][1])
            # Collect all lines with a different number of columns
            errlist = [line for line, row in jobl
                       if any(row) and len(row) != columns]
            if any(errlist):
                self.__output.addMessage(__file__, 3, "EBPARSE",
                    "New Type Batch jobs (see help) should contain the same "
                    "number of columns on every line, please check %i "
                    "line(s): %s" %
                    (len(errlist), makeList(errlist)))

            toolong = [line for line, row in jobl
                       if any(len(col) > max_column_length for col in row)]
            if any(toolong):
                self.__output.addMessage(__file__, 3, "EBPARSE",
                    "Batch input field exceeds %d characters in %i line(s): %s" %
                    (max_column_length, len(toolong), makeList(toolong)))

            ret = []
            for line, job in jobl:
                if not any(job):    #Empty line
                    ret.extend(['~!' for _ in range(columns)])
                    continue
                if line in toolong:
                    #Trim too long
                    ret.append("~!InputFields: " + ('|'.join(job))[:180] + '...')
                    ret.extend(['~!' for _ in range(columns - 1)])
                elif line in errlist:
                    #Dirty Escape BatchEntries
                    ret.append("~!InputFields: " + '|'.join(job))
                    ret.extend(['~!' for _ in range(columns - 1)])
                else:
                    ret.extend([j or '~!' for j in job])
        #else

        if not ret:
            #prevent divide by zero
            return (None, columns)

        err = float(len(errlist))/len(ret)
        if err == 0:
            return (ret, columns)
        elif err < settings.BATCH_JOBS_ERROR_THRESHOLD:
            #allow a 5 (default) percent threshold for errors in batchfiles
            self.__output.addMessage(__file__, 3, "EBPARSE",
                    "There were errors in your batch entry file, they are "
                    "omitted and your batch is started. Please check the "
                    "batch input file help at the top of this page for "
                    "additional information.")
            return (ret, columns)
        else:
            return (None, columns)
    #__checkBatchFormat

    def getMimeType(self, handle) :
        """
        Get the mime type of a stream by inspecting a fixed number of bytes.
        The stream is rewinded after use.

        @arg handle: Stream to be inspected. Must be a seekable binary file
          object.
        @type handle: file object

        @return: The mime type of a file and a textual description.
        @rtype: unicode, unicode
        """
        handle.seek(0)
        buf = handle.read(BUFFER_SIZE)

        MagicInstance = magic.open(magic.MAGIC_MIME)
        MagicInstance.load()
        mimeType = MagicInstance.buffer(buf).decode('utf-8').split(';')[0]
        MagicInstance.close()
        MagicInstance = magic.open(magic.MAGIC_NONE)
        MagicInstance.load()
        description = MagicInstance.buffer(buf).decode('utf-8')
        handle.seek(0)

        return mimeType, description
    #getMimeType

    def parseFileRaw(self, handle) :
        """
        Check which format a stream has and parse it with the appropriate
        parser if the stream is recognised. Does not reset the file handle to
        start.

        @arg handle: Input file to be parsed. Must be a seekable binary file
          object.
        @type handle: file object

        @return: A list of lists, None if an error occured
        @rtype: list
        """

        mimeType = self.getMimeType(handle)

        if mimeType[0] == "text/plain":
            return self.__parseCsvFile(handle)
        if (mimeType[0] in ('application/vnd.ms-excel',
                            'application/vnd.ms-office',
                            'application/msword',
                            'application/zip') or
            mimeType[1] == 'Microsoft OOXML'):
            return self.__parseXlsFile(handle)
        if (mimeType[0] == 'application/vnd.oasis.opendocument.spreadsheet' or
            mimeType[1] in ('OpenDocument Spreadsheet',
                            'OpenOffice.org 1.x Calc spreadsheet')):
            return self.__parseOdsFile(handle)

        return None
    #parseFileRaw

    def parseBatchFile(self, handle) :
        """
        Check which format a stream has and parse it with the appropriate
        parser if the stream is recognised. Does not reset the file handle to
        start.

        @arg handle: Batch job input file. Must be a seekable binary file
          object.
        @type handle: file object

        @return: A sanitised list of lists (without a header or empty lines)
                 (or None if an error occured) and the number of columns.
        @rtype: tuple(list, int)
        """

        job = self.parseFileRaw(handle)
        if job:
            return self.__checkBatchFormat(job)
        return (None, 1)
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
    @rtype: unicode
    """
    ret = ", ".join(str(i) for i in l[:maxlen])
    if len(l)>maxlen:
        return ret+", ..."
    else:
        return ret
#makeList
