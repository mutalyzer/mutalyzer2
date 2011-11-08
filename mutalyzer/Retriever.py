"""
Module for retrieving files from either the cache or the NCBI.

A hash of every retrieved file is stored in the internal database. If a
requested file is not found, but its hash is, we use additional information
to re-download the file.

Public classes:
- Retriever ; Retrieve a record from either the cache or the NCBI.
"""


import os              # path.isfile(), link() path.isdir(), path.mkdir(),
                       # walk(), path.getsize(), path.join(), stat(), remove()
import bz2             # BZ2Compressor(), BZ2File()
import hashlib         # md5(), update(), hexdigest()
import urllib2         # urlopen()
import StringIO        # StringIO()
from Bio import SeqIO  # read()
from Bio import Entrez # efetch(), read(), esearch(), esummary()
from Bio.Seq import UnknownSeq
from Bio.Alphabet import ProteinAlphabet
from xml.dom import DOMException, minidom
from xml.parsers import expat

from mutalyzer import util
from mutalyzer.parsers import lrg
from mutalyzer.parsers import genbank


class Retriever(object) :
    """
    Retrieve a record from either the cache or the NCBI.

    Inherited variables from Db.Output.Config:
        - email     ; The email address which we give to the NCBI.
        - cache     ; The directory where the records are stored.
        - cachesize ; Maximum size of the cache.

    Special methods:
        - __init__(config, output, database) ; Use variables from the
        configuration file to initialise the class private variables.



    Private methods:
        - _foldersize(folder) ; Return the size of a folder.
        - _cleancache()       ; Keep the cache at a maximum size.
        - _nametofile(name)   ; Convert a name to a filename.
        - _write(raw_data, filename, extract) ; Write a record to a file.
        - _calcHash(content)  ; Calculate the md5sum of 'content'.
        - _newUD()            ; Generate a new UD number.

    Public methods:
        - retrieveslice(accno, start, stop, orientation) ; Retrieve a chromosome
        slice from the NCBI.
        - retrievegene(gene, organism, upstream, downstream) ; Retrieve a gene
        from the NCBI.
        - downloadrecord(url)    ; Download a GenBank file.
        - uploadrecord(raw_data) ; Let someone upload a GenBank file.
        - loadrecord(identifier) ; Load a record, store it in the cache, manage
        the cache and return the record.

    Inherited methods from Db.Output:
        - WarningMsg(filename, message) ; Print a warning message.
        - ErrorMsg(filename, message)   ; Print an error message and log it.
        - LogMsg(filename, message)     ; Log a message.
    """

    def __init__(self, config, output, database) :
        """
        Use variables from the configuration file for some simple
        settings. Make the cache directory if it does not exist yet.

        Inherited variables from Db.Output.Config:
            - email     ; The email address which we give to the NCBI.
            - cache     ; The directory where the records are stored.

        @arg config:
        @type config:
        @arg output:
        @type output:
        @arg database:
        @type database:
        """

        self._config = config
        self._output = output
        self._database = database
        if not os.path.isdir(self._config.cache) :
            os.mkdir(self._config.cache)
        Entrez.email = self._config.email
        self.fileType = None
    #__init__

    def _foldersize(self, folder) :
        """
        Return the size of a folder in bytes.

        @arg folder: Name of a directory
        @type folder: string

        @return: The size of the directory
        @rtype: integer
        """

        folder_size = 0
        for (path, dirs, files) in os.walk(folder) :
            for fileName in files :
                folder_size += os.path.getsize(os.path.join(path, fileName))

        return folder_size
    #_foldersize

    def _cleancache(self) :
        """
        Keep removing files until the size of the cache is less than the
        maximum size.
        First, the cache checked for its size, if it exceeds the maximum
        size the ``oldest'' files are deleted. Note that accessing a file
        makes it ``new''.

        Inherited variables from Db.Output.Config:
            - cache     ; Directory under scrutiny.
            - cachesize ; Maximum size of the cache.
        """

        if self._foldersize(self._config.cache) < self._config.cachesize:
            return

        # Build a list of files sorted by access time.
        cachelist = []
        for (path, dirs, files) in os.walk(self._config.cache) :
            for filename in files :
                filepath = os.path.join(path, filename)
                cachelist.append(
                    (os.stat(filepath).st_atime, filepath))
        cachelist.sort()

        # Now start removing pairs of files until the size of the folder is
        # small enough (or until the list is exhausted).
        for i in range(0, len(cachelist)) :
            os.remove(cachelist[i][1])
            if self._foldersize(self._config.cache) < self._config.cachesize:
                break;
        #for
    #_cleancache

    def _nametofile(self, name) :
        """
        Convert an accession number to a filename.

        Inherited variables from Db.Output.Config:
            - cache     ; Name of the cache directory.

        @arg name: The accession number
        @type name: string

        @return: A filename
        @rtype: string
        """

        return self._config.cache + '/' + name + "." + self.fileType + ".bz2"
    #_nametofile

    def _write(self, raw_data, filename) :
        """
        Write raw data to a compressed file.

        @arg raw_data: The raw_data to be compressed and written
        @type raw_data: string
        @arg filename: The intended name of the outfile
        @type filename: string

        @return: outfile ; The full path and name of the file written
        @rtype: string
        """
        # Compress the data to save disk space.
        comp = bz2.BZ2Compressor()
        data = comp.compress(raw_data)
        data += comp.flush()
        out_handle = open(self._nametofile(filename), "w")
        out_handle.write(data)
        out_handle.close()

        # Since we put something in the cache, check if it needs cleaning.
        self._cleancache()

        return out_handle.name      # return the full path to the file
    #_write

    def _calcHash(self, content) :
        """
        Calculate the md5sum of a piece of text.

        @arg content: Arbitrary text
        @type content: string

        @return: The md5sum of 'content'
        @rtype: string
        """

        hashfunc = hashlib.md5()
        hashfunc.update(content)
        md5sum = hashfunc.hexdigest()
        del hashfunc

        return md5sum
    #_calcHash

    def _newUD(self) :
        """
        Make a new UD number based on the current time (seconds since 1970).

        @return: A new UD number
        @rtype: string
        """

        UD = util.generate_id()
        return "UD_" + str(UD)
    #_newUD

    def _updateDBmd5(self, raw_data, name, GI):
        #TODO documentation
        """
        @todo: documentation

        @arg raw_data:
        @type raw_data:
        @arg name:
        @type name:
        @arg GI:
        @type GI:

        @return: filename
        @rtype: string
        """

        currentmd5sum = self._database.getHash(name)
        if currentmd5sum :
            md5sum = self._calcHash(raw_data)
            if md5sum != currentmd5sum :
                self._output.addMessage(__file__, -1, "WHASH",
                    "Warning: Hash of %s changed from %s to %s." % (
                    name, currentmd5sum, md5sum))
                self._database.updateHash(name, md5sum)
            #if
        else :
            self._database.insertGB(name, GI,
                self._calcHash(raw_data), None, 0, 0, 0, None)
        return self._nametofile(name)
    #_updateDBmd5


    def snpConvert(self, rs_id) :
        """
        Search for an rsId in dbSNP and return all annotated HGVS notations of
        it.

        @arg rsId: The rsId of the SNP (example: 'rs9919552').
        @type rsId: string

        @return: A list of HGVS notations.
        @rtype: list(string)
        """
        # A simple input check.
        id = rs_id[2:]
        if rs_id[:2] != 'rs' or not id.isdigit():
            self._output.addMessage(__file__, 4, 'ESNPID',
                                    'This is not a valid dbSNP id.')
            return []

        # Query dbSNP for the SNP.
        try:
            response = Entrez.efetch(db='SNP', id=id, rettype='flt',
                                     retmode='xml')
        except IOError as e:
            # Could not parse XML.
            self._output.addMessage(__file__, 4, 'EENTREZ',
                                    'Error connecting to dbSNP.')
            self._output.addMessage(__file__, -1, 'INFO',
                                    'IOError: %s' % str(e))
            return []

        response_text = response.read()

        try:
            # Parse the output.
            doc = minidom.parseString(response_text)
            exchange_set = doc.getElementsByTagName('ExchangeSet')
            rs = exchange_set[0].getElementsByTagName('Rs')
        except expat.ExpatError as e:
            # Could not parse XML.
            self._output.addMessage(__file__, 4, 'EENTREZ', 'Unknown dbSNP ' \
                                    'error. Error parsing result XML.')
            self._output.addMessage(__file__, -1, 'INFO',
                                    'ExpatError: %s' % str(e))
            self._output.addMessage(__file__, -1, 'INFO',
                                    'Result from dbSNP: %s' % response_text)
            return []
        except IndexError:
            # The expected root element is not present.
            self._output.addMessage(__file__, 4, 'EENTREZ', 'Unkown dbSNP ' \
                                    'error. Result XML was not as expected.')
            self._output.addMessage(__file__, -1, 'INFO',
                                    'Result from dbSNP: %s' % response_text)
            return []

        if len(rs) < 1:
            # No Rs result element.
            text = []
            for node in exchange_set[0].childNodes:
                if node.nodeType == node.TEXT_NODE:
                    text.append(node.data)
            message = ''.join(text)
            if message.find('cannot get document summary') != -1:
                # Entrez does not have this rs ID.
                self._output.addMessage(__file__, 4, 'EENTREZ',
                                        'ID rs%s could not be found in dbSNP.' \
                                        % id)
            else:
                # Something else was wrong (print {message} to see more).
                self._output.addMessage(__file__, 4, 'EENTREZ',
                                        'Unkown dbSNP error. Got no result ' \
                                        'from dbSNP.')
                self._output.addMessage(__file__, -1, 'INFO',
                                        'Message from dbSNP: %s' % message)
            return []

        snps = []
        for i in rs[0].getElementsByTagName('hgvs'):
            snps.append(i.lastChild.data.encode('utf8'))

        return snps
    #snpConvert
#Retriever

class GenBankRetriever(Retriever):
    # TODO documentation
    """
    """

    def __init__(self, config, output, database):
        # TODO documentation
        """
        """

        # Recall init of parent
        Retriever.__init__(self, config, output, database)
        self.fileType = "gb"
        # Child specific init
    #__init__

    def write(self, raw_data, filename, extract) :
        """
        Write raw data to a file. The data is parsed before writing, if a
        parse error occurs an error is returned and the function exits.
        If 'filename' is set and 'extract' is set to 0, then 'filename' is
        used for output.
        If 'extract' is set to 1, then the filename is constructed from the
        id of the GenBank record. Additionally the id and GI number are
        returned for further processing (putting them in the internal
        database).

        @arg raw_data: The data
        @type raw_data: string
        @arg filename: The intended name of the file.
        @type filename: string
        @arg extract: Flag that indicates whether to extract the record ID and
        GI number:
            - 0 ; Do not extract, use 'filename'
            - 1 ; Extract
        @type extract: integer

        @return: tuple ; Depending on the value of 'extract':
            - 0 ; ('filename', None)
            - 1 ; (id, GI)
        @rtype: tuple (string, string)
        """

        if raw_data == "\nNothing has been found\n" :
            self._output.addMessage(__file__, 4, "ENORECORD",
                "The record could not be retrieved.")
            return None
        #if

        fakehandle = StringIO.StringIO() # Unfortunately, BioPython needs a
        fakehandle.write(raw_data)       # file handle.
        fakehandle.seek(0)
        try :
            record = SeqIO.read(fakehandle, "genbank")
        except (ValueError, AttributeError):  # An error occured while parsing.
            self._output.addMessage(__file__, 4, "ENOPARSE",
                "The file could not be parsed.")
            fakehandle.close()
            return None
        #except

        if type(record.seq) == UnknownSeq :
            fakehandle.close()
            self._output.addMessage(__file__, 4, "ENOSEQ",
                "This record contains no sequence. Chromosomal or contig " \
                "records should be uploaded with the GenBank uploader.")
            return None
        #if

        outfile = filename
        GI = None
        if extract :
            outfile = record.id
            GI = record.annotations["gi"]
            if outfile != filename :
                # Add the reference (incl version) to the reference output
                # This differs if the original reference lacks a version
                self._output.addOutput("reference", record.id)
                self._output.addOutput(
                        "BatchFlags", ("A1",(
                            filename,
                            outfile,
                            filename+"[[.period.]]" )))
                self._output.addMessage(__file__, 2, "WNOVER",
                    "No version number is given, using %s. Please use this " \
                    "number to reduce downloading overhead." % record.id)
        #if
        fakehandle.close()

        self._write(raw_data, outfile)

        return outfile, GI
    #write

    def fetch(self, name) :
        """
        Todo: Documentation.
        """
        net_handle = Entrez.efetch(db='nucleotide', id=name, rettype='gb')
        raw_data = net_handle.read()
        net_handle.close()

        if raw_data == '\n' :       # Check if the file is empty or not.
            self._output.addMessage(__file__, 4, 'ERETR',
                                    'Could not retrieve %s.' % name)
            return None

        # This is a hack to detect constructed references, the proper way to
        # do this would be to check the data_file_division attribute of the
        # parsed GenBank file (it would be 'CON').
        if '\nCONTIG' in raw_data:
            try:
                # Get the length in base pairs
                length = int(raw_data[:raw_data.index(' bp', 0, 500)].split()[-1])
            except ValueError, IndexError:
                self._output.addMessage(__file__, 4, 'ERETR',
                                        'Could not retrieve %s.' % name)
                return None
            if length > self._config.maxDldSize:
                self._output.addMessage(__file__, 4, 'ERETR',
                                        'Could not retrieve %s.' % name)
                return None
            net_handle = Entrez.efetch(db='nucleotide', id=name, rettype='gbwithparts')
            raw_data = net_handle.read()
            net_handle.close()

        result = self.write(raw_data, name, 1)
        if not result:
            return None
        name, GI = result
        if name:               # Processing went okay.
            return self._updateDBmd5(raw_data, name, GI)
        else:                  # Parse error in the GenBank file.
            return None
    #fetch

    def retrieveslice(self, accno, start, stop, orientation) :
        """
        Retrieve a slice of a chromosome.
        If the arguments are recognised (found in the internal database),
        we look if the associated file is still present and if so: return
        its UD number.
        If the arguments are recognised but no file was found, we download
        the new slice and update the hash (and log if the hash changes).
        If the arguments are not recognised, we download the new slice and
        make a new UD number.
        The content of the slice is placed in the cache with the UD number
        as filename.

        Inherited variables from Db.Output.Config:
            - maxDldSize ; Maximum size of the slice.

        @arg accno: The accession number of the chromosome
        @type accno: string
        @arg start: Start position of the slice
        @type start: integer
        @arg stop: End position of the slice.
        @type stop: integer
        @arg orientation:
        Orientation of the slice:
            - 1 ; Forward
            - 2 ; Reverse complement
        @type orientation: integer

        @return: An UD number
        @rtype: string
        """

        # Not a valid slice.
        if start >= stop :
            return None

        # The slice can not be too big.
        if stop - start > self._config.maxDldSize :
            return None

        # Check whether we have seen this slice before.
        UD = self._database.getGBFromLoc(accno, start, stop, orientation)
        if UD : # This has been requested before.
            if os.path.isfile(self._nametofile(UD)) : # It's still present.
                return UD

        # It's not present, so download it.
        handle = Entrez.efetch(db = "nucleotide", rettype = "gb",
                               id = accno, seq_start = start, seq_stop = stop,
                               strand = orientation)
        raw_data = handle.read()
        handle.close()

        # Calculate the hash of the downloaded file.
        md5sum = self._calcHash(raw_data)

        if UD : # We have seen this one before.
            currentmd5sum = self._database.getHash(UD)
            if md5sum != currentmd5sum :
                self._output.addMessage(__file__, -1, "WHASH",
                    "Warning: Hash of %s changed from %s to %s." % (
                    UD, currentmd5sum, md5sum))
                self._database.updateHash(UD, md5sum)
            #if
        else : # We haven't seen it before, so give it a name.
            UD = self._newUD()
            self._database.insertGB(UD, None, md5sum, accno, start,
                          stop, orientation, None)
        #else

        return self.write(raw_data, UD, 0) and UD
    #retrieveslice

    def retrievegene(self, gene, organism, upstream, downstream) :
        """
        Query the NCBI for the chromosomal location of a gene and make a
        slice if the gene can be found.

        @arg gene: Name of the gene
        @type gene: string
        @arg organism: The organism in which we search.
        @type organism: string
        @arg upstream: Number of upstream nucleotides for the slice.
        @type upstream: integer
        @arg downstream: Number of downstream nucleotides for the slice.
        @type downstream: integer

        @return: slice
        @rtype:
        """

        # Search the NCBI for a specific gene in an organism.
        query = "%s[Gene] AND %s[Orgn]" % (gene, organism)
        handle = Entrez.esearch(db = "gene", term = query)
        searchresult = Entrez.read(handle)
        handle.close()

        ChrAccVer = None        # We did not find anything yet.
        aliases = []            # A list of aliases in case we find them.
        for i in searchresult["IdList"] :                 # Inspect all results.
            handle = Entrez.esummary(db = "gene", id = i)
            summary = Entrez.read(handle)
            handle.close()
            if summary[0]["NomenclatureSymbol"] == gene : # Found it.
                if not summary[0]["GenomicInfo"] :
                    self._output.addMessage(__file__, 4, "ENOMAPPING",
                        "No mapping information found for gene %s." % gene)
                    return None
                #if
                ChrAccVer = summary[0]["GenomicInfo"][0]["ChrAccVer"]
                ChrLoc = summary[0]["GenomicInfo"][0]["ChrLoc"]
                ChrStart = summary[0]["GenomicInfo"][0]["ChrStart"]
                ChrStop = summary[0]["GenomicInfo"][0]["ChrStop"]
                break;
            #if

            # Collect official symbols that has this gene as alias in case we
            # can not find anything.
            if gene in summary[0]["OtherAliases"] and \
                summary[0]["NomenclatureSymbol"] :
                aliases.append(summary[0]["NomenclatureSymbol"]);
        #for

        if not ChrAccVer : # We did not find any genes.
            if aliases :
                self._output.addMessage(__file__, 4, "ENOGENE",
                    "Gene %s not found, found aliases: %s" % (gene, aliases))
                return None
            #if
            self._output.addMessage(__file__, 4, "ENOGENE",
                "Gene %s not found." % gene)
            return None
        #if

        # Figure out the orientation of the gene.
        orientation = "1"
        if ChrStart > ChrStop :             # Swap start and stop.
            orientation = "2"
            temp = ChrStart
            ChrStart = ChrStop - downstream # Also take care of the flanking
            ChrStop = temp + upstream + 1   # sequences.
        #if
        else :
            ChrStart -= upstream - 1
            ChrStop += downstream + 2
        #else

        # And retrieve the slice.
        return self.retrieveslice(ChrAccVer, ChrStart, ChrStop, orientation)
    #retrievegene

    def downloadrecord(self, url) :
        """
        Download a GenBank record from a URL.
        If the downloaded file is recognised by its hash, the old UD number
        is used.

        Inherited variables from Db.Output.Config:
            - maxDldSize ; Maximum size of the file.
            - minDldSize ; Minimum size of the file.

        @arg url: Location of a GenBank record
        @type url: string

        @return: UD or None
        @rtype: string
        """

        handle = urllib2.urlopen(url)
        info = handle.info()
        if info["Content-Type"] == "text/plain" :
            length = int(info["Content-Length"])
            if length > self._config.minDldSize and \
               length < self._config.maxDldSize :
                raw_data = handle.read()
                md5sum = self._calcHash(raw_data)
                UD = self._database.getGBFromHash(md5sum)
                if UD:  #hash found
                    if not os.path.isfile(self._nametofile(UD)):
                        UD = self.write(raw_data, UD, 0) and UD
                else:
                    UD = self._newUD()
                    if not os.path.isfile(self._nametofile(UD)):
                        UD = self.write(raw_data, UD, 0) and UD
                    if UD:      #Parsing went OK, add to DB
                        self._database.insertGB(UD, None, md5sum,
                                None, 0, 0, 0, url)
                return UD #Returns the UD or None
            #if
            else :
                self._output.addMessage(__file__, 4, "EFILESIZE",
                    "Filesize is not within the allowed boundaries.")
                return None
            #else
        #if
        else :
            self._output.addMessage(__file__, 4, "ERECPARSE",
                                     "This is not a GenBank record.")
            return None
        #else
    #downloadrecord

    def uploadrecord(self, raw_data) :
        """
        Write an uploaded record to a file.
        If the downloaded file is recognised by its hash, the old UD number
        is used.

        @arg raw_data: A GenBank record
        @type raw_data: string

        @return:
        @rtype: string?????
        """
        md5sum = self._calcHash(raw_data)
        UD = self._database.getGBFromHash(md5sum)
        if not UD :
            UD = self._newUD()
            if self.write(raw_data, UD, 0):
                self._database.insertGB(UD, None, md5sum, None, 0, 0, 0, None)
                return UD
        #if
        else:
            if os.path.isfile(self._nametofile(UD)):
                return UD
            else:
                return self.write(raw_data, UD, 0) and UD
    #uploadrecord

    def loadrecord(self, identifier) :
        """
        Load a record and return it.
        If the filename associated with the accession number is not found
        in the cache, try to re-download it.

        @arg identifier: An accession number
        @type identifier: string

        @return: A GenBank.Record record
        @rtype: object
        """
        if (identifier[0].isdigit()) : # This is a GI identifier.
            name = self._database.getGBFromGI(identifier)
        else :
            name = identifier

        # Make a filename based upon the identifier.
        filename = self._nametofile(name)

        if not os.path.isfile(filename) :   # We can't find the file.
            md5 = self._database.getHash(name)
            if md5:   # We have seen it before though.
                Loc = self._database.getLoc(name)  # Try to find the location.
                if not Loc[0]:              # No location found.
                    url = self._database.getUrl(name)   # Try to find an URL.
                    if not url :
                        if self._database.getGI(name) : # It was from NCBI.
                            filename = self.fetch(name)
                        else :
                            self._output.addMessage(__file__, 4, "ERETR",
                                "Please upload this sequence again.")
                            filename = None
                    #if
                    else :                  # This used to be a downloaded seq
                        filename = self.downloadrecord(url) and filename
                #if
                else :                      # This used to be a slice.
                    filename = self.retrieveslice(*Loc) and filename
            #if
            else :                          # Never seen this name before.
                filename = self.fetch(name)
            #else
        #if

        # If filename is None an error occured
        if filename is None:
            #Notify batch to skip all instance of identifier
            self._output.addOutput("BatchFlags", ("S1", identifier))
            return None

        # Now we have the file, so we can parse it.
        GenBankParser = genbank.GBparser()
        record = GenBankParser.create_record(filename)

        # Todo: This will change once we support protein references
        if isinstance(record.seq.alphabet, ProteinAlphabet):
            self._output.addMessage(__file__, 4, 'EPROTEINREF',
                                    'Protein reference sequences are not supported.')
            return None

        return record
    #loadrecord
#GenBankRetriever

class LRGRetriever(Retriever):
    """
    Retrieve a LRG record from either the cache or the web.

    Public methods:
        - loadrecord(identifier) ; Load a record, store it in the cache, manage
                                   the cache and return the record.
    """

    def __init__(self, config, output, database):
        #TODO documentation
        """
        Initialize the class.

        @todo: documentation
        @arg  config:
        @type  config:
        @arg  output:
        @type  output:
        @arg  database:
        @type  database:
        """

        # Recall init of parent
        Retriever.__init__(self, config, output, database)
        self.fileType = "xml"
        # Child specific init
    #__init__

    def loadrecord(self, identifier):
        """
        Load and parse a LRG file based on the identifier

        @arg identifier: The name of the LRG file to read
        @type identifier: string

        @return: record ; GenRecord.Record of LRG file
                   None ; in case of failure
        @rtype:
        """

        # Make a filename based upon the identifier.
        filename = self._nametofile(identifier)

        if not os.path.isfile(filename) :   # We can't find the file.
            filename = self.fetch(identifier)

        if filename is None:                # return None in case of error
            #Notify batch to skip all instance of identifier
            self._output.addOutput("BatchFlags", ("S1", identifier))
            return None

        # Now we have the file, so we can parse it.
        file_handle = bz2.BZ2File(filename, "r")

        #create GenRecord.Record from LRG file
        record = lrg.create_record(file_handle.read())
        file_handle.close()

        return record
    #loadrecord

    def fetch(self, name):
        """
        Fetch the LRG file and store in the cache directory. First try to
        grab the file from the confirmed section, if this fails, get it
        from the pending section.

        Inherited variables from Config.Retriever
            - lrgURL  ; The base url from where LRG files are fetched

        @arg name: The name of the LRG file to fetch
        @type name: string

        @return: the full path to the file; None in case of an error
        @rtype: string
        """

        prefix = self._config.lrgURL
        url        = prefix + "%s.xml"          % name
        pendingurl = prefix + "pending/%s.xml"  % name

        try:
            return self.downloadrecord(url, name)
        except urllib2.URLError: #Catch error: file not found
            pass

        try:                # Try to get the file from the pending section
            filename = self.downloadrecord(pendingurl, name)
            self._output.addMessage(__file__, 2, "WPEND",
                "Warning: LRG file %s is a pending entry." % name)
            return filename
        except urllib2.URLError:
            self._output.addMessage(__file__, 4, "ERETR",
                                 "Could not retrieve %s." % name)
            return None             #Explicit return in case of an Error
    #fetch

    def downloadrecord(self, url, name = None) :
        """
        Download an LRG record from an URL.

        Inherited variables from Db.Output.Config:
            - maxDldSize  ; Maximum size of the file.
            - minDldSize  ; Minimum size of the file.

        @arg url: Location of the LRG record
        @type url: string

        @return:
            - filename    ; The full path to the file
            - None        ; in case of failure
        @rtype: string
        """

        lrgID = name or os.path.splitext(os.path.split(url)[1])[0]
        #if not lrgID.startswith("LRG"):
        #    return None
        filename = self._nametofile(lrgID)

        handle = urllib2.urlopen(url)
        info = handle.info()
        if info["Content-Type"] == "application/xml" and info.has_key("Content-length"):

            length = int(info["Content-Length"])
            if self._config.minDldSize < length < self._config.maxDldSize:
                raw_data = handle.read()
                handle.close()

                #Do an md5 check
                md5sum = self._calcHash(raw_data)
                md5db = self._database.getHash(lrgID)
                if md5db is None:
                    self._database.insertLRG(lrgID, md5sum, url)
                elif md5sum != md5db:       #hash has changed for the LRG ID
                    self._output.addMessage(__file__, -1, "WHASH",
                        "Warning: Hash of %s changed from %s to %s." % (
                        lrgID, md5db, md5sum))
                    self._database.updateHash(lrgID, md5sum)
                else:                       #hash the same as in db
                    pass

                if not os.path.isfile(filename) :
                    return self.write(raw_data, lrgID)
                else:
                    # This can only occur if synchronus calls to mutalyzer are
                    # made to recover a file that did not exist. Still leaves
					# a window in between the check and the write.
                    return filename
            #if
            else :
                self._output.addMessage(__file__, 4, "EFILESIZE",
                    "Filesize is not within the allowed boundaries.")
        #if
        else :
            self._output.addMessage(__file__, 4, "ERECPARSE",
                                     "This is not an LRG record.")
        handle.close()
    #downloadrecord

    def write(self, raw_data, filename) :
        """
        Write raw LRG data to a file. The data is parsed before writing,
        if a parse error occurs None is returned.

        @arg raw_data: The data
        @type raw_data: string
        @arg filename: The intended name of the file
        @type filename: string

        @return:
            - filename ; The full path and name of the file written
            - None     ; In case of an error
        @rtype: string
        """
        # Dirty way to test if a file is valid,
        # Parse the file to see if it's a real LRG file.
        try:
            lrg.create_record(raw_data)
        except DOMException:
            self._output.addMessage(__file__, 4, "ERECPARSE",
                                      "Could not parse file.")
            return None             # Explicit return on Error

        return self._write(raw_data, filename) #returns full path
    #write
#LargeRetriever

if __name__ == "__main__" :
    pass
#if
