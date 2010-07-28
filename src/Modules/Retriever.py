#!/usr/bin/python

"""
    Module for retrieving files from either the cache or the NCBI.

    A hash of every retrieved file is stored in the internal database. If a
    requested file is not found, but its hash is, we use additional information
    to re-download the file.

    Public classes:
        Retriever ; Retrieve a record from either the cache or the NCBI.
"""

import os              # path.isfile(), link() path.isdir(), path.mkdir(),
                       # walk(), path.getsize(), path.join(), stat(), remove()
import bz2             # BZ2Compressor(), BZ2File()
import hashlib         # md5(), update(), hexdigest()
import urllib2         # urlopen()
import StringIO        # StringIO()
import ftplib          # FTP(), all_errors
from Bio import SeqIO  # read()
from Bio import Entrez # efetch(), read(), esearch(), esummary()
from Bio.Seq import UnknownSeq

from Modules import Misc
from Modules import LRGparser
from Modules import GBparser
from xml.dom import DOMException
import xml.dom.minidom

class Retriever(object) :
    """
        Retrieve a record from either the cache or the NCBI.

        Inherited variables from Db.Output.Config:
            email     ; The email address which we give to the NCBI.
            cache     ; The directory where the records are stored.
            cachesize ; Maximum size of the cache.

        Special methods:
            __init__(config,   ; Use variables from the configuration file to
                     output,     initialise the class private variables.
                     database)


        Private methods:
            _foldersize(folder) ; Return the size of a folder.
            _cleancache()       ; Keep the cache at a maximum size.
            _nametofile(name)   ; Convert a name to a filename.
            _write(raw_data,    ; Write a record to a file.
                    filename,
                    extract)
            _calcHash(content)  ; Calculate the md5sum of 'content'.
            _newUD()            ; Generate a new UD number.

        Public methods:
            retrieveslice(accno,   ; Retrieve a chromosome slice from the NCBI.
                          start,
                          stop,
                          orientation)
            retrievegene(gene,     ; Retrieve a gene from the NCBI.
                         organism,
                         upstream,
                         downstream)
            downloadrecord(url)    ; Download a GenBank file.
            uploadrecord(raw_data) ; Let someone upload a GenBank file.
            loadrecord(identifier) ; Load a record, store it in the
                                     cache, manage the cache and return
                                     the record.

        Inherited methods from Db.Output:
            WarningMsg(filename, message) ; Print a warning message.
            ErrorMsg(filename, message)   ; Print an error message and log it.
            LogMsg(filename, message)     ; Log a message.
    """

    def __init__(self, config, output, database) :
        """
            Use variables from  the configuration file for some simple
            settings. Make the cache directory if it does not exist yet.

            Inherited variables from Db.Output.Config:
                email     ; The email address which we give to the NCBI.
                cache     ; The directory where the records are stored.
        """

        self._config = config
        self._output = output
        self._database = database
        if not os.path.isdir(self._config.cache) :
            os.mkdir(self._config.cache)
        Entrez.email = self._config.email
    #__init__

    def _foldersize(self, folder) :
        """
            Return the size of a folder in bytes.

            Arguments:
                folder ; Name of a directory.

            Returns:
                integer ; The size of the directory.
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
                cache     ; Directory under scrutiny.
                cachesize ; Maximum size of the cache.
        """

        if self._foldersize(self._config.cache) < self._config.cachesize:
            return

        # Build a list of files sorted by access time.
        cachelist = []
        for (path, dirs, files) in os.walk(self._config.cache) :
            for filename in files :
                cachelist.append(
                    (os.stat(os.path.join(path, filename)).st_atime, filename))
        cachelist.sort()

        # Now start removing pairs of files until the size of the folder is
        # small enough (or until the list is exhausted).
        for i in range(0, len(cachelist)) :
            os.remove(os.path.join(path, cachelist[i][1]))
            if self._foldersize(self._config.cache) < self._config.cachesize:
                break;
        #for
    #_cleancache

    def _nametofile(self, name) :
        """
            Convert an accession number to a filename.

            Arguments:
                name ; The accession number.

            Inherited variables from Db.Output.Config:
                cache     ; Name of the cache directory.

            Returns:
                string ; A filename.
        """

        return self._config.cache + '/' + name + ".gb.bz2"
    #_nametofile

    def _write(self, raw_data, filename) :
        """
            Write raw data to a compressed file.

            Arguments:
                raw_data    ; The raw_data to be compressed and written
                filename    ; The intended name of the outfile

            Returns:
                outfile     ; The full paht and name of the file written
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

            Arguments:
                content ; Arbitrary text.

            Returns:
                string ; The md5sum of 'content'.
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

            Returns:
                string ; A new UD number.
        """

        M = Misc.Misc()
        UD = M.ID()
        del M
        return "UD_" + str(UD)
    #_newUD

    def _updateDBmd5(self, raw_data, name, GI):
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


    def snpConvert(self, rsId) :
        x = Entrez.efetch(db = "SNP", id = rsId, rettype = "flt", 
            retmode = "xml")
        
        doc = xml.dom.minidom.parseString(x.read())
        for i in doc.getElementsByTagName("hgvs") :
            self._output.addOutput("snp", i.lastChild.data.encode("utf8"))
    #snpConvert
#Retriever

class GenBankRetriever(Retriever):
    """
        TODO: Update docstring
    """
    def __init__(self, config, output, database):
        # Recall init of parent
        Retriever.__init__(self, config, output, database)
        # Child specific init

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

            Arguments:
                raw_data ; The data.
                filename ; The intended name of the file.
                extract  ; Flag that indicates whether to extract the record ID
                           and GI number:
                           0 ; Do not extract, use 'filename'.
                           1 ; Extract.

            Returns:
                tuple ; Depending on the value of 'extract':
                        0 ; ('filename', None)
                        1 ; (id, GI)

        """
        fakehandle = StringIO.StringIO() # Unfortunately, BioPython needs a
        fakehandle.write(raw_data)       # file handle.
        fakehandle.seek(0)
        try :
            record = SeqIO.read(fakehandle, "genbank")
        except (ValueError, AttributeError):  # An error occured while parsing.
            self._output.addMessage(__file__, 2, "ENOPARSE",
                    "The file could not be parsed.")
            fakehandle.close()
            return None, None
        #except

        if type(record.seq) == UnknownSeq :
            fakehandle.close()
            self._output.addMessage(__file__, 4, "ENOSEQ",
                "This record contains no sequence. Chromosomal or contig " \
                "records should be uploaded with the GenBank uploader.")
            return None, None
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
        net_handle = Entrez.efetch(db = "nucleotide", id = name, rettype = "gb")
        raw_data = net_handle.read()
        net_handle.close()

        if raw_data == "\n" :       # Check if the file is empty or not.
            self._output.addMessage(__file__, 4, "ERETR",
                                     "Could not retrieve %s." % name)
            return None
        #if
        else :                      # Something is present in the file.
            name, GI = self.write(raw_data, name, 1)
            if name :               # Processing went okay.
                return self._updateDBmd5(raw_data, name, GI)
            else :                  # Parse error in the GenBank file.
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

            Arguments:
                accno       ; The accession number of the chromosome.
                start       ; Start position of the slice.
                stop        ; End position of the slice.
                orientation ; Orientatiion of the slice:
                              1 ; Forward.
                              2 ; Reverse complement.

            Inherited variables from Db.Output.Config:
                maxDldSize ; Maximum size of the slice.

            Returns:
                string ; An UD number.
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

            Arguments:
                gene       ; Name of the gene.
                organism   ; The organism in which we search.
                upstream   ; Number of upstream nucleotides for the slice.
                downstream ; Number of downstream nucleotides for the slice.

            Returns:
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

            Arguments:
                url ; Location of a GenBank record.

            Inherited variables from Db.Output.Config:
                maxDldSize ; Maximum size of the file.
                minDldSize ; Minimum size of the file.
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

            Arguments:
                raw_data ; A GenBank record.
        """

        md5sum = self._calcHash(raw_data)
        UD = self._database.getGBFromHash(md5sum)
        if not UD :
            UD = self._newUD()
            self._database.insertGB(UD, None, md5sum, None, 0, 0, 0, None)
        #if
        else :
            if os.path.isfile(self._nametofile(UD)) :
                return UD

        return self.write(raw_data, UD, 0) and UD
    #uploadrecord

    def loadrecord(self, identifier) :
        """
            Load a record and return it.
            If the filename associated with the accession number is not found
            in the cache, try to re-download it.

            Arguments:
                identifier ; An accession number.

            Returns:
                record ; A GenBank.Record record
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
                    else :                  # This used to be a downloaded seq.
                        filename = self.downloadrecord(url) and filename
                #if
                else :                      # This used to be a slice.
                    filename = self.retrieveslice(*Loc) and filename
            #if
            else :                          # Never seen this name before.
                filename = self.fetch(name)
            #else
        #if

        # TODO: If filename is None an error occured
        if filename is None:
            #Notify batch to skip all instance of identifier
            self._output.addOutput("BatchFlags", ("S1", identifier))
            return None

        # Now we have the file, so we can parse it.
        GenBankParser = GBparser.GBparser()
        return GenBankParser.createGBRecord(filename)

    #loadrecord
#GenBankRetriever

class LargeRetriever(Retriever):
    """
    LargeRetriever Docstring
    """
    def __init__(self, config, output, database):
        # Recall init of parent
        Retriever.__init__(self, config, output, database)
        # Child specific init

    def loadrecord(self, identifier):
        """
            Load and parse a LRG file based on the identifier

            returns a GenRecord.Record instance on succes, None on failure
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
        record = LRGparser.createLrgRecord(file_handle.read())
        file_handle.close()

        return record
    #loadrecord

    def fetch(self, name):
        """
            Fetch the LRG file from the ebi FTP domain

            first try to grab the file from the confirmed section, if this
            fails, get it from the pending section.

            returns the full path to the file or None in case of an error
        """
        prefix = "ftp://ftp.ebi.ac.uk/pub/databases/lrgex/"
        url =        prefix + "%s.xml"          % name
        pendingurl = prefix + "pending/%s.xml"  % name

        try:
            return self.downloadrecord(url, name)
        except urllib2.URLError, e: #Catch error: file not found
            pass

        try:                # Try to get the file from the pending section
            filename = self.downloadrecord(pendingurl, name)
            self._output.addMessage(__file__, 2, "WPEND",
                "Warning: LRG file %s is a pending entry." % name)
            return filename
        except urllib2.URLError, e:
            self._output.addMessage(__file__, 4, "ERETR",
                                 "Could not retrieve %s." % name)
            return None             #Explicit return in case of an Error

    def downloadrecord(self, url, name = None) :
        """
            Download an LRG record from an URL.

            Arguments:
                url ; Location of the LRG record.

            Inherited variables from Db.Output.Config:
                maxDldSize  ; Maximum size of the file.
                minDldSize  ; Minimum size of the file.

            Returns:
                filename    ; The full path to the file
                None        ; in case of failure
        """
        lrgID = name or os.path.splitext(os.path.split(url)[1])[0]
        #if not lrgID.startswith("LRG"):
        #    return None
        filename = self._nametofile(lrgID)

        handle = urllib2.urlopen(url)
        info = handle.info()
        if info["Content-Type"] == "application/xml" :

            if not info.has_key("Content-length") : # FIXME
                raise urllib2.URLError("hallo")

            length = int(info["Content-Length"])
            if self._config.minDldSize < length < self._config.maxDldSize:
                raw_data = handle.read()
                handle.close()

                #Do an md5 check
                md5sum = self._calcHash(raw_data)
                md5db = self._database.getHash(lrgID)
                if md5db is None:
                    self._database.insertLRG(lrgID, md5sum, url)
                elif md5sum != md5db:           #hash has changed for the LRG ID
                    self._output.addMessage(__file__, -1, "WHASH",
                        "Warning: Hash of %s changed from %s to %s." % (
                        lrgID, md5db, md5sum))
                    self._database.updateHash(lrgID, md5sum)
                else:                           #hash the same as in db
                    pass

                if not os.path.isfile(filename) :
                    return self.write(raw_data, lrgID)
                else:
                    # This can only occur if synchronus calls to mutalyzer are 
                    # made to recover a file that did not exist. Still leaves a
                    # window in between the check and the write.
                    return filename
            #if
            else :
                self._output.addMessage(__file__, 4, "EFILESIZE",
                    "Filesize is not within the allowed boundaries.")
        #if
        else :
            self._output.addMessage(__file__, 4, "ERECPARSE",
                                     "This is not a LRG record.")
        handle.close()
    #downloadrecord

    def write(self, raw_data, filename) :
        """
            Write raw LRG data to a file. The data is parsed before writing,
            if a parse error occurs an error is returned and the function exits.

            Arguments:
                raw_data ; The data.
                filename ; The intended name of the file.

            Returns:
                filename ; The full path and name of the file written
                None     ; In case of an error

        """
        # Dirty way to test if a file is valid,
        # Parse the file to see if it's a real LRG file.
        try:
            LRGparser.createLrgRecord(raw_data)
        except DOMException:
            self._output.addMessage(__file__, 4, "ERECPARSE",
                                      "Could not parse file.")
            return None             # Explicit return on Error

        return self._write(raw_data, filename) #returns full path
    #write



if __name__ == "__main__" :
    pass
#if
