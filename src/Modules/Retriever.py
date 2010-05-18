#!/usr/bin/python

import os              # path.isfile(), link() path.isdir(), path.mkdir(),
                       # walk(), path.getsize(), path.join(), stat(), remove()
import bz2             # BZ2Compressor(), BZ2File()
import hashlib         # md5(), update(), hexdigest()
import urllib2         # urlopen()
import StringIO        # StringIO()
from Bio import SeqIO  # read()
from Bio import Entrez # efetch(), read(), esearch(), esummary()

import Misc

#from Output import Output
#from Db import Db

class Retriever() :
    """
        Retrieve a record from either the cache or the NCBI.

        Inherited variables from Db.Output.Config:
            email     ; The email address which we give to the NCBI.
            cache     ; The directory where the records are stored.
            cachesize ; Maximum size of the cache.

        Special methods:
            __init__(config) ; Use variables from the configuration file to 
                               initialise the class private variables.

        Private methods:
            __foldersize(folder) ; Return the size of a folder.
            __cleancache()       ; Keep the cache at a maximum size.
            __nametofile(name)   ; Convert a name to a filename.
            __write(raw_data,    ; Write a record to a file.
                    filename, 
                    extract) 
            __calcHash(content)  ; Calculate the md5sum of 'content'.
            __newUD()            ; Generate a new UD number.

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

        #Db.__init__(self, "local", "mutalyzer")

        self.__config = config
        self.__output = output
        self.__database = database
        if not os.path.isdir(self.__config.cache) :
            os.mkdir(self.__config.cache)
        Entrez.email = self.__config.email
    #__init__
    
    def __foldersize(self, folder) :
        """
            Return the size of a folder in bytes.

            Arguments:
                folder ; Name of a directory.

            Returns: 
                integer ; The size of the directory.
        """

        folder_size = 0
        for (path, dirs, files) in os.walk(folder) :
            for file in files :
                folder_size += os.path.getsize(os.path.join(path, file))

        return folder_size
    #__foldersize
    
    def __cleancache(self) :
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

        if self.__foldersize(self.__config.cache) < self.__config.cachesize:
            return
    
        # Build a list of files sorted by access time.
        cachelist = []
        for (path, dirs, files) in os.walk(self.__config.cache) :
            for filename in files :
                cachelist.append(
                    (os.stat(os.path.join(path, filename)).st_atime, filename))
        cachelist.sort()
    
        # Now start removing pairs of files until the size of the folder is
        # small enough (or until the list is exhausted).
        for i in range(0, len(cachelist)) :
            os.remove(os.path.join(path, cachelist[i][1]))
            if self.__foldersize(self.__config.cache) < self.__config.cachesize:
                break;
        #for
    #__cleancache

    def __nametofile(self, name) :
        """
            Convert an accession number to a filename.
            
            Arguments:
                name ; The accession number.

            Inherited variables from Db.Output.Config: 
                cache     ; Name of the cache directory.

            Returns:
                string ; A filename.
        """

        return self.__config.cache + '/' + name + ".gb.bz2"
    #__nametofile

    def __write(self, raw_data, filename, extract) :
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

        # Parse the file to see if it's a real GenBank file.
        fakehandle = StringIO.StringIO() # Unfortunately, BioPython needs a
        fakehandle.write(raw_data)       # file handle.
        fakehandle.seek(0)
        try :
            record = SeqIO.read(fakehandle, "genbank")
        except ValueError :              # An error occured while parsing.
            fakehandle.close()
            self.__output.addMessage(__file__, 4, "ERECPARSE", 
                                     "Could not parse file.")
            return None
        #except
        
        outfile = filename
        GI = None
        if extract :
            outfile = record.id
            GI = record.annotations["gi"]
            if outfile != filename :
                self.__output.addMessage(__file__, 2, "WNOVER", 
                    "No version number is given, using %s. Please use this " \
                    "number to reduce downloading overhead." % record.id)
        #if
        fakehandle.close()

        # Compress the data to save disk space.
        comp = bz2.BZ2Compressor()
        data = comp.compress(raw_data)
        data += comp.flush()
        out_handle = open(self.__nametofile(outfile), "w")
        out_handle.write(data)
        out_handle.close()
        
        # Since we put something in the cache, check if it needs cleaning.
        self.__cleancache()

        return outfile, GI
    #__write

    def __calcHash(self, content) :
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
    #__calcHash

    def __newUD(self) :
        """
            Make a new UD number based on the current time (seconds since 1970).

            Returns:
                string ; A new UD number.
        """

        M = Misc.Misc()
        UD = M.ID()
        del M
        return "UD_" + str(UD)
    #__newUD

    def __eFetch(self, name) :
        net_handle = Entrez.efetch(db = "nucleotide", id = name, rettype = "gb")
        raw_data = net_handle.read()
        net_handle.close()
        
        if raw_data == "\n" :       # Check if the file is empty or not.
            self.__output.addMessage(__file__, 4, "ERETR", 
                                     "Could not retrieve %s." % name)
            return None
        #if
        else :                      # Something is present in the file.
            name, GI = self.__write(raw_data, name, 1)
            if name :               # Processing went okay.
                currentmd5sum = self.__database.getHash(name)
                md5sum = self.__calcHash(raw_data)
                if md5sum != currentmd5sum :
                    self.__output.addMessage(__file__, -1, "WHASH", 
                        "Warning: Hash of %s changed from %s to %s." % (
                        name, currentmd5sum, md5sum))
                    self.__database.updateHash(name, md5sum)
                #if
                else :
                    self.__database.insertGB(name, GI, 
                        self.__calcHash(raw_data), None, 0, 0, 0, None)
                return self.__nametofile(name)
            #if
            else :                  # Parse error in the GenBank file.
                return None
         #else
    #__eFetch                        

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

        # The slice can not be too big.
        if stop - start > self.__config.maxDldSize :
            return None

        # Check whether we have seen this slice before.
        UD = self.__database.getGBFromLoc(accno, start, stop, orientation) 
        if UD : # This has been requested before.
            if os.path.isfile(self.__nametofile(UD)) : # It's still present.
                return UD

        # It's not present, so download it.
        handle = Entrez.efetch(db = "nucleotide", rettype = "gb", 
                               id = accno, seq_start = start, seq_stop = stop, 
                               strand = orientation)
        raw_data = handle.read()                               
        handle.close()                               

        # Calculate the hash of the downloaded file.
        md5sum = self.__calcHash(raw_data)

        if UD : # We have seen this one before.
            currentmd5sum = self.__database.getHash(UD)
            if md5sum != currentmd5sum :
                self.__output.addMessage(__file__, -1, "WHASH", 
                    "Warning: Hash of %s changed from %s to %s." % (
                    UD, currentmd5sum, md5sum))
                self.__database.updateHash(UD, md5sum)
            #if
        else : # We haven't seen it before, so give it a name.
            UD = self.__newUD()
            self.__database.insertGB(UD, None, md5sum, accno, start, 
                          stop, orientation, None)
        #else
                               
        self.__write(raw_data, UD, 0)

        return UD
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
        """

        # Search the NCBI for a specific gene in an organism.
        query = "%s[Gene] AND %s[Orgn]" % (gene, organism)
        handle = Entrez.esearch(db = "gene", term = query)
        searchresult = Entrez.read(handle)
        handle.close()

        # FIXME
        if len(searchresult["IdList"]) > 1 :
            print "Hmmmmm."
            return None
        #if

        # Get summary information for the first search hit.
        handle = Entrez.esummary(db = "gene", id = searchresult["IdList"][0])
        summary = Entrez.read(handle)
        handle.close()
        if not len(summary[0]["GenomicInfo"]) :
            print "No mapping information found."
            #FIXME Output and stuff.
            return
        #if
        ChrAccVer = summary[0]["GenomicInfo"][0]["ChrAccVer"] # Extract the mapping
        ChrLoc = summary[0]["GenomicInfo"][0]["ChrLoc"]       # information.
        ChrStart = summary[0]["GenomicInfo"][0]["ChrStart"]
        ChrStop = summary[0]["GenomicInfo"][0]["ChrStop"]

        # Figure out the orientation of the gene.
        orientation = "1"
        if ChrStart > ChrStop :             # Swap start and stop.
            orientation = "2"
            temp = ChrStart
            ChrStart = ChrStop - downstream # Also take care of the flanking
            ChrStop = temp + upstream       # sequences.
        #if
        else :
            ChrStart -= upstream
            ChrStop += downstream
        #else

        # And retrieve the slice.
        self.retrieveslice(ChrAccVer, ChrStart, ChrStop, orientation)
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
            if length > self.__config.minDldSize and \
               length < self.__config.maxDldSize :
                raw_data = handle.read()
                md5sum = self.__calcHash(raw_data)
                UD = self.__database.getGBFromHash(md5sum)
                if not UD :
                    UD = self.__newUD()
                    self.__database.insertGB(UD, None, md5sum, None, 0, 0, 0, url)
                #if
                if not os.path.isfile(self.__nametofile(UD)) :
                    self.__write(raw_data, UD, 0)
            #if
            else :
                self.__output.addMessage(__file__, 4, "EFILESIZE", 
                    "Filesize is not within the allowed boundaries.")
        #if
        else :
            self.__output.addMessage(__file__, 4, "ERECPARSE", 
                                     "This is not a GenBank record.")
        handle.close()
    #downloadrecord

    def uploadrecord(self, raw_data) :
        """
            Write an uploaded record to a file.
            If the downloaded file is recognised by its hash, the old UD number
            is used.

            Arguments:
                raw_data ; A GenBank record.
        """

        md5sum = self.__calcHash(raw_data)
        UD = self.getGBFromHash(md5sum)
        if not UD :
            UD = self.__newUD()
            self.__database.insertGB(UD, None, md5sum, None, 0, 0, 0, None)
        #if
        else :
            if os.path.isfile(self.__nametofile(UD)) :
                return

        self.__write(raw_data, UD, 0)
    #uploadrecord
    
    def loadrecord(self, identifier) :
        """
            Load a record and return it.
            If the identifier is a GI number, try to find its record ID in the
            local database.
            If the filename associated with the accession number is not found
            in the cache, first try to re-download it by checking the following:
            - If slicing information can be found, make a new slice.
            - If an URL can be found, download it again.
            - If nothing can be found and the GI number can not be found, ask
              the user to upload the file again.
            In other cases, download the RefSeq sequence, extract the record ID
            and GI number and put them in the local database.

            Arguments:
                identifier ; An accession number.

            Returns:
                record ; A Genbank record.
        """

        if (identifier[0].isdigit()) : # This is a GI identifier.
            name = self.__database.getGBFromGI(identifier)
        else :
            name = identifier
        
        # Make a filename based upon the identifier.
        filename = self.__nametofile(name)

        if not os.path.isfile(filename) :   # We can't find the file.
            md5 = self.__database.getHash(name)
            if md5 :                        # We have seen it before though.
                Loc = self.__database.getLoc(name)  # Try to find the location.
                if not Loc[0]:              # No location found.
                    url = self.__database.getUrl(name)   # Try to find an URL.
                    if not url :
                        if self.__database.getGI(name) : # It was from NCBI.
                            filename = self.__eFetch(name)
                        else :
                            self.__output.addMessage(__file__, 4, "ERETR", 
                                "Please upload this sequence again.")
                            return None
                    #if
                    else :                  # This used to be a downloaded seq.
                        self.downloadrecord(url)
                #if
                else :                      # This used to be a slice.
                    self.retrieveslice(*Loc)
            #if
            else :                          # Never seen this name before.
                filename = self.__eFetch(name)
                """
                net_handle = Entrez.efetch(db = "nucleotide", id = name, 
                                           rettype = "gb")
                raw_data = net_handle.read()
                net_handle.close()
                
                if raw_data == "\n" :       # Check if the file is empty or not.
                    self.__output.addMessage(__file__, 4, "ERETR", 
                        "Could not retrieve %s." % name)
                    return None
                #if
                else :                      # Something is present in the file.
                    name, GI = self.__write(raw_data, name, 1)
                    if name :               # Processing went okay.
                        self.__database.insertGB(name, GI, 
                            self.__calcHash(raw_data), None, 0, 0, 0, None)
                        filename = self.__nametofile(name)
                    #if
                    else :                  # Parse error in the GenBank file.
                        return None
                #else
                """
            #else
        #if
                
        # Now we have the file, so we can parse it.
        file_handle = bz2.BZ2File(filename, "r")
        record = SeqIO.read(file_handle, "genbank")
        file_handle.close()

        return record
    #loadrecord

    #'''
    #def loadrecord_old(self, identifier) :
    #    """
    #        Return a record from either the cache or the NCBI. 
    #        If a file is retrieved from the NCBI, a hard link is made to its
    #        alternative name (GI when an accession number is given and vice
    #        versa). If no version is given, it will be retrieved and the
    #        record will be renamed.
    #        After downloading a file, the cache is checked for overflows by
    #        calling the __cleancache() function.
    #        The files are stored in compressed format in the cache.

    #        Variables: 
    #            identifier ; Either an accession number or a GI number.

    #        Inherited variables from Config:
    #            cache      ; The directory where the record is stored.
    #            email      ; The email address which we give to the NCBI.
    #            output     ; The output object.

    #        Returns:
    #            SeqRecord ; The record that was requested.
    #    """

    #    # If a GI is given, remove the "GI" or "GI:" part.
    #    if (identifier[:2] == "GI") :
    #        if (identifier[2] == ':') :
    #            name = identifier[3:]
    #        else :
    #            name = identifier[2:]
    #    #if
    #    else :
    #        name = identifier
    #
    #    # Make a filename based upon the identifier.
    #    filename = self.__nametofile(name)
    #
    #    # If the filename is not present, retrieve it from the NCBI.
    #    if not os.path.isfile(filename) :
    #        Loc = self.__database.getLoc(name) # Look in the UD database
    #        if not Loc :            # Never seen this name before.
    #            net_handle = \
    #                Entrez.efetch(db = "nucleotide", id = name, rettype = "gb")
    #            raw_data = net_handle.read()
    #            net_handle.close()
    #            # Check if the record is empty or not.
    #            if raw_data != "\n" :
    #                self.__write(raw_data, filename, 1)
    #            else :
    #                self.__output.addMessage(__file__, 4, "ERETR", 
    #                                         "Could not retrieve %s." % name)
    #                return None
    #        #if                    
    #        else :                  # We know the name.
    #            self.retrieveslice(*Loc)
    #    #if
    #    
    #    # Now we have the file, so we can parse it.
    #    file_handle = bz2.BZ2File(filename, "r")
    #    try :
    #        record = SeqIO.read(file_handle, "genbank")
    #    except ValueError :
    #        self.__output.addMessage(__file__, 4, "ERECPARSE", 
    #                                 "Could not parse %s, purging." % filename)
    #        os.remove(filename)
    #        file_handle.close()
    #        return None
    #    #except

    #    file_handle.close()

    #    if name[:3] == "UD_" : # No renaming is needed.
    #        return record

    #    # If a GI is supplied, find out the accession number (plus version)
    #    #   and vice versa.
    #    if name != record.annotations["gi"] :
    #        altfilename = self.__nametofile(record.annotations["gi"])
    #        altfilename2 = self.__nametofile(record.id)
    #    #if
    #    else :
    #        altfilename = self.__nametofile(record.id)
    #
    #    # If the alternative filename is not present yet, make a hard link.
    #    if not os.path.isfile(altfilename) :
    #        os.link(filename, altfilename)

    #    # If the other alternative filename is not present (no version was 
    #    #   given), rename the file. If it already exists, remove the file.
    #    if filename != altfilename2 :
    #        if not os.path.isfile(altfilename2) :
    #            os.rename(filename, altfilename2)
    #        else :
    #            os.remove(filename)

    #        self.__output.addMessage(__file__, 2, "WNOVRE", 
    #        "No version number is given, using %s. Please use version numbers to reduce " \
    #              "downloading overhead." % record.id)
    #    #if

    #    return record
    ##loadrecord
    #'''
#Retriever

#
# Unit test.
#
if __name__ == "__main__" :
    # Get the location of the cache, the cachesize and the email address from 
    #   the config file.
    R = Retriever()

    R.loadrecord("AB026906.1") # Retrieve a GenBank record.
    del R
#if
