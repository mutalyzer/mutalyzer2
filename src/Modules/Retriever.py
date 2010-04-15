#!/usr/bin/python

from Output import Output

class Retriever(Output) :
    """
        Retrieve a record from either the cache or the NCBI.

        Inherited variables from Output.Config:
            email     ; The email address which we give to the NCBI.
            cache     ; The directory where the records are stored.
            cachesize ; Maximum size of the cache.

        Inherited functions from Output:
            output    ; The output object.

        Special methods:
            __init__(config) ; Use variables from the configuration file to 
                               initialise the class private variables.

        Private methods:
            __foldersize(folder) ; Return the size of a folder.
            __cleancache()       ; Keep the cache at a maximum size.

        Public methods:
            loadrecord(identifier) ; Load a record, store it in the 
                                     cache, manage the cache and return 
                                     the record.
    """

    def __init__(self) :
        """
            Use variables from  the configuration file for some simple
            settings. Make the cache directory if it does not exist yet.

            Arguments:
                config ; The configuration object.
                output ; The output object (for logging and error messages).

            Inherited variables from Config:
                email     ; The email address which we give to the NCBI.
                cache     ; The directory where the records are stored.
                cachesize ; Maximum size of the cache.
        """

        import os # os.path.isdir(), os.path.mkdir()

        Output.__init__(self, __file__)

        if not os.path.isdir(self.cache) :
            os.mkdir(self.cache)
    #__init__
    
    def __foldersize(self, folder) :
        """
            Return the size of a folder in megabytes.

            Arguments:
                folder ; Name of a directory.

            Returns: 
                integer ; The size of the directory.
        """

        import os # walk(), path.getsize(), path.join()

        folder_size = 0
        for (path, dirs, files) in os.walk(folder) :
            for file in files :
                folder_size += os.path.getsize(os.path.join(path, file))

        return folder_size / 1048576.0 
    #__foldersize
    
    def __cleancache(self) :
        """
            Keep removing files until the size of the cache is less than the
            given size.
            First, the cache checked for its size, if it exceeds the maximum
            size the ``oldest'' files are deleted. Note that accessing a file
            makes it ``new''.
            
            Inherited variables from Config: 
                cache     ; Directory under scrutiny.
                cachesize ; Maximum size of the cache.
        """

        import os # walk(), stat(), path.join(), remove()

        if self.__foldersize(self.cache) < self.cachesize :
            return
    
        # Build a list of files sorted by access time.
        cachelist = []
        for (path, dirs, files) in os.walk(self.cache) :
            for filename in files :
                cachelist.append(
                    (os.stat(os.path.join(path, filename)).st_atime, filename))
        cachelist.sort()
    
        # Now start removing pairs of files until the size of the folder is
        # small enough (or until the list is exhausted).
        for i in range(0, len(cachelist)) :
            os.remove(os.path.join(path, cachelist[i][1]))
            if self.foldersize(self.cache) < self.cachesize :
                break;
        #for
    #__cleancache
    
    def loadrecord(self, identifier) :
        """
            Return a record from either the cache or the NCBI. 
            If a file is retrieved from the NCBI, a hard link is made to its
            alternative name (GI when an accession number is given and vice
            versa). If no version is given, it will be retrieved and the
            record will be renamed.
            After downloading a file, the cache is checked for overflows by
            calling the __cleancache() function.
            The files are stored in compressed format in the cache.

            Variables: 
                identifier ; Either an accession number or a GI number.

            Inherited variables from Config:
                cache      ; The directory where the record is stored.
                email      ; The email address which we give to the NCBI.
                output     ; The output object.

            Returns:
                SeqRecord ; The record that was requested.
        """

        import os              # path.isfile(), link()
        import bz2             # BZ2Compressor(), BZ2File()
        from Bio import SeqIO  # read()
        from Bio import Entrez # efetch()

        # If a GI is given, remove the "GI" or "GI:" part.
        if (identifier[:2] == "GI") :
            if (identifier[2] == ':') :
                name = identifier[3:]
            else :
                name = identifier[2:]
        #if
        else :
            name = identifier
    
        # Make a filename based upon the identifier.
        filename = self.cache + '/' + name + ".gb.bz2"
    
        # If the filename is not present, retrieve it from the NCBI.
        if not os.path.isfile(filename) :
            Entrez.email = self.email
            net_handle = \
                Entrez.efetch(db = "nucleotide", id = name, rettype = "gb")
            raw_data = net_handle.read()
            net_handle.close()

            # Check if the record is empty or not.
            if raw_data != "\n" :
                # Compress it to save disk space.
                comp = bz2.BZ2Compressor()
                data = comp.compress(raw_data)
                data += comp.flush()
                out_handle = open(filename, "w")
                out_handle.write(data)
                out_handle.close()
    
                # Since we put something in the cache, check if it needs 
                # cleaning.
                self.__cleancache()
            #if
            else :
                self.ErrorMsg(__file__, "Could not retrieve %s." % name)
                return None
        #if
        
        # Now we have the file, so we can parse it.
        file_handle = bz2.BZ2File(filename, "r")
        try :
            record = SeqIO.read(file_handle, "genbank")
        except ValueError :
            self.ErrorMsg(__file__, "Could not parse %s, purging." % filename)
            os.remove(filename)
            file_handle.close()
            return None
        #except

        file_handle.close()

        # If a GI is supplied, find out the accession number (plus version)
        #   and vice versa.
        if name != record.annotations["gi"] :
            altfilename = \
                self.cache + '/' + record.annotations["gi"] + ".gb.bz2"
            altfilename2 = \
                self.cache + '/' + record.id + ".gb.bz2"
        #if
        else :
            altfilename = self.cache + '/' + record.id + ".gb.bz2"
    
        # If the alternative filename is not present yet, make a hard link.
        if not os.path.isfile(altfilename) :
            os.link(filename, altfilename)

        # If the other alternative filename is not present (no version was 
        #   given), rename the file. If it already exists, remove the file.
        if filename != altfilename2 :
            if not os.path.isfile(altfilename2) :
                os.rename(filename, altfilename2)
            else :
                os.remove(filename)

            self.WarningMsg(__file__, "No version number is given, " \
                  "using %s. Please use version numbers to reduce " \
                  "downloading overhead." % record.id)
        #if

        return record
    #loadrecord
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
