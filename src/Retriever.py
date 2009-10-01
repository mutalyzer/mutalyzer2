#!/usr/bin/python

class Retriever() :
    """
        Retrieve a record from either the cache or the NCBI.

        Private variables:
            __email     ; The email address which we give to the NCBI.
            __cache     ; The directory where the records are stored.
            __cachesize ; Maximum size of the cache.

        Special methods:
            __init__()             ; Read the configuration file and 
                                     initialise the class global 
                                     variables.

        Private methods:
            __foldersize(folder)     ; Return the size of a folder.
            __cleancache()           ; Keep the cache at a maximum size.

       Public methods:
            loadrecord(identifier) ; Load a record, store it in the 
                                     cache, manage the cache and return 
                                     the record.
    """

    def __init__(self) :
        """
            Read the configuration file for some simple settings.

            Private variables (altered):
                __email     ; The email address which we give to the NCBI.
                __cache     ; The directory where the records are stored.
                __cachesize ; Maximum size of the cache.
        """
        from configobj import ConfigObj # ConfigObj()

        config = ConfigObj("./mutalyzer.conf")
        self.__email = config["email"]
        self.__cache = config["cache"]
        self.__cachesize = int(config["cachesize"])
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
            
            Private variables: 
                __cache     ; Directory under scrutiny.
                __cachesize ; Maximum size of the cache.
        """
        import os # walk(), stat(), path.join(), remove()

        if self.__foldersize(self.__cache) < self.__cachesize :
            return
    
        # Build a list of files sorted by access time.
        list = []
        for (path, dirs, files) in os.walk(self.__cache) :
            for file in files :
                list.append((os.stat(os.path.join(path, file)).st_atime, file))
        list.sort()
    
        # Now start removing pairs of files until the size of the folder is
        # small enough (or until the list is exhausted).
        for i in range(0, len(list)) :
            os.remove(os.path.join(path, list[i][1]))
            if self.__foldersize(self.__cache) < self.__cachesize :
                break;
        #for
    #__cleancache
    
    def loadrecord(self, identifier) :
        """
            Return a record from either the cache or the NCBI. 
            If a file is retrieved from the NCBI, a hard link is made to its
            alternative name (GI when an accession number is given and vice
            versa).
            After downloading a file, the cache is checked for overflows by
            calling the __cleancache() function.
            The files are stored in compressed format in the cache.

            Variables: 
                identifier ; Either an accession number or a GI number.

            Private variables:
                __cache      ; The directory where the record is stored.
                __email      ; The email address which we give to the NCBI.

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
        filename = self.__cache + '/' + name + ".gb.bz2"
    
        # If the filename is not present, retrieve it from the NCBI.
        if not os.path.isfile(filename) :
            Entrez.email = self.__email
            net_handle = \
                Entrez.efetch(db = "nucleotide", id = name, rettype = "gb")
            # Compress it to save disk space.
            comp = bz2.BZ2Compressor()
            data = comp.compress(net_handle.read())
            data += comp.flush()
            out_handle = open(filename, "w")
            out_handle.write(data)
            out_handle.close()
            net_handle.close()
    
            # Since we put something in the cache, check if it needs cleaning.
            self.__cleancache()
        #if
        
        # Now we have the file, so we can parse it.
        file_handle = bz2.BZ2File(filename, "r")
        record = SeqIO.read(file_handle, "genbank")
        file_handle.close()

        # If a GI is supplied, find out the accession number and vice versa.
        if name == record.annotations["gi"] :
            altfilename = self.__cache + '/' + record.id + ".gb.bz2"
        else :
            altfilename = \
                self.__cache + '/' + record.annotations["gi"] + ".gb.bz2"
    
        # If the alternative filename is not present yet, make a hard link.
        if not os.path.isfile(altfilename) :
            os.link(filename, altfilename)

        return record
    #loadrecord
#Retriever
