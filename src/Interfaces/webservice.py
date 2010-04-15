#!/usr/bin/python

from soaplib.wsgi_soap import SimpleWSGISoapApp
from soaplib.service import soapmethod
from soaplib.serializers.primitive import String, Integer, Array
from ZSI.fault import Fault

from Modules import Output

class MutalyzerService(SimpleWSGISoapApp) :
    """
        Mutalyzer webservices.
    
        These methods are made public via a SOAP interface.
    
        Private methods:
            __checkBuild(L, D, build) ; Check if the build is supported.
            __checkChrom(L, D, chrom) ; Check if the chromosome is in our 
                                        database.
            __checkPos(L, pos)        ; Check if the posision is valid.

        Public methods:
            getTranscripts(build, chrom, pos)    ; Get all transcripts that 
                                                   overlap with a chromosomal 
                                                   position.
            getTranscriptsRange(build, chrom,    ; Get all transcripts that 
                                pos1, pos2,        overlap with a range on a 
                                method)            chromosome.                 
            getGeneName(build, accno)            ; Find the gene name associated
                                                   with a transcript.
            varInfo(LOVD_ver, build, accno, var) ; Convert g. to c. and vice 
                                                   versa.
    """
    
    def __checkBuild(self, L, D, build) :
        """
            Check if the build is supported (hg18 or hg19).

            Arguments:
                L     ; An output object for logging.
                D     ; A handle to the database.
                build ; The build name that needs to be checked.

            Returns:
                Nothing (but raises an EARG exception).
        """

        if not D.opened :
            L.LogMsg(__file__, "EARG %s" % build)
            raise Fault(Fault.Client, "EARG", 
                detail = "The build argument (%s) was not a valid " \
                         "build name." % build)
        #if                         
    #__checkBuild

    def __checkChrom(self, L, D, chrom) :
        """
            Check if the chromosome is in our database.

            Arguments:
                L     ; An output object for logging.
                D     ; A handle to the database.
                chrom ; The name of the chromosome.

            Returns:
                Nothing (but raises an EARG exception).
        """

        if not D.isChrom(chrom) :
            L.LogMsg(__file__, "EARG %s" % chrom)
            raise Fault(Fault.Client, "EARG", 
                detail = "The chrom argument (%s) was not a valid " \
                         "chromosome name." % chrom)
        #if                         
    #__checkChrom

    def __checkPos(self, L, pos) :
        """
            Check if the posision is valid.

            Arguments:
                L   ; An output object for logging.
                pos ; The position.

            Returns:
                Nothing (but raises an ERANGE exception).
        """

        if pos < 1 :
            L.LogMsg(__file__, "ERANGE %i" % pos)
            raise Fault(Fault.Client, "ERANGE", 
                detail = "The pos argument (%i) is out of range." % pos)
        #if                         
    #__checkPos

    @soapmethod(String, String, Integer, _returns = Array(String))
    def getTranscripts(self, build, chrom, pos) :
        """
            Get all the transcripts that overlap with a chromosomal position.
    
            Arguments:
                string build ; The build name encoded as "hg18" or "hg19".
                string chrom ; A chromosome encoded as "chr1", ..., "chrY".
                int    pos   ; A postion on the chromosome.
    
            Returns:
                string ; A list of transcripts.

            On error an exception is raised:
                detail       ; Human readable description of the error.
                faultstring: ; A code to indicate the type of error.
                    EARG   ; The argument was not valid.
                    ERANGE ; An invalid range was given.
        """
    
        from Modules import Db
        from Modules import Output
    
        L = Output.Output(__file__)
    
        L.LogMsg(__file__, "Received request getTranscripts(%s %s %s)" % (
            build, chrom, pos))

        D = Db.Db("local", build)
        self.__checkBuild(L, D, build)

        self.__checkChrom(L, D, chrom)
        self.__checkPos(L, pos)
        
        ret = str(D.get_Transcripts(chrom, pos, pos, True))
        L.LogMsg(__file__, "Finished processing getTranscripts(%s %s %s)" % (
                 build, chrom, pos))
    
        del L
        del D
        return ret
    #getTranscripts
    
    @soapmethod(String, String, Integer, Integer, Integer, 
                _returns = Array(String))
    def getTranscriptsRange(self, build, chrom, pos1, pos2, method) :
        """
            Get all the transcripts that overlap with a range on a chromosome.
    
            Arguments:
                string build ; The build name encoded as "hg18" or "hg19".
                string chrom  ; A chromosome encoded as "chr1", ..., "chrY".
                int    pos1   ; The first postion of the range.
                int    pos2   ; The last postion of the range.
                int    method ; The method of determining overlap:
                                0 ; Return all hit transcripts.
                                1 ; Return only the transcripts that completely
                                    fall in the range [pos1, pos2].
    
            Returns:
                string ; A list of transcripts.
        """
    
        from Modules import Db
        from Modules import Output
    
        L = Output.Output(__file__)
    
        L.LogMsg(__file__, 
            "Received request getTranscriptsRange(%s %s %s %s %s)" % (
            build, chrom, pos1, pos2, method))

        D = Db.Db("local", build)
        self.__checkBuild(L, D, build)

        ret = str(D.get_Transcripts(chrom, pos1, pos2, method))
        L.LogMsg(__file__, 
            "Finished processing getTranscriptsRange(%s %s %s %s %s)" % (
            build, chrom, pos1, pos2, method))
    
        del D
        del L
        return ret
    #getTranscriptsRange

    @soapmethod(String, String, _returns = String)
    def getGeneName(self, build, accno) :
        """
            Find the gene name associated with a transcript.
    
            Arguments:
                string build ; The build name encoded as "hg18" or "hg19".
                string accno ; The identifier of a transcript.
    
            Returns:
                string ; The name of the associated gene.
        """
    
        from Modules import Db
        from Modules import Output
    
        
        L = Output.Output(__file__)
        L.LogMsg(__file__, "Received request getGeneName(%s %s)" % (build, 
            accno))

        D = Db.Db("local", build)
        self.__checkBuild(L, D, build)
    
        ret = str(D.get_GeneName(accno.split('.')[0]))
        L.LogMsg(__file__, "Finished processing getGeneName(%s %s)" % (build,
            accno))
    
        del L
        del D
        return ret
    #getGeneName

    @soapmethod(String, String, String, String, _returns = String)
    def varInfo(self, LOVD_ver, build, accno, var) :
        """
            Search for an NM number in the MySQL database, if the version
            number matches, get the start and end positions in a variant and
            translate these positions to g. notation if the variant is in c.
            notation and vice versa.

            - If no end position is present, the start position is assumed to
              be the end position. 
            - If the version number is not found in the database, an error
              message is generated and a suggestion for an other version is
              given.
            - If the reference sequence is not found at all, an error is
              returned.
            - If no variant is present, the transcription start and end and CDS
              end in c. notation is returned.
            - If the variant is not accepted by the nomenclature parser, a
              parse error will be printed.

            
            Arguments:
                LOVD_ver ; The LOVD version.
                build ; The human genome build (ignored for now, hg19 assumed).
                accno ; The NM accession number and version.
                var ; The variant, or empty.
             
            Returns:
                string:
                    start_main   ; The main coordinate of the start position 
                                   in c. (non-star) notation.
                    start_offset ; The offset coordinate of the start position
                                   in c. notation (intronic position).
                    end_main     ; The main coordinate of the end position in 
                                   c. (non-star) notation.
                    end_offset   ; The offset coordinate of the end position in
                                   c. notation (intronic position).
                    start_g      ; The g. notation of the start position.
                    end_g        ; The g. notation of the end position.
                    type         ; The mutation type.

            Returns (alternative):
                string:
                    trans_start  ; Transcription start in c. notation.
                    trans_stop   ; Transcription stop in c. notation.
                    CDS_stop     ; CDS stop in c. notation.
        """

        import Variant_info as VI
        from Modules import Web
        from Modules import Db
        from Modules import Output
    
        L = Output.Output(__file__)

        L.LogMsg(__file__, "Received request varInfo(%s %s %s %s)" % (
                 LOVD_ver, build, accno, var))

        D = Db.Db("local", build)
        self.__checkBuild(L, D, build)

        W = Web.Web()
        result = W.run(VI.main, LOVD_ver, build, accno, var)
        del W

        L.LogMsg(__file__, "Finished processing varInfo(%s %s %s %s)" % (
                 LOVD_ver, build, accno, var))

        del L
        del D
        return str(result.split("\n")[:-1])
    #varInfo
#MutalyzerService
