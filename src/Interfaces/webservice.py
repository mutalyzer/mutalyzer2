#!/usr/bin/python

from soaplib.wsgi_soap import SimpleWSGISoapApp
from soaplib.service import soapmethod
from soaplib.serializers.primitive import String, Integer, Array

class MutalyzerService(SimpleWSGISoapApp) :
    """
        Mutalyzer webservices.
    
        These methods are made public via a SOAP interface.
    
        Public methods:
            getTranscripts(v1, v2)              ; Get all transcripts that 
                                                  overlap with a chromosomal 
                                                  position.
            getTranscriptsRange(v1, v2, v3, v4) ; Get all transcripts that 
                                                  overlap with a range on a 
                                                  chromosome.
            getGeneName(v1)                     ; Find the gene name associated 
                                                  with a transcript.
            varInfo(v1, v2, v3, v4)             ; Convert g. to c. and vice 
                                                  versa.
    """
    
    @soapmethod(String, Integer, _returns = Array(String))
    def getTranscripts(self, v1, v2) :
        """
            Get all the transcripts that overlap with a chromosomal position.
    
            Arguments:
                string v1 ; A chromosome encoded as "chr1", ..., "chrY".
                int    v2 ; A postion on the chromosome.
    
            Returns:
                string ; A list of transcripts.
        """
    
        from Modules import Config
        from Modules import Db
        from Modules import Output
    
        C = Config.Config()
        D = Db.Db(C, "local")
        L = Output.Output(C, __file__)
    
        L.LogMsg(__file__, "Reveived request getTranscripts(%s %s)" % (v1, v2))
        ret = str(D.get_Transcripts(v1, v2, v2, True))
        L.LogMsg(__file__, "Finished processing getTranscripts(%s %s)" % (
                 v1, v2))
    
        del L
        del D
        del C
        return ret
    #getTranscripts
    
    @soapmethod(String, Integer, Integer, Integer, _returns = Array(String))
    def getTranscriptsRange(self, v1, v2, v3, v4) :
        """
            Get all the transcripts that overlap with a range on a chromosome.
    
            Arguments:
                string v1 ; A chromosome encoded as "chr1", ..., "chrY".
                int    v2 ; The first postion of the range.
                int    v3 ; The last postion of the range.
                int    v4 ; The method of determining overlap:
                            0 ; Return all hit transcripts.
                            1 ; Return only the transcripts that completely fall
                                in the range [v2, v3].
    
            Returns:
                string ; A list of transcripts.
        """
    
        from Modules import Config
        from Modules import Db
        from Modules import Output
    
        C = Config.Config()
        D = Db.Db(C, "local")
        L = Output.Output(C, __file__)
    
        L.LogMsg(__file__, 
            "Reveived request getTranscriptsRange(%s %s %s %s)" % (
            v1, v2, v3, v4))
        ret = str(D.get_Transcripts(v1, v2, v3, v4))
        L.LogMsg(__file__, 
            "Finished processing getTranscriptsRange(%s %s %s %s)" % (
            v1, v2, v3, v4))
    
        del L
        del D
        del C
        return ret
    #getTranscriptsRange
    
    @soapmethod(String, _returns = String)
    def getGeneName(self, v1) :
        """
            Find the gene name associated with a transcript.
    
            Arguments:
                string v1 ; The identifier of a transcript.
    
            Returns:
                string ; The name of the associated gene.
        """
    
        from Modules import Config
        from Modules import Db
        from Modules import Output
    
        C = Config.Config()
        D = Db.Db(C, "local")
        L = Output.Output(C, __file__)
    
        L.LogMsg(__file__, "Reveived request getGeneName(%s)" % v1)
        ret = str(D.get_GeneName(v1.split('.')[0]))
        L.LogMsg(__file__, "Finished processing getGeneName(%s)" % v1)
    
        del L
        del D
        del C
        return ret
    #getGeneName

    @soapmethod(String, String, String, String, _returns = String)
    def varInfo(self, v1, v2, v3, v4) :
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
                v1 ; The LOVD version.
                v2 ; The human genome build (ignored for now, hg19 assumed).
                v3 ; The NM accession number and version.
                v4 ; The variant, or empty.
             
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
        from Modules import Config
        from Modules import Db
        from Modules import Output
    
        C = Config.Config()
        D = Db.Db(C, "local")
        L = Output.Output(C, __file__)
    
        L.LogMsg(__file__, "Reveived request varInfo(%s %s %s %s)" % (
                 v1, v2, v3, v4))

        W = Web.Web()
        result = W.run(VI.main, v1, v2, v3, v4)
        del W

        L.LogMsg(__file__, "Finished processing varInfo(%s %s %s %s)" % (
                 v1, v2, v3, v4))

        del L
        del D
        del C
        return str(result.split("\n")[:-1])
    #varInfo
#MutalyzerService
