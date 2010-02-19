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
