#!/usr/bin/python

def getTranscripts(v1, v2) :
    from Modules import Config
    from Modules import Db
    from Modules import Output

    C = Config.Config()
    D = Db.Db(C)
    L = Output.Output(C, __file__)

    L.LogMsg(__file__, "Reveived request getTranscripts(%s %s)" % (v1, v2))
    ret = str(D.get_Transcripts(v1, v2))
    L.LogMsg(__file__, "Finished processing getTranscripts(%s %s)" % (v1, v2))

    return ret
#getTranscripts

def getGeneName(v1) :
    from Modules import Config
    from Modules import Db
    from Modules import Output

    C = Config.Config()
    D = Db.Db(C)
    L = Output.Output(C, __file__)

    L.LogMsg(__file__, "Reveived request getGeneName(%s)" % v1)
    ret = str(D.get_GeneName(v1.split('.')[0]))
    L.LogMsg(__file__, "Finished processing getGeneName(%s)" % v1)

    return ret
#getGeneName
