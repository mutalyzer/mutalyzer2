#!/usr/bin/python

import os
import sys
from mod_python import apache
from cStringIO import StringIO

#directory = os.path.dirname(__file__)
#VI = apache.import_module('Variant_info', path=[directory + "../"])
import Variant_info as VI

def __run(func, *args) :
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    func(*args)
    reply = sys.stdout.getvalue()
    sys.stdout = old_stdout

    return reply
#__run

def getTranscripts(v1, v2) :
    """
    #T = apache.import_module('getTranscripts', path=[directory + "/src"])
    from Services import getTranscripts as T

    return __run(T.main, v1, v2)
    """
    from Modules import Config
    from Modules import Db

    C = Config.Config()
    D = Db.Db(C)

    return str(D.get_Transcripts(v1, v2))
#getTranscripts

def getGeneName(v1) :
    """
    #G = apache.import_module('getGeneName', path=[directory + "/src"])
    from Services import getGeneName as G

    return __run(G.main, v1)
    """
    from Modules import Config
    from Modules import Db

    C = Config.Config()
    D = Db.Db(C)

    return str(D.get_GeneName(v1.split('.')[0]))
#getGeneName

def Variant_info(v1) :
    return __run(VI.main, v1)
