#!/usr/bin/python

import os
import sys
from mod_python import apache
from cStringIO import StringIO
#from simpletal import simpleTAL, simpleTALES

directory = os.path.dirname(__file__)
VI = apache.import_module('Variant_info', path=[directory + "/src"])
os.chdir(directory)

def __run(func, *args) :
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    func(*args)
    reply = sys.stdout.getvalue()
    sys.stdout = old_stdout

    return reply
#__run

def getTranscripts(v1, v2) :
    T = apache.import_module('getTranscripts', path=[directory + "/src"])

    return __run(T.main, v1, v2)
#getTranscripts

def getGeneName(v1) :
    G = apache.import_module('getGeneName', path=[directory + "/src"])

    return __run(G.main, v1)
#getGeneName

def Variant_info(v1) :
    return __run(VI.main, v1)
