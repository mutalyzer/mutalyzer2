from ZSI import dispatch
from mod_python import apache, publisher
from simpletal import simpleTAL, simpleTALES
from cStringIO import StringIO

import sys
import os

myPath = os.path.dirname(__file__) + "/.."
if myPath not in sys.path :
    sys.path.append(myPath)
os.chdir(myPath + "/..")

from Interfaces import webservice

#mod = __import__('encodings.utf_8', globals(), locals(), '*')
#mod = __import__('encodings.utf_16_be', globals(), locals(), '*')

Mut_ver = "2.0&alpha;"

def __read(path, req) :
    handle = open(path + req.uri.split('/')[-1], "r")
    s = handle.read()
    handle.close

    return s
#__read

def __tal(filename, args) :
    context = simpleTALES.Context()

    for i in args :
        context.addGlobal(i, args[i])

    templateFile = open(filename, 'r')
    template = simpleTAL.compileXMLTemplate(templateFile)
    templateFile.close()

    string = StringIO()
    template.expand(context, string)

    return string.getvalue()
#__tal


def handler(req):
    if "services" in req.uri :
        dispatch.AsHandler(modules=(webservice,), request=req, rpc=True)
        return apache.OK
    #if
    if ".js" in req.uri :
        req.content_type = 'text/html'
        req.write(__read("templates/", req))
        return apache.OK
    #if
    if ".wsdl" in req.uri :
        req.content_type = 'text/xml'
        args = {
            "path"    : "localhost/mutalyzer2/services",
            "version" : Mut_ver
        }
        req.write(__tal("templates/service.wsdl", args))
        return apache.OK
    #if
    req.filename = myPath + "/Interfaces/" + req.filename.split('/')[-1]
    return publisher.handler(req)
#handler
