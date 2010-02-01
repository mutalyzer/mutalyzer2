from ZSI import dispatch
from mod_python import apache, publisher
from simpletal import simpleTAL, simpleTALES
from cStringIO import StringIO

import webservice
#import index

#mod = __import__('encodings.utf_8', globals(), locals(), '*')
#mod = __import__('encodings.utf_16_be', globals(), locals(), '*')

Mut_ver = "2.0&alpha;"

def __read(path, req) :
    handle = open(req.uri.split('/')[-1], "r")
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
    if "html" in req.uri :
        req.content_type = 'text/html'
        req.write(__read("./html/", req))
        return apache.OK
    #if
    if "wsdl" in req.uri :
        req.content_type = 'text/xml'
        args = {
            "path"    : "localhost/mutalyzer2/services",
            "version" : Mut_ver
        }
        req.write(__tal("service.wsdl", args))
        return apache.OK
    #if
    return publisher.handler(req)
#handler
