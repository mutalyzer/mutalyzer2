#!/usr/bin/python

import os
from ZSI import dispatch
from mod_python import apache, publisher

from Modules import Web
from Interfaces import webservice

def handler(req):
    myPath = os.path.dirname(__file__) + "/.."
    os.chdir(myPath + "/..")

    W = Web.Web()

    if "services" in req.uri :
        dispatch.AsHandler(modules=(webservice,), request=req, rpc=True)
        return apache.OK
    #if
    if ".js" in req.uri :
        req.content_type = 'text/html'
        req.write(W.read("templates/", req))
        return apache.OK
    #if
    if ".wsdl" in req.uri :
        req.content_type = 'text/xml'
        args = {
            "path"    : req.hostname + req.uri.rsplit('/', 1)[0] + "/services",
            "version" : W.version
        }
        req.write(W.tal("XML", "templates/service.wsdl", args))
        return apache.OK
    #if
    req.filename = myPath + "/Interfaces/" + req.filename.split('/')[-1]
    return publisher.handler(req)
#handler
