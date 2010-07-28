#!/usr/bin/python

"""
    General handler for mod_python.

    This handler calls the webservices handler, the HTML publisher or handles
    a request itself, depending on keywords in the URI of the request.

    Public methods:
        handler(req) ; The handler.
"""

import os
import bz2
from ZSI import dispatch
from soaplib.client import make_service_client
from mod_python import apache, publisher

from Modules import Web
from Modules import Config
from Modules import File

import webservice

def handler(req):
    """
        Handle a request passed to us by mod_python.

        Keywords in the URI of the request are used to decide what to do:
            "services" ; Dispatch the webservices handler.
            ".js"      ; Return the raw content of the file (to include
                         JavaScript from an HTML file).
            ".py"      ; Return the content as a downloadable file after it has
                         been processed by TAL (to generate webservice client
                         files).
            ".wsdl"    ; Return the content of the file after it has been
                         processed by TAL (to generate a WSDL file that refers
                         to the correct server).
        By default, the HTML publisher is used for normal HTML files and other
        unhandled requests.

        Arguments:
            req ; The request.

        Returns:
            int ; An apache return code, either generated in this function
                  itself, or by the publisher handler which handles normal HTML
                  requests.
    """

    # Figure out where this program is located and go to the parent directory.
    myPath = os.path.dirname(__file__)
    os.chdir(myPath + "/..")

    reqPath = req.hostname + req.uri.rsplit('/', 1)[0]

    W = Web.Web()

    # Dispatch the webservices handler.
    if "services" in req.uri :
        Export = webservice.MutalyzerService()
        dispatch.AsHandler(modules=(Export,), request=req, rpc=True)
        return apache.OK
    #if

    # Return raw content (for includes from an HTML file).
    if req.uri.endswith(".js") or "base" in req.uri :
        #req.content_type = 'text/html'
        req.write(W.read("templates/", req))
        return apache.OK
    #if

    # Make all txt files in downloads directory downloadable
    if "downloads/" in req.uri :
        reqFile = req.uri.rsplit('/')[-1]
        req.headers_out["Content-Disposition"] = \
            "attachment; filename = \"%s\"" % reqFile # To force downloading.
        handle = open("templates/downloads/" + reqFile)
        C = Config.Config()
        F = File.File(C.File, None)
        req.content_type = F.getMimeType(handle)[0]
        req.write(handle.read())
        handle.close()
        return apache.OK
    #if

    # Process the file with TAL and return the content as a downloadable file.
    if req.uri.endswith(".py") :
        reqFile = req.uri.split('/')[-1]
        req.content_type = 'text/plain'
        req.headers_out["Content-Disposition"] = \
            "attachment; filename = \"%s\"" % reqFile # To force downloading.
        req.write(W.tal("TEXT", "templates/" + reqFile, {"path": reqPath}))
        return apache.OK
    #if

    # Return raw content (for batch checker results).
    if "Results" in req.uri:
        reqFile = req.uri.split('/')[-1]
        C = Config.Config()
        req.content_type = 'text/plain'
        req.headers_out["Content-Disposition"] = \
            "attachment; filename = \"%s\"" % reqFile # To force downloading.
        fh = open("%s/%s" % (C.Scheduler.resultsDir, reqFile))
        req.write(fh.read())
        fh.close()
        return apache.OK
    #if

    # Return uncompressed GenBank files from the cache.
    if "GenBank" in req.uri:
        reqFile = req.uri.split('/')[-1]
        C = Config.Config()
        fileName = "%s/%s.bz2" % (C.Retriever.cache, reqFile)
        if os.path.isfile(fileName) :
            handle = bz2.BZ2File("%s/%s.bz2" % (C.Retriever.cache, reqFile),
                                 "r")
            req.content_type = 'text/plain'
            req.headers_out["Content-Disposition"] = \
                "attachment; filename = \"%s\"" % reqFile # Force downloading.
            req.write(handle.read())
            handle.close()
            del C
        #if
        else :
            #raise fileName
            return apache.HTTP_FORBIDDEN
        return apache.OK
    #if

    # Generate the WSDL file from the MutalyzerService class.
    if req.uri.endswith(".wsdl"):
        servicepath = "http://" + reqPath + "/services"
        client = make_service_client(servicepath, webservice.MutalyzerService())
        req.content_type = 'text/xml'
        req.write(client.server.wsdl(servicepath))
        return apache.OK
    #if

    # By default, use the HTML publisher.
    req.filename = myPath + "/" + req.filename.split('/')[-1]
    return publisher.handler(req)
#handler
