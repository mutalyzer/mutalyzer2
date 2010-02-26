#!/usr/bin/python

"""
    General handler for mod_python.

    This handler calls the webservices handler, the HTML publisher or handles
    a request itself, depending on keywords in the URI of the request.

    Public methods:
        handler(req) ; The handler.
"""

import os
from ZSI import dispatch

# This is a workaround for pydoc. Aparently it crashes when a module can not be
# imported.
try : 
    from mod_python import apache, publisher
except ImportError :
    pass

from Modules import Web
from Interfaces import webservice

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
    myPath = os.path.dirname(__file__) + "/.."
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
    if ".js" in req.uri or "base" in req.uri :
        req.content_type = 'text/html'
        req.write(W.read("templates/", req))
        return apache.OK
    #if

    # Process the file with TAL and return the content as a downloadable file.
    if ".py" in req.uri :
        reqFile = req.uri.split('/')[-1]
        req.content_type = 'text/plain'
        req.headers_out["Content-Disposition"] = \
          "attachment; filename = \"%s\"" % reqFile # To force downloading.
        args = {"path" : reqPath}                   # Replace the path variable.
        req.write(W.tal("HTML", "templates/" + reqFile, args))
        return apache.OK
    #if

    # Generate the WSDL file from the MutalyzerService class.
    if ".wsdl" in req.uri :
        from soaplib.client import make_service_client

        servicepath = "http://" + reqPath + "/services"
        client = make_service_client(servicepath, webservice.MutalyzerService())
        req.content_type = 'text/xml'
        req.write(client.server.wsdl(servicepath))
        return apache.OK
    #if

    # By default, use the HTML publisher.
    req.filename = myPath + "/Interfaces/" + req.filename.split('/')[-1]
    return publisher.handler(req)
#handler
