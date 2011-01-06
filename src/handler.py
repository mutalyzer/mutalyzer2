#!/usr/bin/python

"""
General handler for mod_python.

This handler calls the webservices handler, the HTML publisher or handles
a request itself, depending on keywords in the URI of the request.

Public methods:
   -  handler(req) ; The handler.
   
@requires: os
@requires: bz2
@requires: mod_python.apache
@requires: mod_python.publisher
@requires: Modules.Web
@requires: Modules.Config
@requires: Modules.File
"""

import os
import bz2
from mod_python import apache, publisher

from Modules import Web
from Modules import Config
from Modules import File

def handler(req):
    """
    Handle a request passed to us by mod_python.

    Keywords in the URI of the request are used to decide what to do:
        - ".js"      ; Return the raw content of the file (to include JavaScript
                       from an HTML file).
        - ".py"      ; Return the content as a downloadable file after it has
          ".cs"        been processed by TAL (to generate webservice client
                       files).

    By default, the HTML publisher is used for normal HTML files and other
    unhandled requests.

    @arg req: The request
    @type req: string

    @return: An apache return code, either generated in this function
    itself, or by the publisher handler which handles normal HTML requests
    @rtype: integer
    """

    # Figure out where this program is located and go to the parent directory.
    myPath = os.path.dirname(__file__)
    os.chdir(myPath + "/..")

    reqPath = req.hostname + req.uri.rsplit('/', 1)[0]

    W = Web.Web()

    # Make all txt files in downloads directory downloadable
    if "downloads/" in req.uri :
        reqFile = req.uri.rsplit('/')[-1]
        req.headers_out["Content-Disposition"] = \
            'attachment; filename="%s"' % reqFile # To force downloading.
        handle = open("templates/downloads/" + reqFile)
        C = Config.Config()
        F = File.File(C.File, None)
        req.content_type = F.getMimeType(handle)[0]
        req.write(handle.read())
        handle.close()
        return apache.OK
    #if
