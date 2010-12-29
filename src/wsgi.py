#!/usr/bin/python

"""
General WSGI interface.

This handler calls the webservices handler, the HTML publisher or handles
a request itself, depending on keywords in the URI of the request.

Public fields:
   -  application ; The WSGI application

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

import sys
sys.path.append(os.path.dirname(__file__))

import web

from Modules import Web
from Modules import Config
from Modules import File

urls = (
    '/',               'Index',
    '.*\.(?:py|cs)',   'Download'
)

app = web.application(urls, globals(), autoreload=False)

# Figure out where this program is located and go to the parent directory.
myPath = os.path.dirname(__file__)
os.chdir(myPath + "/..")

W = Web.Web()

# todo: "downloads/" in req.uri

class Index:
    def GET(self):
        return W.tal("HTML", "templates/index.html", {})

class Download:
    def GET(self):
        # Process the file with TAL and return the content as a downloadable file.
        file = web.ctx.path.split('/')[-1]
        #req.content_type = 'application/octet-stream'
        # Force downloading.
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)
        return W.tal("TEXT", "templates/" + file, {'path': web.ctx.homedomain + web.ctx.homepath})

# todo: "Results" in req.uri:

# todo: "Reference" in req.uri:

# todo: publisher -> index.py

application = app.wsgifunc()
