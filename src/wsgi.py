#!/usr/bin/python

"""
General WSGI interface.

This handler calls the webservices handler, the HTML publisher or handles
a request itself, depending on keywords in the URI of the request.

Public fields:
   -  application ; The WSGI application
"""

import os
import bz2
import web
import site

# Add /src to Python path
site.addsitedir(os.path.dirname(__file__))

# todo: fix Mutalyzer to not depend on working directory
os.chdir(os.path.split(os.path.dirname(__file__))[0])

import Mutalyzer
from Modules import Web
from Modules import Config
from Modules import Output
from Modules import Parser

urls = (
    '/(?:index)?',                'Index',
    '/about',                     'About',
    '/help',                      'Help',
    '/nameGenerator',             'Generator',
    '/webservices',               'Webservices',
    '/check',                     'Check',
    '/syntaxCheck',               'SyntaxCheck',
    '/download/(.*\.(?:py|cs))',  'Download',
    '/Reference/(.*)',            'Reference'
)

app = web.application(urls, globals(), autoreload=False)

# todo: we should probably get rid of Web alltogether
W = Web.Web()
C = Config.Config()
#O = Output.Output(__file__, C.Output)

class Download:
    def GET(self, file):
        # Process the file with TAL and return the content as a downloadable file.
        #file = web.ctx.path.split('/')[-1]
        if not os.path.isfile("templates/" + file):
            raise web.notfound()
        #req.content_type = 'application/octet-stream'
        # Force downloading.
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)
        return W.tal("TEXT", "templates/" + file, {'path': web.ctx.homedomain + web.ctx.homepath})

class Reference:
    def GET(self, file):
        fileName = "%s/%s.bz2" % (C.Retriever.cache, file)
        if not os.path.isfile(fileName):
            raise web.notfound()
        handle = bz2.BZ2File(fileName, 'r')
        web.header('Content-Type', 'text/plain')
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)
        return handle.read()

class SyntaxCheck:
    def GET(self):
        args = {
            "variant"       : '',
            "messages"      : [],
            "parseError"    : None,
            "debug"         : ""
        }
        return W.tal("HTML", "templates/parse.html", args)
    def POST(self):
        O = Output.Output(__file__, C.Output)
        i = web.input()
        variant = i.variant
        if variant.find(',') >= 0:
            O.addMessage(__file__, 2, "WCOMMASYNTAX", "Comma's are not allowed in the syntax, autofixed")
            variant = variant.replace(',', '')
            #args["variant"]=variant
        P = Parser.Nomenclatureparser(O)
        parsetree = P.parse(variant)
        pe = O.getOutput("parseError")
        if pe: pe[0] = pe[0].replace('<', "&lt;")
        args = {
            "variant"       : variant,
            "messages"      : O.getMessages(),
            "parseError"    : pe,
            "debug"         : ""
        }
        del O
        return W.tal("HTML", "templates/parse.html", args)

class Check:
    def GET(self):
        args = {
            "lastpost"           : None,
            "messages"           : [],
            "summary"            : '',
            "parseError"         : None,
            "errors"             : [],
            "genomicDescription" : '',
            "chromDescription"   : '',
            "genomicDNA"         : '',
            "visualisation"      : '',
            "descriptions"       : '',
            "protDescriptions"   : '',
            "oldProtein"         : '',
            "altStart"           : '',
            "altProtein"         : '',
            "newProtein"         : '',
            "exonInfo"           : '',
            "cdsStart_g"         : '',
            "cdsStart_c"         : '',
            "cdsStop_g"          : '',
            "cdsStop_c"          : '',
            "restrictionSites"   : '',
            "legends"            : '',
            "reference"          : ''
        }
        args["interactive"] = True
        return W.tal("HTML", "templates/check.html", args)

    def POST(self):
        O = Output.Output(__file__, C.Output)
        i = web.input()
        name = i.mutationName
        # todo: load from session
        O.addMessage(__file__, -1, "INFO", "Received variant %s" % name)
        RD = Mutalyzer.process(name, C, O)
        O.addMessage(__file__, -1, "INFO", "Finished processing variant %s" % \
                     name)

        errors, warnings, summary = O.Summary()
        recordType = O.getIndexedOutput("recordType",0)
        reference = O.getIndexedOutput("reference", 0)
        if recordType == "LRG" :
            reference = reference + ".xml" if reference else ""
        else :
            reference = reference + ".gb" if reference else ""

        pe = O.getOutput("parseError")
        if pe :
            pe[0] = pe[0].replace('<', "&lt;")

        genomicDNA = True
        if O.getIndexedOutput("molType", 0) == 'n' :
            genomicDNA = False

        genomicDescription = O.getIndexedOutput("genomicDescription", 0)

        args = {
            "lastpost"           : name,
            "messages"           : O.getMessages(),
            "summary"            : summary,
            "parseError"         : pe,
            "errors"             : errors,
            "genomicDescription" : W.urlEncode([genomicDescription])[0] if genomicDescription else "",
            "chromDescription"   : O.getIndexedOutput("genomicChromDescription", 0),
            "genomicDNA"         : genomicDNA,
            "visualisation"      : O.getOutput("visualisation"),
            "descriptions"       : W.urlEncode(O.getOutput("descriptions")),
            "protDescriptions"   : O.getOutput("protDescriptions"),
            "oldProtein"         : O.getOutput("oldProteinFancy"),
            "altStart"           : O.getIndexedOutput("altStart", 0),
            "altProtein"         : O.getOutput("altProteinFancy"),
            "newProtein"         : O.getOutput("newProteinFancy"),
            "exonInfo"           : O.getOutput("exonInfo"),
            "cdsStart_g"         : O.getIndexedOutput("cdsStart_g", 0),
            "cdsStart_c"         : O.getIndexedOutput("cdsStart_c", 0),
            "cdsStop_g"          : O.getIndexedOutput("cdsStop_g", 0),
            "cdsStop_c"          : O.getIndexedOutput("cdsStop_c", 0),
            "restrictionSites"   : O.getOutput("restrictionSites"),
            "legends"            : O.getOutput("legends"),
            "reference"          : reference
        }

        # todo: this shouldn't really be necessary
        del O

        # todo: there was support for non-interactive usage?? (by GET)
        args["interactive"] = True
        return W.tal("HTML", "templates/check.html", args)

# todo: merge the static pages below

class Index:
    def GET(self):
        return W.tal("HTML", "templates/index.html", {})

class About:
    def GET(self):
        return W.tal("HTML", "templates/about.html", {})

class Help:
    def GET(self):
        return W.tal("HTML", "templates/help.html", {})

class Generator:
    def GET(self):
        return W.tal("HTML", "templates/generator.html", {})

class Webservices:
    def GET(self):
        return W.tal("HTML", "templates/webservices.html", {})

# todo:
#   "downloads/" in req.uri
#   "Results" in req.uri:
#   publisher -> index.py
#   everything in index

application = app.wsgifunc()
