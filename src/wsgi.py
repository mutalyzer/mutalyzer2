#!/usr/bin/python

"""
General WSGI interface.

@todo: Remove handler.py and index.py
@todo: Integrate webservice.py

Public fields:
   -  application ; The WSGI application
"""

import os
import bz2
import web
import site

from cStringIO import StringIO
from simpletal import simpleTALES
from simpletal import simpleTAL

import logging; logging.basicConfig()

# Add /src to Python path
site.addsitedir(os.path.dirname(__file__))

root_dir = os.path.split(os.path.dirname(__file__))[0]
# todo: fix Mutalyzer to not depend on working directory
os.chdir(root_dir)

import Mutalyzer
from Modules import Web
from Modules import Config
from Modules import Output
from Modules import Parser


C = Config.Config()


class render_tal:
    """
    Rendering interface to TAL templates.

    Example:

        render = render_tal('templates')
        render.hello('alice')
    """
    def __init__(self, path, globals={}):
        self.path = path
        self.globals = globals

    def __getattr__(self, name):

        filename = name

        def template(args={}, scheme='html', standalone=False):
            """
            @arg args: Arguments for template.
            @kwarg scheme: One of 'html', 'file'.
            @kwarg standalone: Includes HTML site layout and sets interactive
                               argument for template.
            """
            file = filename
            if scheme == 'html':
                file += '.html'
            path = os.path.join(self.path, file)

            context = simpleTALES.Context()

            context.addGlobal("interactive", not standalone)

            for name, value in self.globals.items():
                context.addGlobal(name, value)

            for name, value in args.items():
                context.addGlobal(name, value)

            templateFile = open(path, 'r')
            template = simpleTAL.compileHTMLTemplate(templateFile)
            templateFile.close()

            # Wrap in site layout
            if scheme == 'html' and not standalone:
                context.addGlobal('sitemacros', template)
                templateFile = open(os.path.join(self.path, 'menu.html'), 'r')
                template = simpleTAL.compileHTMLTemplate(templateFile)
                templateFile.close()

            io = StringIO()
            template.expand(context, io)
            return io.getvalue()

        return template


urls = (
    '/(?:index)?',                'Index',
    '/about',                     'About',
    '/help',                      'Help',
    '/nameGenerator',             'Generator',
    '/webservices',               'Webservices',
    '/check',                     'Check',
    '/syntaxCheck',               'SyntaxCheck',
    '/checkForward',              'CheckForward',
    '/download/(.*\.(?:py|cs))',  'Download',
    '/Reference/(.*)',            'Reference'
)


render = render_tal(os.path.join(root_dir, 'templates'),
                    globals={'version': '2.0&nbsp;&beta;-5',
                             'nomenclatureVersion': '2.0',
                             'releaseDate': '10 Dec 2010',
                             'contactEmail': C.Retriever.email})
app = web.application(urls, globals(), autoreload=False)
session = web.session.Session(app,
                              web.session.DiskStore(os.path.join(root_dir, 'var', 'sessions')),
                              initializer={'variant': None})

# todo: we should probably get rid of Web alltogether
W = Web.Web()
#O = Output.Output(__file__, C.Output)


class Download:
    """
    @todo: This is a potential security hole.
    """
    def GET(self, file):
        # Process the file with TAL and return the content as a downloadable file.
        #file = web.ctx.path.split('/')[-1]
        if not os.path.isfile("templates/" + file):
            raise web.notfound()
        #req.content_type = 'application/octet-stream'
        # Force downloading.
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)
        return getattr(render, file)({'path': web.ctx.homedomain + web.ctx.homepath},
                                     scheme='file')

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
        return render.parse(args)
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
        return render.parse(args)

class Check:
    """
    @todo: These handlers need some documentation.
    """
    def GET(self):
        interactive = True
        i = web.input(mutationName=None)
        if i.mutationName:
            interactive = False
            variant = i.mutationName
        else:
            variant = session.variant
            session.variant = None
        return self.check(variant, interactive=interactive)

    def POST(self):
        i = web.input(mutationName=None)
        return self.check(i.mutationName)

    def check(self, name=None, interactive=True):
        O = Output.Output(__file__, C.Output)

        if name:
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

        return render.check(args, standalone=not interactive)


class CheckForward:
    def GET(self):
        i = web.input()
        session.variant = i.mutationName
        raise web.seeother('check')

# todo: merge the static pages below

class Index:
    def GET(self):
        return render.index()

class About:
    def GET(self):
        return render.about()

class Help:
    def GET(self):
        return render.help()

class Generator:
    def GET(self):
        return render.generator()

class Webservices:
    def GET(self):
        return render.webservices()

# todo:
#   "downloads/" in req.uri
#   "Results" in req.uri:
#   publisher -> index.py
#   everything in index

application = app.wsgifunc()
