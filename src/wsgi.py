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
import urllib2
from lxml import etree
import site
#import pydoc

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
#import webservice
from Modules import Web
from Modules import Config
from Modules import Output
from Modules import Parser
from Modules import Mapper
from Modules import Db
from Modules import Scheduler
from Modules import Retriever
from Modules import File


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
    '/documentation',             'Documentation',
    '/snp',                       'Snp',
    '/positionConverter',         'PositionConverter',
    '/check',                     'Check',
    '/syntaxCheck',               'SyntaxCheck',
    '/checkForward',              'CheckForward',
    '/batch([a-zA-Z]+)?',         'BatchChecker',
    '/progress',                  'BatchProgress',
    '/Results_(\d+)\.txt',        'BatchResult',
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


class Snp:
    """
    @todo: Some documentation.
    """
    def GET(self):
        return self.snp()

    def POST(self):
        i = web.input(rsId=None)
        return self.snp(i.rsId)

    def snp(self, rsId=None):
        """
        @todo: documentation

        @arg req: the HTTP request
        @type req: object
        @return: compiled TAL template
        @rtype: object
        """
        O = Output.Output(__file__, C.Output)

        if rsId :
            O.addMessage(__file__, -1, "INFO", "Received rs%s" % rsId)
            R = Retriever.Retriever(C.Retriever, O, None)
            R.snpConvert(rsId)
            O.addMessage(__file__, -1, "INFO", "Finished processing rs%s" % rsId)
        #if

        args = {
            "snp"      : O.getOutput("snp"),
            "messages" : O.getMessages(),
            "summary"  : O.Summary()[2],
            "lastpost" : rsId
        }

        return render.snp(args)


class PositionConverter:
    """
    @todo: Some documentation.
    """
    def GET(self):
        return self.position_converter()

    def POST(self):
        i = web.input(build='', variant='')
        # We stringify the variant, because a unicode string crashes
        # Bio.Seq.reverse_complement in Mapper.py:607.
        return self.position_converter(i.build, str(i.variant))

    def position_converter(self, build='', variant=''):
        """
        @arg req:
        @type req:

        @todo: documentation
        """
        O = Output.Output(__file__, C.Output)

        avail_builds = C.Db.dbNames[::-1]

        if build :
            avail_builds.remove(build)
            avail_builds.insert(0, build)
        #if

        attr = {
            "avail_builds" : avail_builds,
            "variant"      : variant,
            "gName"        : "",
            "cNames"       : [],
            "messages"     : [],
            "errors"       : [],
            "debug"        : []
        }

        if build and variant:
            converter = Mapper.Converter(build, C, O)

            #Conver chr accNo to NC number
            variant = converter.correctChrVariant(variant)

            if variant :
                if not(":c." in variant or ":g." in variant):
                    #Bad name
                    P = Parser.Nomenclatureparser(O)
                    parsetree = P.parse(variant)
                #if

                if ":c." in variant:
                    # Do the c2chrom dance
                    variant = converter.c2chrom(variant)
                if variant and ":g." in variant:
                    # Do the g2c dance
                    variants = converter.chrom2c(variant, "dict")
                    if variants:
                        attr["gName"] = variant
                        output = ["%-10s:\t%s" % (key[:10], "\n\t\t".join(value))\
                                  for key, value in variants.items()]
                        attr["cNames"].extend(output)
                    #if
                #if
            #if

            attr["errors"].extend(O.getMessages())
        return render.converter(attr)


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


class BatchProgress:
    def GET(self):
        """
        Progress page for batch runs
        @todo: documentation
        """
        O = Output.Output(__file__, C.Output)

        attr = {"percentage"    : 0}

        i = web.input(ajax=None)
        try:
            jobID = int(i.jobID)
            total = int(i.totalJobs)
        except Exception, e:
            return
        D = Db.Batch(C.Db)
        left = D.entriesLeftForJob(jobID)
        percentage = int(100 - (100 * left / float(total)))
        if i.ajax:
            if percentage == 100:
                #download url, check if file still exists
                ret = "OK"
            else:
                ret = percentage
            return ret
        else:
            #Return progress html page
            return render.progress(attr)


class BatchChecker:
    def GET(self, batchType=None):
        return self.batch(batchType=batchType)

    def POST(self, bt=None):
        i = web.input(batchEmail=None, batchFile={}, arg1='',
                      batchType=None)
        return self.batch(email=i.batchEmail, inFile=i.batchFile, arg1=i.arg1,
                          batchType=i.batchType)

    def batch(self, email=None, inFile=None, arg1='', batchType=None):
        """
        Batch function to add batch jobs to the Database

        @arg batchType: Type of the batch job
        @type batchType: string
        """
        O = Output.Output(__file__, C.Output)

        attr = {"messages"      : [],
                "errors"        : [],
                "debug"         : [],
                "batchTypes"    : ["NameChecker",
                                   "SyntaxChecker",
                                   "PositionConverter"],
                "hideTypes"     : batchType and 'none' or '',
                "selected"      : "0",
                "batchType"     : batchType or "",
                "avail_builds"  : C.Db.dbNames[::-1],
                "jobID"         : None,
                "totalJobs"     : None
        }

        #Make sure the correct page is displayed for an entrypoint
        if not batchType: batchType = 'NameChecker'

        if batchType in attr["batchTypes"]:
            attr["selected"] = str(attr["batchTypes"].index(batchType))

        # Note: A FieldStorage instance (like inFile) seems to always test
        # to the truth value False, so 'if inFile: ...' is not useful.

        if email and W.isEMail(email) and not inFile == None and inFile.file:
            D = Db.Batch(C.Db)
            S = Scheduler.Scheduler(C.Scheduler, D)
            FileInstance = File.File(C.File, O)

            # Generate the fromhost URL from which the results can be fetched
            fromHost = web.ctx.homedomain + web.ctx.homepath + '/'
            #fromHost = "http://%s%s" % (
            #    req.hostname, req.uri.rsplit("/", 1)[0]+"/")

            job = FileInstance.parseBatchFile(inFile.file)
            if job is None:
                O.addMessage(__file__, 4, "PRSERR", "Could not parse input"
                             " file, please check your file format.")
            else:
                #TODO: Add Binair Switches to toggle some events
                attr["jobID"] =\
                              S.addJob("BINSWITHCES", email, job, fromHost, batchType, arg1)
                attr["totalJobs"] = len(job) or 1
                attr["messages"].append("Your file has been parsed and the job"
                                        " is scheduled, you will receive an email when the job is "
                                        "finished.")

            attr["errors"].extend(O.getMessages())

        return render.batch(attr)


class BatchResult:
    def GET(self, result):
        """
        Return raw content (for batch checker results).
        """
        file = 'Results_%s.txt' % result
        handle = open(os.path.join(C.Scheduler.resultsDir, file))
        web.header('Content-Type', 'text/plain')
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)
        return handle.read()


class Documentation:
    def GET(self):
        """
        HTML documentation for the webservice.

        Generate the documentation by a XSL transform of the WSDL document.
        The XSL transformation used is from Tomi Vanek:

          http://tomi.vanek.sk/index.php?page=wsdl-viewer

        We apply a small patch to this transformation to show newlines in
        the SOAP method docstrings:

          Around line 1195, the description <div>, replace
          '<div class="value">' by '<div class="value documentation">'.

          In the style sheet, add:
            .documentation { white-space: pre-line; }

        @todo: Use some configuration setting for the location of the
               webservice.
        @todo: Cache this transformation.
        """
        wsdl_url = web.ctx.homedomain + web.ctx.homepath + '/service?wsdl'
        wsdl_handle = urllib2.urlopen(wsdl_url)
        xsl_handle = open('wsdl-viewer.xsl', 'r')
        wsdl_doc = etree.parse(wsdl_handle)
        xsl_doc = etree.parse(xsl_handle)
        transform = etree.XSLT(xsl_doc)
        return str(transform(wsdl_doc))


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
#   publisher -> index.py
#   everything in index

application = app.wsgifunc()
