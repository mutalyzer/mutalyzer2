#!/usr/bin/env python

"""
General WSGI interface.

The WSGI interface is exposed through the module variable 'application'.
Static files are not handled by this interface and should be served through
the '/base' url prefix separately.

Example Apache/mod_wsgi configuration:

  WSGIScriptAlias / /var/www/mutalyzer/src/wsgi.py
  Alias /base /var/www/mutalyzer/templates/base

You can also use the built-in HTTP server by running this file directly.
Note, however, that static files are not served by this server. A common
pattern is to use Nginx as a proxy and static file server.

Example Nginx configuration (assumes the built-in HTTP server is running on
port 8080):

  server {
    listen 80;
    location /base/ {
      root /var/www/mutalyzer/templates/base;
      if (-f $request_filename) {
        rewrite ^/base/(.*)$  /base/$1 break;
      }
    }
    location / {
      proxy_read_timeout 300;  # 5 minutes
      proxy_pass http://127.0.0.1:8080;
    }
  }

@todo: Integrate webservice.py (http://webpy.org/cookbook/webservice/).
@todo: Move /templates/base to /static for web.py compatibility.
"""


VERSION = '2.0&nbsp;&beta;-7'
NOMENCLATURE_VERSION = '2.0'
RELEASE_DATE = '18 Jan 2011'
WEBSERVICE_LOCATION = '/services'
WSDL_VIEWER = 'templates/wsdl-viewer.xsl'


# Log exceptions to stdout
import logging; logging.basicConfig()

import re
import os
import bz2
import web
import urllib
import site

from lxml import etree
from cStringIO import StringIO
from simpletal import simpleTALES
from simpletal import simpleTAL

# Add /src to Python path
site.addsitedir(os.path.dirname(__file__))

# Todo: Get this from the configuration file
root_dir = os.path.split(os.path.dirname(__file__))[0]
# Todo: Fix Mutalyzer to not depend on working directory
if not __name__ == '__main__':
    os.chdir(root_dir)

import webservice
import Mutalyzer
import VarInfo
from Modules import Config
from Modules import Output
from Modules import Parser
from Modules import Mapper
from Modules import Db
from Modules import Scheduler
from Modules import Retriever
from Modules import File


web.config.debug = True


# Load configuration from configuration file
C = Config.Config()


# URL dispatch table
urls = (
    '/(index)?',            'Static',
    '/(about)',             'Static',
    '/(help)',              'Static',
    '/(faq)',               'Static',
    '/(exercise)',          'Static',
    '/(disclaimer)',        'Static',
    '/(nameGenerator)',     'Static',
    '/(webservices)',       'Static',
    '/documentation',       'Documentation',
    '/snp',                 'Snp',
    '/positionConverter',   'PositionConverter',
    '/Variant_info',        'VariantInfo',
    '/getGS',               'GetGS',
    '/check',               'Check',
    '/syntaxCheck',         'SyntaxCheck',
    '/checkForward',        'CheckForward',
    '/batch([a-zA-Z]+)?',   'BatchChecker',
    '/progress',            'BatchProgress',
    '/Results_(\d+)\.txt',  'BatchResult',
    '/download/([a-zA-Z-]+\.(?:py|cs))',  'Download',
    '/downloads/([a-zA-Z\._-]+)',         'Downloads',
    '/Reference/([\da-zA-Z\._-]+)',       'Reference',
    '/upload',              'Uploader'
)


class render_tal:
    """
    Render interface to TAL templates.

    Example to render /templates/hello.html with parameter 'alice':

        render = render_tal('templates')
        render.hello('alice')
    """
    def __init__(self, path, globals={}):
        """
        @arg path: Path to templates directory.
        @kwarg globals: Dictionary of global template variables.
        """
        self.path = path
        self.globals = globals

    def __getattr__(self, name):
        """
        Returns a template. Call the template to get a render.

        @arg name: Template name (usually a HTML filename without '.html').
        @return: Template render function.
        """
        filename = name

        def template(args={}, scheme='html', standalone=False):
            """
            Template render function.

            If a scheme of 'html' is choosen, the template name is assumed
            to be a filename without its '.html' suffix. Otherwise it is
            assumed to be a full filename.

            The render of the template is wrapped in the HTML site layout
            with menu if scheme is 'html' and standalone is False.

            @arg args: Arguments for template.
            @kwarg scheme: One of 'html', 'file'.
            @kwarg standalone: Includes HTML site layout and sets interactive
                               argument for template.
            @return: Render of template.
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

            # Wrap in site layout with menu
            if scheme == 'html' and not standalone:
                context.addGlobal('sitemacros', template)
                templateFile = open(os.path.join(self.path, 'menu.html'), 'r')
                template = simpleTAL.compileHTMLTemplate(templateFile)
                templateFile.close()

            if scheme == 'html':
                web.header('Content-Type', 'text/html')

            io = StringIO()
            template.expand(context, io)
            return io.getvalue()

        return template
#render_tal


# TAL template render
render = render_tal(os.path.join(root_dir, 'templates'),
                    globals={'version': VERSION,
                             'nomenclatureVersion': NOMENCLATURE_VERSION,
                             'releaseDate': RELEASE_DATE,
                             'contactEmail': C.Retriever.email})

# web.py application
app = web.application(urls, globals(), autoreload=False)

# Sessions are only used by CheckForward (as a hack)
session = web.session.Session(app,
                              web.session.DiskStore(os.path.join(root_dir, 'var', 'sessions')),
                              initializer={'variant': None})


class Download:
    """
    Download file from template directory, processing it with TAL first.
    """
    def GET(self, file):
        """
        @arg file: Filename to download.
        @type file: string

        Be very careful to not call this with anything but an ordinary
        filename. A possible security issue is allowing this method to be
        called with file='../mutalyzer.conf' for example.

        The url routing currently makes sure to only call this with filenames
        of the form [a-zA-Z-]+\.(?:py|cs).
        """
        # Process the file with TAL and return the content as a downloadable file.
        if not os.path.isfile("templates/" + file):
            raise web.notfound()
        # Force downloading
        web.header('Content-Type', 'text/plain')
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)
        return getattr(render, file)({'path': web.ctx.homedomain + web.ctx.homepath},
                                     scheme='file')
#Download


class Downloads:
    """
    Download plain text files from /templates/downloads directory.
    """
    def GET(self, file):
        """
        @arg file: Filename to download.
        @type file: string

        Be very careful to not call this with anything but an ordinary
        filename. A possible security issue is allowing this method to be
        called with file='../../mutalyzer.conf' for example.

        The url routing currently makes sure to only call this with filenames
        of the form [a-zA-Z\._-]+.
        """
        if not os.path.isfile("templates/downloads/" + file):
            raise web.notfound()
        handle = open("templates/downloads/" + file)
        F = File.File(C.File, None)
        web.header('Content-Type', F.getMimeType(handle)[0])
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)
        return handle.read()
#Downloads


class Reference:
    """
    Download reference file from cache.
    """
    def GET(self, file):
        """
        @arg file: Filename to download from cache.
        @type file: string

        Be very careful to not call this with anything but an ordinary
        filename. A possible security issue is allowing this method to be
        called with file='../../mutalyzer.conf' for example.

        The url routing currently makes sure to only call this with filenames
        of the form [a-zA-Z\._-]+.
        """
        fileName = "%s/%s.bz2" % (C.Retriever.cache, file)
        if not os.path.isfile(fileName):
            raise web.notfound()
        handle = bz2.BZ2File(fileName, 'r')
        web.header('Content-Type', 'text/plain')
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)
        return handle.read()
#Reference


class GetGS:
    """
    LOVD bypass to get the correct GeneSymbol incl Transcript variant.

    Used by LOVD to get the correct transcript variant out of a genomic
    record. LOVD uses a genomic reference (NC_?) in combination with a gene
    symbol to pass variant info to mutalyzer. Mutalyzer 1.0 was only using
    the first transcript. LOVD supplies the NM of the transcript needed but
    this was ignored. This helper allows LOVD to get the requested
    transcript variant from a genomic reference.

    @todo: Test this.
    """
    def GET(self):
        """
        Parameters:
        - mutationName: The mutationname without gene symbol.
        - variantRecord: The NM reference of the variant.
        - forward: If set this forwards the request to the name checker.

        @return: Output of name checker if forward is set, otherwise the
                 GeneSymbol with the variant notation as string.
        """
        O = Output.Output(__file__, C.Output)

        i = web.input(mutationName=None, variantRecord=None, forward=None)

        # We are only interested in the legend
        Mutalyzer.process(i.mutationName, C, O)

        legends = O.getOutput("legends")

        # Filter the transcript from the legend
        legends = [l for l in legends if "_v" in l[0]]
        for l in legends:
            if l[1] == i.variantRecord:
                if i.forward:
                    p,a = i.mutationName.split(':')
                    return Check.check(p+'('+l[0]+'):'+a)
                else:
                    web.header('Content-Type', 'text/plain')
                    return l[0]

        web.header('Content-Type', 'text/plain')
        return "Transcript not found"#+`legends`
#GetGS


class SyntaxCheck:
    """
    Syntax checker.
    """
    def GET(self):
        """
        Render syntax checker HTML form.
        """
        args = {
            "variant"       : '',
            "messages"      : [],
            "parseError"    : None,
            "debug"         : ""
        }
        return render.parse(args)

    def POST(self):
        """
        Parse the given variant and render the syntax checker HTML form.

        Parameters:
        - variant: Variant name to check.
        """
        O = Output.Output(__file__, C.Output)
        i = web.input()
        variant = i.variant
        if variant.find(',') >= 0:
            O.addMessage(__file__, 2, "WCOMMASYNTAX",
                         "Comma's are not allowed in the syntax, autofixed")
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
#SyntaxCheck


class Snp:
    """
    SNP converter.

    Convert a dbSNP rs number to HGVS description(s) of the SNP specified on
    the reference sequence(s) used by dbSNP.
    """
    def GET(self):
        """
        Render SNP converter HTML form.
        """
        return self.snp()

    def POST(self):
        """
        Convert to HGVS description(s) and render SNP converter HTML form.

        Parameters:
        - rsId: The dbSNP rs number.
        """
        i = web.input(rsId=None)
        return self.snp(i.rsId)

    def snp(self, rsId=None):
        """
        Convert to HGVS description(s) and render SNP converter HTML form.

        @kwarg rsId: The dbSNP rs number.
        """
        O = Output.Output(__file__, C.Output)

        if rsId :
            O.addMessage(__file__, -1, "INFO", "Received rs%s" % rsId)
            R = Retriever.Retriever(C.Retriever, O, None)
            R.snpConvert(rsId)
            O.addMessage(__file__, -1, "INFO",
                         "Finished processing rs%s" % rsId)
        #if

        args = {
            "snp"      : O.getOutput("snp"),
            "messages" : O.getMessages(),
            "summary"  : O.Summary()[2],
            "lastpost" : rsId
        }

        return render.snp(args)
#Snp


class PositionConverter:
    """
    Convert a variant between genomic and coding positions.
    """
    def GET(self):
        """
        Render position converter HTML form.
        """
        return self.position_converter()

    def POST(self):
        """
        Convert a variant and render position converter HTML form.

        Parameters:
        - build: Human genome build (currently 'hg18' or 'hg19').
        - variant: Variant to convert.
        """
        i = web.input(build='', variant='')
        # Todo: The following is probably a problem elsewhere too.
        # We stringify the variant, because a unicode string crashes
        # Bio.Seq.reverse_complement in Mapper.py:607.
        return self.position_converter(i.build, str(i.variant))

    def position_converter(self, build='', variant=''):
        """
        Convert a variant and render position converter HTML form.

        @kwarg build: Human genome build (currently 'hg18' or 'hg19').
        @kwarg variant: Variant to convert.
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
#PositionConverter


class VariantInfo:
    """
    The I{g.} to I{c.} and vice versa interface for LOVD.

    @todo: Tests.
    """
    def GET(self):
        """
        Run VarInfo and return the result as plain text.

        Parameters:
        - LOVD_ver: The version of the calling LOVD.
        - build: The human genome build (hg19 assumed).
        - acc: The accession number (NM number).
        - var: A description of the variant.
        """
        i = web.input(var='')
        LOVD_ver = i.LOVD_ver
        build = i.build
        acc = i.acc
        var = i.var

        result = VarInfo.main(LOVD_ver, build, acc, var)

        web.header('Content-Type', 'text/plain')

        if LOVD_ver == "2.0-23" : # Obsoleted error messages, remove when possible.
            import re
            return re.sub("^Error \(.*\):", "Error:", result)
        #if
        return result
#VariantInfo


class Check:
    """
    The variant checker.
    """
    def GET(self):
        """
        Render the variant checker HTML form.

        There are two modes of invoking the checker with a GET request:
        1. Provide the 'mutationName' parameter. In this case, the checker is
           called non-interactively, meaning the result is rendered without
           the HTML form, site layout, and menu.
        2. By having a 'variant' value in the session. The value is removed.

        Parameters:
        - mutationName: Variant to check.
        """
        interactive = True
        i = web.input(mutationName=None)
        if i.mutationName:
            # Run checker non-interactively
            interactive = False
            variant = i.mutationName
        else:
            # Run checker if session.variant is not None
            variant = session.variant
            session.variant = None
        return self.check(variant, interactive=interactive)

    def POST(self):
        """
        Run the name checker and render the variant checker HTML form.

        Parameters:
        - mutationName: Variant to check.
        """
        i = web.input(mutationName=None)
        return self.check(i.mutationName)

    @staticmethod
    def check(name=None, interactive=True):
        """
        Render the variant checker HTML form. If the name argument is given,
        run the name checker.

        @kwarg name: Variant to check.
        @kwarg interactive: Run interactively, meaning we wrap the result in
                            the site layout and include the HTML form.
        """
        O = Output.Output(__file__, C.Output)

        if name:
            O.addMessage(__file__, -1, "INFO", "Received variant %s" % name)
            # Todo: The following is probably a problem elsewhere too.
            # We stringify the variant, because a unicode string crashes
            # Bio.Seq.reverse_complement in Mapper.py:607.
            RD = Mutalyzer.process(str(name), C, O)
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

        def urlEncode(descriptions):
            """
            @todo: This should probably be done in the template.

            @arg descriptions:
            @type descriptions: list

            @return: urlEncode descriptions???????????????
            @rtype: list
            """
            newDescr = []
            for i in descriptions :
                newDescr.append([i, urllib.quote(i)])
            return newDescr

        args = {
            "lastpost"           : name,
            "messages"           : O.getMessages(),
            "summary"            : summary,
            "parseError"         : pe,
            "errors"             : errors,
            "genomicDescription" : urlEncode([genomicDescription])[0] if genomicDescription else "",
            "chromDescription"   : O.getIndexedOutput("genomicChromDescription", 0),
            "genomicDNA"         : genomicDNA,
            "visualisation"      : O.getOutput("visualisation"),
            "descriptions"       : urlEncode(O.getOutput("descriptions")),
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

        # Todo: This shouldn't really be necessary
        del O

        return render.check(args, standalone=not interactive)
#Check


class CheckForward:
    """
    Set the given variant in the session and redirect to the name checker.

    @todo: Cleaner solution (one without using a session variable).
    """
    def GET(self):
        """
        Set the 'variant' session value to the given variant and redirect
        to the name checker (where we will arrive by a GET request).

        Parameters:
        - mutationName: Variant to set in the session.
        """
        i = web.input(mutationName=None)
        session.variant = i.mutationName
        raise web.seeother('check')
#CheckForward


class BatchProgress:
    """
    Batch jobs progress viewer.

    Used from the 'batch' template by AJAX to get the progress of a batch
    job.
    """
    def GET(self):
        """
        Progress for a batch job.

        Parameters:
        - jobID: ID of the job to return progress for.
        - totalJobs: Total number of entries in this job.
        - ajax: If set, return plain text result.

        @todo: The 'progress' template does not exist.
        """
        O = Output.Output(__file__, C.Output)

        attr = {"percentage": 0}

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
            web.header('Content-Type', 'text/plain')
            return ret
        else:
            #Return progress html page
            return render.progress(attr)
#BatchProgress


class BatchChecker:
    """
    Run batch jobs.
    """
    def GET(self, batchType=None):
        """
        Render batch checker HTML form.

        @kwarg batchType: Type of batch job.
        """
        return self.batch(batchType=batchType)

    def POST(self, bt=None):
        """
        Run batch jobs and render batch checker HTML form.

        @kwarg bt: Not used, batch type in URL.

        Parameters:
        - batchEmail: Email address to mail results to.
        - batchFile: Uploaded file with batch job entries.
        - arg1: Additional argument. Currently only used if batchType is
                'PositionConverter', denoting the human genome build.
        - batchType: Type of batch job to run. One of 'NameChecker' (default),
                     'SyntaxChecker', or 'PositionChecker'.
        """
        i = web.input(batchEmail=None, batchFile={}, arg1='',
                      batchType=None)
        return self.batch(email=i.batchEmail, inFile=i.batchFile, arg1=i.arg1,
                          batchType=i.batchType)

    def batch(self, email=None, inFile=None, arg1='', batchType=None):
        """
        Run batch jobs and render batch checker HTML form. The batch jobs are
        added to the database by the scheduler and ran by the BatchChecker
        daemon.

        @kwarg email: Email address to mail results to.
        @kwarg inFile: Uploaded file with batch job entries.
        @kwarg arg1: Additional argument. Currently only used if batchType is
                     'PositionConverter', denoting the human genome build.
        @kwarg batchType: Type of batch job to run. One of 'NameChecker'
                          (default), 'SyntaxChecker', or 'PositionChecker'.
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

        # Todo: I think this test is kindof bogus
        def isEMail(a):
            return bool(
                re.match("^[a-zA-Z0-9._%-]+@[a-zA-Z0-9._%-]+.[a-zA-Z]{2,6}$",
                         a))

        # Note: A FieldStorage instance (like inFile) seems to always test
        # to the truth value False, so 'if inFile: ...' is not useful.

        if email and isEMail(email) and not inFile == None and inFile.file:
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
#BatchChecker


class BatchResult:
    """
    Download result from the batch checker.
    """
    def GET(self, result):
        """
        Return raw content (for batch checker results).

        @arg result: Result identifier.
        @type result: string

        Be very careful to not call this with anything but an ordinary
        filename. A possible security issue is allowing this method to be
        called with file='../../mutalyzer.conf' for example.

        The url routing currently makes sure to only call this with filenames
        of the form \d+.
        """
        file = 'Results_%s.txt' % result
        handle = open(os.path.join(C.Scheduler.resultsDir, file))
        web.header('Content-Type', 'text/plain')
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)
        return handle.read()
#BatchResult


def _checkInt(inpv, refname):
    """
    Remove [,.-] from inpv and try to convert the result to an integer value.
    Raise InputException if the conversion fails.

    This private function is used by Uploader.

    @arg inpv: Input value to convert.
    @arg refname: Name of input value.

    @return: Input value converted to an integer value.

    @raise InputException: If the converting to an integer value fails.
    """
    inpv = inpv.replace(',','').replace('.','').replace('-','')
    try:
        return int(inpv)
    except ValueError, e:
        raise InputException("Expected an integer in field: %s" % refname)
#_checkInt


class InputException(Exception):
    """
    This exception is raised by Uploader.
    """
    pass
#InputException


class Uploader:
    """
    Reference sequence uploader.

    Upload or retrieve a reference sequence.

    @todo: Test this class.
    """
    def GET(self):
        """
        Render reference sequence uploader form.
        """
        maxUploadSize = C.Retriever.maxDldSize
        UD, errors = "", []
        args = {
            "UD"      : UD,
            "maxSize" : float(maxUploadSize) / 1048576,
            "errors"  : errors
        }
        return render.gbupload(args)

    def POST(self):
        """
        Render reference sequence uploader form and handle a reference
        sequence upload or retrieval.

        This handler has four methods:
        1. The reference sequence file is a local file.
        2. The reference sequence file can be found at the following URL.
        3. Retrieve part of the reference genome for a (HGNC) gene symbol.
        4. Retrieve a range of a chromosome.

        Parameters:
        - invoermethode: Input method. One of 'file', 'url', 'gene', 'chr'.

        Depending on the input method, additional parameters are expected.

        Parameters (method 'file'):
        - bestandsveld: Reference sequence file to upload.

        Parameters (method 'url'):
        - urlveld: URL of reference sequence file to upload.

        Parameters (method 'gene'):
        - genesymbol: Gene symbol.
        - organism: Organism.
        - fiveutr: Number of 5' flanking nucleotides.
        - threeutr: Number of 3' flanking nucleotides.

        Parameters (method 'chr'):
        - chracc: Chromosome Accession Number.
        - start: Start position.
        - stop: Stop position.
        - orientation: Orientation.
        """
        maxUploadSize = C.Retriever.maxDldSize

        O = Output.Output(__file__, C.Output)
        D = Db.Cache(C.Db)
        R = Retriever.GenBankRetriever(C.Retriever, O, D)

        UD, errors = "", []

        i = web.input(invoermethode='', bestandsveld={}, urlveld='',
                      genesymbol='', organism='', fiveutr='', threeutr='',
                      chracc='', start='', stop='', orientation='')

        try:
            if i.invoermethode == "file" :
                if not 'Content-Length' in web.ctx.environ:
                    web.ctx.status = '411 Length required'
                    return 'Content length required.'
                #if
                if int(web.ctx.environ['Content-Length']) > maxUploadSize :
                    web.ctx.status = '413 Request entity too large'
                    return 'Upload limit exceeded.'
                #if
                if not i.bestandsveld == None and i.bestandsveld.file:
                    UD = R.uploadrecord(i.bestandsveld.file.read())
                else:
                    raise web.badrequest()
            #if
            elif i.invoermethode == "url" :
                UD = R.downloadrecord(i.urlveld)
            #if
            elif i.invoermethode == "gene" :
                geneName = i.genesymbol
                organism = i.organism
                upStream = _checkInt(i.fiveutr,
                        "5' flanking nucleotides")
                downStream = _checkInt(i.threeutr,
                        "3' flanking nucleotides")
                UD = R.retrievegene(geneName, organism, upStream, downStream)
            #if
            elif i.invoermethode == "chr" :
                accNo = i.chracc
                start = _checkInt(i.start,
                        "Start position")
                stop = _checkInt(i.stop,
                        "Stop position")
                orientation = int(i.orientation)
                UD = R.retrieveslice(accNo, start, stop, orientation)
            #if
            else:
                #unknown "invoermethode"
                raise InputException("Wrong method selected")
        except InputException, e:
            #DUMB USERS
            errors.append(e)
        finally:
            if not UD:
                #Something went wrong
                errors += ["The request could not be completed"]
                errors.extend(O.getMessages())

        args = {
            "UD"      : UD,
            "maxSize" : float(maxUploadSize) / 1048576,
            "errors"  : errors
        }
        return render.gbupload(args)
#Uploader


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
        @todo: Use configuration value for .xsl location.
        @todo: Cache this transformation.
        """
        url = web.ctx.homedomain + web.ctx.homepath + WEBSERVICE_LOCATION
        wsdl_handle = StringIO(webservice.soap_application.get_wsdl(url))
        xsl_handle = open(WSDL_VIEWER, 'r')
        wsdl_doc = etree.parse(wsdl_handle)
        xsl_doc = etree.parse(xsl_handle)
        transform = etree.XSLT(xsl_doc)
        web.header('Content-Type', 'text/html')
        return str(transform(wsdl_doc))
#Documentation


class Static:
    """
    Static page, just render a TAL template on GET.
    """
    def GET(self, page=None):
        """
        Render a TAL template as HTML.

        @kwarg page: Page name to render. A TAL template with this name must
                     exist. Special case is a page of None, having the same
                     effect as 'index'.
        @type page: string

        Be careful to only call this method with an argument that is a simple
        template name. For example, make sure this is not called with page
        value '../forbidden'. This check is implemented in the url routing.
        """
        if not page:
            page = 'index'
        return getattr(render, page)()


if __name__ == '__main__':
    # Todo: Setting the working directory probably doesn't work
    # Usage:
    #   ./src/wsgi.py [port]
    app.run()
else:
    # WSGI application
    application = app.wsgifunc()
