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


WEBSERVICE_LOCATION = '/services'
WSDL_VIEWER = 'templates/wsdl-viewer.xsl'


# WSGI applications should never print anything to stdout. We redirect to
# stderr, but eventually Mutalyzer should be fixed to never just 'print'
# anything.
# http://code.google.com/p/modwsgi/wiki/DebuggingTechniques
import sys
sys.stdout = sys.stderr

# Log exceptions to stdout
import logging; logging.basicConfig()

import re
import os
import bz2
import web
import urllib

from lxml import etree
from cStringIO import StringIO
from simpletal import simpleTALES
from simpletal import simpleTAL

import mutalyzer
from mutalyzer import util
from mutalyzer.config import Config
from mutalyzer.grammar import Grammar
from mutalyzer import webservice
from mutalyzer import variantchecker
from mutalyzer.output import Output
from mutalyzer import Mapper
from mutalyzer import Db
from mutalyzer import Scheduler
from mutalyzer import Retriever
from mutalyzer import File


web.config.debug = True


# Load configuration from configuration file
config = Config()


# URL dispatch table
urls = (
    '',                     'RedirectHome',
    '/(index)?',            'Static',
    '/(about)',             'Static',
    '/(help)',              'Static',
    '/(faq)',               'Static',
    '/(exercise)',          'Static',
    '/(disclaimer)',        'Static',
    '/(nameGenerator)',     'Static',
    '/(webservices)',       'Static',
    '/(webservdoc)',        'Static',
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
    '/download/([a-zA-Z-]+\.(?:py|cs|php|rb))', 'Download',
    '/downloads/([a-zA-Z\._-]+)',               'Downloads',
    '/Reference/([\da-zA-Z\._-]+)',             'Reference',
    '/upload',              'Uploader'
)


class render_tal:
    """
    Render interface to TAL templates.

    Example to render /templates/hello.html with parameter 'alice':

        >>> render = render_tal('templates')
        >>> render.hello('alice')
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
render = render_tal(os.path.join(mutalyzer.package_root(), 'templates'),
                    globals={
    'version': mutalyzer.__version__,
    'nomenclatureVersion': mutalyzer.NOMENCLATURE_VERSION,
    'releaseDate': mutalyzer.__date__,
    'release': mutalyzer.RELEASE,
    'contactEmail': config.Retriever.email})

# web.py application
app = web.application(urls, globals(), autoreload=False)


class RedirectHome:
    """
    Permanent redirect to the homepage.
    """
    def GET(self):
        """
        Redirect to / and include the query string.
        """
        raise web.redirect('/' + web.ctx.query)

    def POST(self):
        """
        Redirect to / and include the query string.
        """
        raise web.redirect('/' + web.ctx.query)


class Download:
    """
    Download file from template directory, formatting it first.
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
        file_path = os.path.join(mutalyzer.package_root(), 'templates', file)
        if not os.path.isfile(file_path):
            raise web.notfound()
        content = open(file_path, 'r').read()
        # Force downloading
        web.header('Content-Type', 'text/plain')
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)
        # We use new style string formatting (available from Python 2.6)
        # http://www.python.org/dev/peps/pep-3101/
        return content.format(path=web.ctx.homedomain + web.ctx.homepath)
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
        file_path = os.path.join(mutalyzer.package_root(),
                                 'templates', 'downloads', file)
        if not os.path.isfile(file_path):
            raise web.notfound()
        handle = open(file_path)
        F = File.File(config.File, None)
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
        file_path = os.path.join(config.Retriever.cache, '%s.bz2' % file)
        if not os.path.isfile(file_path):
            raise web.notfound()
        handle = bz2.BZ2File(file_path, 'r')
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
        O = Output(__file__, config.Output)

        i = web.input(mutationName=None, variantRecord=None, forward=None)

        # Todo: The following is probably a problem elsewhere too.
        # We stringify the variant, because a unicode string crashes
        # Bio.Seq.reverse_complement in Mapper.py:607.

        # We are only interested in the legend
        #Mutalyzer.process(str(i.mutationName), config, O)
        variantchecker.check_variant(str(i.mutationName), config, O)

        legends = O.getOutput("legends")

        # Filter the transcript from the legend
        legends = [l for l in legends if "_v" in l[0]]
        for l in legends:
            if l[1] == i.variantRecord:
                if i.forward:
                    p,a = i.mutationName.split(':')
                    return Check.check(p+'('+l[0]+'):'+a, interactive=False)
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
        output = Output(__file__, config.Output)
        i = web.input()
        variant = i.variant
        if variant.find(',') >= 0:
            output.addMessage(__file__, 2, "WCOMMASYNTAX",
                         "Comma's are not allowed in the syntax, autofixed.")
            variant = variant.replace(',', '')
            #args["variant"]=variant
        grammar = Grammar(output)
        parsetree = grammar.parse(variant)
        pe = output.getOutput("parseError")
        if pe: pe[0] = pe[0].replace('<', "&lt;")
        args = {
            "variant"       : variant,
            "messages"      : map(util.message_info, output.getMessages()),
            "parseError"    : pe,
            "debug"         : ""
        }
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

    def snp(self, rs_id=None):
        """
        Convert {rs_id} to HGVS description(s) and render SNP converter HTML
        form.

        @kwarg rs_id: The dbSNP rs number (including 'rs' prefix).
        @type rs_id: string
        """
        output = Output(__file__, config.Output)

        descriptions = []

        if rs_id:
            output.addMessage(__file__, -1, 'INFO', 'Received %s' % rs_id)
            retriever = Retriever.Retriever(config.Retriever, output, None)
            descriptions = retriever.snpConvert(rs_id)
            output.addMessage(__file__, -1, 'INFO',
                              'Finished processing %s' % rs_id)

        args = {
            'snp'      : descriptions,
            'messages' : map(util.message_info, output.getMessages()),
            'summary'  : output.Summary()[2],
            'lastpost' : rs_id
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
        output = Output(__file__, config.Output)

        avail_builds = config.Db.dbNames[::-1]

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
            "posted"       : build and variant
        }

        if build and variant:
            converter = Mapper.Converter(build, config, output)

            #Convert chr accNo to NC number
            variant = converter.correctChrVariant(variant)

            if variant :
                if not(":c." in variant or ":g." in variant):
                    #Bad name
                    grammar = Grammar(output)
                    parsetree = grammar.parse(variant)
                #if

                if ":c." in variant:
                    # Do the c2chrom dance
                    variant = converter.c2chrom(variant)

                attr["gName"] = variant

                if variant and ":g." in variant:
                    # Do the g2c dance
                    variants = converter.chrom2c(variant, "dict")
                    if variants:
                        out = ["%-10s:\t%s" % (key[:10], "\n\t\t".join(value))\
                               for key, value in variants.items()]
                        attr["cNames"].extend(out)
                    #if
                #if
            #if

            attr['messages'] = map(util.message_info, output.getMessages())
        return render.converter(attr)
#PositionConverter


class VariantInfo:
    """
    The I{g.} to I{c.} and vice versa interface for LOVD.

    Search for an NM number in the MySQL database, if the version number
    matches, get the start and end positions in a variant and translate these
    positions to I{g.} notation if the variant is in I{c.} notation and vice
    versa.
    - If no end position is present, the start position is assumed to be the
      end position.
    - If the version number is not found in the database, an error message is
      generated and a suggestion for an other version is given.
    - If the reference sequence is not found at all, an error is returned.
    - If no variant is present, the transcription start and end and CDS end
      in I{c.} notation is returned.
    - If the variant is not accepted by the nomenclature parser, a parse error
      will be printed.
    """
    def GET(self):
        """
        Get variant info and return the result as plain text.

        Parameters:
        - LOVD_ver: The version of the calling LOVD.
        - build: The human genome build (hg19 assumed).
        - acc: The accession number (NM number).
        - var: A description of the variant.

        Returns:
        - start_main   ; The main coordinate of the start position in I{c.}
                         (non-star) notation.
        - start_offset ; The offset coordinate of the start position in I{c.}
                         notation (intronic position).
        - end_main     ; The main coordinate of the end position in I{c.}
                         (non-star) notation.
        - end_offset   ; The offset coordinate of the end position in I{c.}
                         notation (intronic position).
        - start_g      ; The I{g.} notation of the start position.
        - end_g        ; The I{g.} notation of the end position.
        - type         ; The mutation type.

        Returns (alternative):
        - trans_start  ; Transcription start in I{c.} notation.
        - trans_stop   ; Transcription stop in I{c.} notation.
        - CDS_stop     ; CDS stop in I{c.} notation.
        """
        i = web.input(var='')
        LOVD_ver = i.LOVD_ver
        build = i.build
        acc = i.acc
        var = i.var

        output = Output(__file__, config.Output)

        output.addMessage(__file__, -1, 'INFO',
                          'Received %s:%s (LOVD_ver %s, build %s)' \
                          % (acc, var, LOVD_ver, build))

        Converter = Mapper.Converter(build, config, output)

        result = ''

        # If no variant is given, return transcription start, transcription
        # end and CDS stop in c. notation.
        if var:
            ret = Converter.mainMapping(acc, var)
        else:
            ret = Converter.giveInfo(acc)
            if ret:
                result = '%i\n%i\n%i' % ret

        if not result and not getattr(ret, 'startmain', None):
            out = output.getOutput('LOVDERR')
            if out:
                result = out[0]
            else:
                result = 'Unknown error occured'

        output.addMessage(__file__, -1, 'INFO',
                          'Finished processing %s:%s (LOVD_ver %s, build %s)' \
                          % (acc, var, LOVD_ver, build))

        result = '%i\n%i\n%i\n%i\n%i\n%i\n%s' \
                 % (ret.startmain, ret.startoffset, ret.endmain,
                    ret.endoffset, ret.start_g, ret.end_g, ret.mutationType)

        web.header('Content-Type', 'text/plain')

        if LOVD_ver == "2.0-23" : # Obsoleted error messages, remove when possible.
            return re.sub("^Error \(.*\):", "Error:", result)

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
        2. By having a 'variant' value in the cookie. The value is removed.

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
            # Run checker if cookie variant is not None
            variant = web.cookies().get('variant')
            web.setcookie('variant', '', 60)
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
        output = Output(__file__, config.Output)

        args = {
            'lastpost' : name
        }

        if not name:
            return render.check(args, standalone=not interactive)

        output.addMessage(__file__, -1, 'INFO', 'Received variant %s' % name)
        # Todo: The following is probably a problem elsewhere too.
        # We stringify the variant, because a unicode string crashes
        # Bio.Seq.reverse_complement in Mapper.py:607.
        variantchecker.check_variant(str(name), config, output)
        output.addMessage(__file__, -1, 'INFO',
                          'Finished processing variant %s' % name)

        errors, warnings, summary = output.Summary()
        record_type = output.getIndexedOutput('recordType', 0, '')
        reference = output.getIndexedOutput('reference', 0, '')

        if reference:
            if record_type == 'LRG':
                reference = reference + '.xml'
            else:
                reference = reference + '.gb'

        # This is a tuple (variant, position)
        parse_error = output.getOutput('parseError')
        if parse_error:
            parse_error[0] = parse_error[0].replace('<', '&lt;')

        genomic_dna = output.getIndexedOutput('molType', 0) != 'n'

        genomic_description = output.getIndexedOutput('genomicDescription', 0, '')

        # Todo: Generate the fancy HTML views for the proteins here instead
        # of in mutalyzer/variantchecker.py.
        args = {
            'lastpost'           : name,
            'messages'           : map(util.message_info, output.getMessages()),
            'summary'            : summary,
            'parseError'         : parse_error,
            'errors'             : errors,
            'genomicDescription' : (genomic_description, urllib.quote(genomic_description)),
            'chromDescription'   : output.getIndexedOutput('genomicChromDescription', 0),
            'genomicDNA'         : genomic_dna,
            'visualisation'      : output.getOutput('visualisation'),
            'descriptions'       : map(lambda d: (d, urllib.quote(d)), output.getOutput('descriptions')),
            'protDescriptions'   : output.getOutput('protDescriptions'),
            'oldProtein'         : output.getOutput('oldProteinFancy'),
            'altStart'           : output.getIndexedOutput('altStart', 0),
            'altProtein'         : output.getOutput('altProteinFancy'),
            'newProtein'         : output.getOutput('newProteinFancy'),
            'transcriptInfo'     : output.getIndexedOutput('hasTranscriptInfo', 0, False),
            'transcriptCoding'   : output.getIndexedOutput('transcriptCoding', 0, False),
            'exonInfo'           : output.getOutput('exonInfo'),
            'cdsStart_g'         : output.getIndexedOutput('cdsStart_g', 0),
            'cdsStart_c'         : output.getIndexedOutput('cdsStart_c', 0),
            'cdsStop_g'          : output.getIndexedOutput('cdsStop_g', 0),
            'cdsStop_c'          : output.getIndexedOutput('cdsStop_c', 0),
            'restrictionSites'   : output.getOutput('restrictionSites'),
            'legends'            : output.getOutput('legends'),
            'reference'          : reference
        }

        return render.check(args, standalone=not interactive)
#Check


class CheckForward:
    """
    Set the given variant in the cookie and redirect to the name checker.

    @todo: Cleaner solution (one without using a cookie).
    """
    def GET(self):
        """
        Set the 'variant' cookie value to the given variant and redirect
        to the name checker (where we will arrive by a GET request).

        Parameters:
        - mutationName: Variant to set in the cookie.
        """
        i = web.input(mutationName=None)
        web.setcookie('variant', i.mutationName, 5 * 60)  # Five minutes
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
        O = Output(__file__, config.Output)

        attr = {"percentage": 0}

        i = web.input(ajax=None)
        try:
            jobID = int(i.jobID)
            total = int(i.totalJobs)
        except ValueError:
            return
        D = Db.Batch(config.Db)
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
                     'SyntaxChecker', 'PositionConverter', or 'SnpConverter'.
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
                          (default), 'SyntaxChecker', 'PositionConverter', or
                          'SnpConverter'.
        """
        O = Output(__file__, config.Output)

        maxUploadSize = config.Batch.batchInputMaxSize

        attr = {"messages"      : [],
                "errors"        : [],
                "debug"         : [],
                "maxSize"       : float(maxUploadSize) / 1048576,
                "batchTypes"    : ["NameChecker",
                                   "SyntaxChecker",
                                   "PositionConverter",
                                   "SnpConverter"],
                "hideTypes"     : batchType and 'none' or '',
                "selected"      : "0",
                "batchType"     : batchType or "",
                "avail_builds"  : config.Db.dbNames[::-1],
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

            # Todo: These error messages could be delivered trough a template
            if not 'CONTENT_LENGTH' in web.ctx.environ.keys():
                web.header('Content-Type', 'text/plain')
                web.ctx.status = '411 Length required'
                return 'Content length required'
            if int(web.ctx.environ.get('CONTENT_LENGTH')) > maxUploadSize:
                web.header('Content-Type', 'text/plain')
                web.ctx.status = '413 Request entity too large'
                return 'Sorry, only files up to %s megabytes are accepted.' % (float(maxUploadSize) / 1048576)

            D = Db.Batch(config.Db)
            S = Scheduler.Scheduler(config.Scheduler, D)
            FileInstance = File.File(config.File, O)

            # Generate the fromhost URL from which the results can be fetched
            fromHost = web.ctx.homedomain + web.ctx.homepath + '/'
            #fromHost = "http://%s%s" % (
            #    req.hostname, req.uri.rsplit("/", 1)[0]+"/")

            job, columns = FileInstance.parseBatchFile(inFile.file)
            if job is None:
                O.addMessage(__file__, 4, "PRSERR", "Could not parse input"
                             " file, please check your file format.")
            else:
                #TODO: Add Binair Switches to toggle some events
                attr["jobID"] = S.addJob("BINSWITHCES", email, job, columns,
                                         fromHost, batchType, arg1)
                attr["totalJobs"] = len(job) or 1
                attr["messages"].append("Your file has been parsed and the job"
                                        " is scheduled, you will receive an email when the job is "
                                        "finished.")

            attr["errors"].extend(map(util.message_info, O.getMessages()))

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
        handle = open(os.path.join(config.Scheduler.resultsDir, file))
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
        maxUploadSize = config.Retriever.maxDldSize
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
        maxUploadSize = config.Retriever.maxDldSize

        O = Output(__file__, config.Output)
        D = Db.Cache(config.Db)
        R = Retriever.GenBankRetriever(config.Retriever, O, D)

        UD, errors = "", []

        i = web.input(invoermethode='', bestandsveld={}, urlveld='',
                      genesymbol='', organism='', fiveutr='', threeutr='',
                      chracc='', start='', stop='', orientation='')

        try:
            if i.invoermethode == "file" :
                if not 'CONTENT_LENGTH' in web.ctx.environ.keys():
                    web.header('Content-Type', 'text/plain')
                    web.ctx.status = '411 Length required'
                    return 'Content length required.'
                #if
                if int(web.ctx.environ.get('CONTENT_LENGTH')) > maxUploadSize :
                    web.header('Content-Type', 'text/plain')
                    web.ctx.status = '413 Request entity too large'
                    return 'Upload limit exceeded.'
                #if
                # Non-conforming clients (read: LOVD) might send the form
                # request urlencoded (and not as the requested multipart/
                # form-data). We try to support this anyway.
                if web.ctx.env.get('CONTENT_TYPE', '') \
                       == 'application/x-www-form-urlencoded' \
                       and isinstance(i.bestandsveld, str):
                    UD = R.uploadrecord(i.bestandsveld)
                elif not i.bestandsveld == None and i.bestandsveld.file:
                    # Todo: actually we should check if .file exists
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
                errors.extend(map(lambda m: str(m), O.getMessages()))

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
        xsl_handle = open(os.path.join(mutalyzer.package_root(), WSDL_VIEWER),
                          'r')
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


# WSGI application
application = app.wsgifunc()
