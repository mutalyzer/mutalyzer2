"""
General Mutalyzer website interface.
"""


SERVICE_SOAP_LOCATION = '/services'
SERVICE_JSON_LOCATION = '/json'
WSDL_VIEWER = 'templates/wsdl-viewer.xsl'
GENOME_BROWSER_URL = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={chromosome}:{start}-{stop}&hgt.customText={bed_file}'


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
from collections import defaultdict

from lxml import etree
from cStringIO import StringIO
from simpletal import simpleTALES
from simpletal import simpleTAL
from spyne.interface.wsdl import Wsdl11

import mutalyzer
from mutalyzer import util
from mutalyzer import config
from mutalyzer.grammar import Grammar
from mutalyzer.services import soap
from mutalyzer import variantchecker
from mutalyzer.output import Output
from mutalyzer.mapping import Converter
from mutalyzer import Db
from mutalyzer import Scheduler
from mutalyzer import Retriever
from mutalyzer import File
from mutalyzer import describe


# Show web.py debugging information.
web.config.debug = config.get('debug')


# URL dispatch table
urls = (
    '',                                         'RedirectHome',
    '/(index)?',                                'Static',
    '/(about)',                                 'Static',
    '/(help)',                                  'Static',
    '/(faq)',                                   'Static',
    '/(exercise)',                              'Static',
    '/(disclaimer)',                            'Static',
    '/(nameGenerator)',                         'Static',
    '/(webservices)',                           'Static',
    '/checkForward',                            'CheckForward',
    '/check',                                   'Check',
    '/bed',                                     'Bed',
    '/syntaxCheck',                             'SyntaxCheck',
    '/positionConverter',                       'PositionConverter',
    '/snp',                                     'Snp',
    '/descriptionExtract',                      'DescriptionExtractor',
    '/upload',                                  'Uploader',
    '/batch([a-zA-Z]+)?',                       'BatchChecker',
    '/progress',                                'BatchProgress',
    '/Results_(\d+)\.txt',                      'BatchResult',
    '/soap-api',                                'SoapApi',
    '/Variant_info',                            'VariantInfo',
    '/getGS',                                   'GetGS',
    '/download/([a-zA-Z-]+\.(?:py|cs|php|rb))', 'Download',
    '/downloads/([a-zA-Z\._-]+)',               'Downloads',
    '/Reference/([\da-zA-Z\._-]+)',             'Reference'
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
    #__init__

    def __getattr__(self, name):
        """
        Returns a template. Call the template to get a render.

        @arg name: Template name (usually a HTML filename without '.html').
        @return: Template render function.
        """
        filename = name

        def template(args={}, scheme='html', standalone=False,
                     prevent_caching=False, page=None):
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
            page = page or filename
            file = filename
            if scheme == 'html':
                file += '.html'
            path = os.path.join(self.path, file)

            context = simpleTALES.Context()

            context.addGlobal('interactive', not standalone)

            context.addGlobal('location', web.ctx.homedomain + web.ctx.homepath)

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
                # The following three lines are a hack to get class="active"
                # on the active menu item, working around TAL.
                active = defaultdict(lambda: 'menu')
                active[page] = 'menu active'
                context.addGlobal('active', active)
                templateFile = open(os.path.join(self.path, 'menu.html'), 'r')
                template = simpleTAL.compileHTMLTemplate(templateFile)
                templateFile.close()

            if scheme == 'html':
                web.header('Content-Type', 'text/html')

            if prevent_caching:
                web.header('Cache-Control', 'no-cache')
                web.header('Expires', '-1')

            io = StringIO()
            template.expand(context, io)

            return io.getvalue()
        #template

        return template
    #__getattr__
#render_tal


# TAL template render
render = render_tal(os.path.join(mutalyzer.package_root(), 'templates'),
    globals = {
    'version'             : mutalyzer.__version__,
    'nomenclatureVersion' : mutalyzer.NOMENCLATURE_VERSION,
    'releaseDate'         : mutalyzer.__date__,
    'release'             : mutalyzer.RELEASE,
    'copyrightYears'      : mutalyzer.COPYRIGHT_YEARS,
    'contactEmail'        : config.get('email'),
    'serviceSoapLocation' : SERVICE_SOAP_LOCATION,
    'serviceJsonLocation' : SERVICE_JSON_LOCATION,
    'piwik'               : config.get('piwik'),
    'piwikBase'           : config.get('piwikBase'),
    'piwikSite'           : config.get('piwikSite')
})

# web.py application
app = web.application(urls, globals(), autoreload = False)


class RedirectHome:
    """
    Permanent redirect to the homepage.
    """
    def GET(self):
        """
        Redirect to / and include the query string.
        """
        raise web.redirect('/' + web.ctx.query)
    #GET

    def POST(self):
        """
        Redirect to / and include the query string.
        """
        raise web.redirect('/' + web.ctx.query)
    #POST
#RedirectHome


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
        called with file = '../mutalyzer.conf' for example.

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
        return content.format(path = web.ctx.homedomain + web.ctx.homepath)
    #GET
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
        called with file = '../../mutalyzer.conf' for example.

        The url routing currently makes sure to only call this with filenames
        of the form [a-zA-Z\._-]+.
        """
        file_path = os.path.join(mutalyzer.package_root(), 'templates',
                                 'downloads', file)

        if not os.path.isfile(file_path):
            raise web.notfound()

        handle = open(file_path)
        F = File.File(None)
        web.header('Content-Type', F.getMimeType(handle)[0])
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)

        return handle.read()
    #GET
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
        called with file = '../../mutalyzer.conf' for example.

        The url routing currently makes sure to only call this with filenames
        of the form [a-zA-Z\._-]+.
        """
        file_path = os.path.join(config.get('cache'), '%s.bz2' % file)

        if not os.path.isfile(file_path):
            raise web.notfound()

        handle = bz2.BZ2File(file_path, 'r')
        web.header('Content-Type', 'text/plain')
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)

        return handle.read()
    #GET

    def HEAD(self, file):
        """
        Do the same as in the GET case, but don't actually bunzip and send the
        file, just check if it exists.

        @arg file: Filename to download from cache.
        @type file: string

        This is used by LOVD to quickly check if a reference file is in the
        cache. If it isn't, it will resubmit it.
        Of course a more proper solution here would be to have some web
        service method which checks if the GenBank file is in the cache *or*
        can be reconstructed from the information in the database. Because if
        the latter is the case, Mutalyzer will add it to the cache on the fly.
        """
        file_path = os.path.join(config.get('cache'), '%s.bz2' % file)

        if not os.path.isfile(file_path):
            # The following is a hack to return a 404 not found status with
            # empty body (as is checked by our unit test framework, WebTest).
            # Just passing nothing, or the empty string, causes web.py to
            # insert some default 'not found' message.
            class TrueEmptyString(object):
                def __str__(self):
                    return ''
                def __nonzero__( self):
                    return True
            #TrueEmptyString

            raise web.notfound(message = TrueEmptyString())
        #if

        web.header('Content-Type', 'text/plain')
        web.header('Content-Disposition', 'attachment; filename="%s"' % file)

        return ''
    #HEAD
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
        output = Output(__file__)
        IP = web.ctx["ip"]

        i = web.input(mutationName=None, variantRecord=None, forward=None)

        output.addMessage(__file__, -1, 'INFO',
            'Received request getGS(%s, %s, %s) from %s' % (i.mutationName,
            i.variantRecord, i.forward, IP))

        # Todo: The following is probably a problem elsewhere too.
        # We stringify the variant, because a unicode string crashes
        # Bio.Seq.reverse_complement in mapping.py:607.

        variantchecker.check_variant(str(i.mutationName), output)

        output.addMessage(__file__, -1, 'INFO',
            'Finished request getGS(%s, %s, %s)' % (i.mutationName,
            i.variantRecord, i.forward))

        legends = output.getOutput("legends")

        # Filter the transcript from the legend
        legends = [l for l in legends if "_v" in l[0]]
        for l in legends:
            if l[1] == i.variantRecord:
                if i.forward:
                    p, a = i.mutationName.split(':')
                    return Check.check(p+'('+l[0]+'):'+a, standalone=True)
                else:
                    web.header('Content-Type', 'text/plain')
                    return l[0]

        web.header('Content-Type', 'text/plain')

        return "Transcript not found"#+`legends`
    #GET
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
    #GET

    def POST(self):
        """
        Parse the given variant and render the syntax checker HTML form.

        Parameters:
        - variant: Variant name to check.
        """
        output = Output(__file__)
        IP = web.ctx["ip"]
        i = web.input()

        output.addMessage(__file__, -1, 'INFO',
            'Received request syntaxCheck(%s) from %s' % (i.variant, IP))

        variant = i.variant
        if variant.find(',') >= 0:
            output.addMessage(__file__, 2, "WCOMMASYNTAX",
                "Comma's are not allowed in the syntax, autofixed.")
            variant = variant.replace(',', '')
            #args["variant"]=variant

        grammar = Grammar(output)
        grammar.parse(variant)

        pe = output.getOutput("parseError")
        if pe:
            pe[0] = pe[0].replace('<', "&lt;")

        args = {
            "variant"       : variant,
            "messages"      : map(util.message_info, output.getMessages()),
            "parseError"    : pe,
            "debug"         : ""
        }

        output.addMessage(__file__, -1, 'INFO',
            'Finished request syntaxCheck(%s)' % i.variant)

        return render.parse(args)
    #POST
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
    #GET

    def POST(self):
        """
        Convert to HGVS description(s) and render SNP converter HTML form.

        Parameters:
          - rsId: The dbSNP rs number.
        """
        i = web.input(rsId=None)
        return self.snp(i.rsId)
    #POST

    def snp(self, rs_id=None):
        """
        Convert {rs_id} to HGVS description(s) and render SNP converter HTML
        form.

        @kwarg rs_id: The dbSNP rs number (including 'rs' prefix).
        @type rs_id: string
        """
        output = Output(__file__)

        IP = web.ctx["ip"]

        descriptions = []

        if rs_id:
            output.addMessage(__file__, -1, 'INFO',
                'Received request snpConvert(%s) from %s' % (rs_id, IP))
            retriever = Retriever.Retriever(output, None)
            descriptions = retriever.snpConvert(rs_id)
            output.addMessage(__file__, -1, 'INFO',
                'Finished request snpConvert(%s)' % rs_id)

        args = {
            'snp'      : descriptions,
            'messages' : map(util.message_info, output.getMessages()),
            'summary'  : output.Summary()[2],
            'lastpost' : rs_id
        }

        return render.snp(args)
    #snp
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
    #GET

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
        # Bio.Seq.reverse_complement in mapping.py:607.
        return self.position_converter(i.build, str(i.variant))
    #POST

    def position_converter(self, build='', variant=''):
        """
        Convert a variant and render position converter HTML form.

        @kwarg build: Human genome build (currently 'hg18' or 'hg19').
        @kwarg variant: Variant to convert.
        """
        output = Output(__file__)
        IP = web.ctx["ip"]

        avail_builds = config.get('dbNames')[::-1]

        if build :
            avail_builds.remove(build)
            avail_builds.insert(0, build)

        attr = {
            "avail_builds" : avail_builds,
            "variant"      : variant,
            "gName"        : "",
            "cNames"       : [],
            "messages"     : [],
            "posted"       : build and variant
        }

        if build and variant:

            output.addMessage(__file__, -1, 'INFO',
                'Received request positionConverter(%s, %s) from %s' % (
                build, variant, IP))

            # Todo: check for correct build.
            converter = Converter(build, output)

            #Convert chr accNo to NC number
            variant = converter.correctChrVariant(variant)

            if variant:
                if not(":c." in variant or ":n." in variant or ":g." in variant):
                    #Bad name
                    grammar = Grammar(output)
                    grammar.parse(variant)

                if ":c." in variant or ":n." in variant:
                    # Do the c2chrom dance
                    variant = converter.c2chrom(variant)

                attr["gName"] = variant

                if variant and ":g." in variant:
                    # Do the g2c dance
                    variants = converter.chrom2c(variant, "dict")
                    if variants is None:
                        attr['gName'] = None
                    elif variants:
                        out = ["%-10s:\t%s" % (key[:10], "\n\t\t".join(value))
                               for key, value in variants.items()]
                        attr["cNames"].extend(out)

            attr['messages'] = map(util.message_info, output.getMessages())

            output.addMessage(__file__, -1, 'INFO',
                'Finished request positionConverter(%s, %s)' % (build,
                variant))

        return render.converter(attr)
    #position_converter
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
        i = web.input(var = '')
        LOVD_ver = i.LOVD_ver
        build = i.build
        acc = i.acc
        var = i.var

        output = Output(__file__)
        IP = web.ctx["ip"]

        output.addMessage(__file__, -1, 'INFO',
            'Received request variantInfo(%s:%s (LOVD_ver %s, build %s))'
            ' from %s' % (acc, var, LOVD_ver, build, IP))

        converter = Converter(build, output)

        result = ''

        # If no variant is given, return transcription start, transcription
        # end and CDS stop in c. notation.
        if var:
            ret = converter.mainMapping(acc, var)
        else:
            ret = converter.giveInfo(acc)
            if ret:
                result = '%i\n%i\n%i' % ret

        if not result and not getattr(ret, 'startmain', None):
            out = output.getOutput('LOVDERR')
            if out:
                result = out[0]
            else:
                result = 'Unknown error occured'

        output.addMessage(__file__, -1, 'INFO',
            'Finished request variantInfo(%s:%s (LOVD_ver %s, build %s))' % (
            acc, var, LOVD_ver, build))

        if not result and getattr(ret, 'startmain', None):
            result = '%i\n%i\n%i\n%i\n%i\n%i\n%s' % (ret.startmain,
            ret.startoffset, ret.endmain, ret.endoffset, ret.start_g,
            ret.end_g, ret.mutationType)

        web.header('Content-Type', 'text/plain')

        if LOVD_ver == "2.0-23": # Obsoleted error messages, remove soon.
            return re.sub("^Error \(.*\):", "Error:", result)

        return result
    #GET
#VariantInfo


class CheckForward:
    """
    Old entrypoint to the namechecker. We keep it to not break existing
    bookmarks (but this could also be done with an Apache rewrite rule).
    """
    def GET(self):
        """
        Permanently redirect to the name checker.

        Parameters:
        - mutationName: Variant to check.
        """
        i = web.input(mutationName=None)
        raise web.redirect('/check?name=' + urllib.quote(i.mutationName))
    #GET
#CheckForward


class Check:
    """
    The variant checker.
    """
    def GET(self):
        """
        Render the variant checker HTML form.

        For backwards compatibility with older LOVD versions, we support the
        'mutationName' argument. If present, we redirect and add standalone=1.

        Parameters:
        - name: Variant to check.
        """
        i = web.input(name=None, mutationName=None, standalone=False)

        if i.mutationName:
            raise web.redirect('/check?name=%s&standalone=1'
                               % urllib.quote(i.mutationName))

        return self.check(i.name, standalone=bool(i.standalone))
    #GET

    def POST(self):
        """
        For now we also accept POST requests with a permanent redirect.
        """
        i = web.input(name=None, mutationName=None, standalone=False)
        raise web.redirect('/check?name=%s'
                           % urllib.quote(i.name or i.mutationName))
    #POST

    @staticmethod
    def check(name=None, standalone=False):
        """
        Render the variant checker HTML form. If the name argument is given,
        run the name checker.

        @kwarg name: Variant to check.
        @kwarg interactive: Run interactively, meaning we wrap the result in
            the site layout and include the HTML form.
        """
        if not name:
            return render.check(dict(name=None), standalone=standalone)

        output = Output(__file__)
        output.addMessage(__file__, -1, 'INFO', 'Received variant %s from %s'
                          % (name, web.ctx['ip']))

        # Todo: The following is probably a problem elsewhere too.
        # We stringify the variant, because a unicode string crashes
        # Bio.Seq.reverse_complement in mapping.py:607.
        variantchecker.check_variant(str(name), output)

        errors, warnings, summary = output.Summary()
        record_type = output.getIndexedOutput('recordType', 0, '')
        reference = output.getIndexedOutput('reference', 0, '')

        if reference:
            if record_type == 'LRG':
                reference = reference + '.xml'
            else :
                reference = reference + '.gb'

        # This is a tuple (variant, position) if we have a parse error
        parse_error = output.getOutput('parseError')
        if parse_error:
            parse_error[0] = parse_error[0].replace('<', '&lt;')

        genomic_dna = output.getIndexedOutput('molType', 0) != 'n'
        genomic_description = output.getIndexedOutput('genomicDescription', 0, '')

        # Create a tuple (description, link) from a description
        def description_to_link(description):
            link = None
            if description[-1] != '?':
                link = urllib.quote(description)
            return description, link

        # Create a link to the UCSC Genome Browser
        browser_link = None
        raw_variants = output.getIndexedOutput('rawVariantsChromosomal', 0)
        if raw_variants:
            positions = [pos
                for descr, (first, last) in raw_variants[2]
                for pos in (first, last)]
            bed_url = web.ctx.homedomain + web.ctx.homepath + \
                '/bed?name=' + urllib.quote(name)
            browser_link = GENOME_BROWSER_URL.format(
                chromosome=raw_variants[0], start=min(positions) - 10,
                stop=max(positions) + 10, bed_file=urllib.quote(bed_url))

        allele = describe.describe(output.getIndexedOutput("original", 0),
                                   output.getIndexedOutput("mutated", 0))
        prot_allele = describe.describe(output.getIndexedOutput("oldprotein", 0),
                                        output.getIndexedOutput("newprotein", 0, default=""), DNA=False)

        extracted = extractedProt = '(skipped)'

        if allele:
            extracted = describe.alleleDescription(allele)
        if prot_allele:
            extractedProt = describe.alleleDescription(prot_allele)


        # Todo: Generate the fancy HTML views for the proteins here instead
        # of in mutalyzer/variantchecker.py.
        args = {
            'name'               : name,
            'messages'           : map(util.message_info, output.getMessages()),
            'summary'            : summary,
            'parseError'         : parse_error,
            'errors'             : errors,
            'genomicDescription' : (genomic_description, urllib.quote(genomic_description)),
            'chromDescription'   : output.getIndexedOutput('genomicChromDescription', 0),
            'genomicDNA'         : genomic_dna,
            'visualisation'      : output.getOutput('visualisation'),
            'descriptions'       : map(description_to_link, output.getOutput('descriptions')),
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
            'reference'          : reference,
            'browserLink'        : browser_link,
            'extractedDescription' : (extracted, urllib.quote(extracted)),
            'extractedProtein'   : (extractedProt, urllib.quote(extractedProt))
        }

        output.addMessage(__file__, -1, 'INFO', 'Finished variant %s' % name)

        return render.check(args, standalone=standalone, prevent_caching=True)
    #check
#Check


class DescriptionExtractor:
    """
    The Variant Description Extractor.
    """
    def GET(self):
        """
        Render the description extractor HTML form.
        """
        return self.descriptionExtract()
    #GET

    def POST(self):
        """
        Run the description extractor and render the description extractor HTML
        form.

        Parameters:
        - referenceSeq:
        - variantSeq:
        """
        i = web.input(referenceSeq=None, variantSeq=None)
        return self.descriptionExtract(i.referenceSeq, i.variantSeq)
    #POST

    @staticmethod
    def descriptionExtract(referenceSeq=None, variantSeq=None):
        """
        Render the description extractor HTML form. If the referenceSeq and
        variantSeq argument are given, run the description extractor.

        @kwarg referenceSeq: The reference sequence.
        @type referenceSeq: string
        @kwarg variantSeq: The observed sequence.
        @type variantSeq: string
        """
        output = Output(__file__)
        IP = web.ctx["ip"]

        args = {
            'lastReferenceSeq' : referenceSeq,
            'lastVariantSeq'   : variantSeq
        }

        if not (referenceSeq and variantSeq):
            return render.descriptionExtract(args)

        output.addMessage(__file__, -1, 'INFO',
            "Received Description Extract request from %s" % IP)

        # Move this to the describe module.
        if not util.is_dna(referenceSeq):
            output.addMessage(__file__, 3, "ENODNA",
                "Reference sequence is not DNA.")
        if not util.is_dna(variantSeq):
            output.addMessage(__file__, 3, "ENODNA",
                "Variant sequence is not DNA.")

        result = describe.describe(referenceSeq, variantSeq)
        description = describe.alleleDescription(result)

        errors, warnings, summary = output.Summary()

        visualisation = []
        for i in result:
            visualisation.append([i.start, i.end, i.type, i.deleted,
                i.inserted, i.shift, i.hgvs])

        args = {
            'lastReferenceSeq' : referenceSeq,
            'lastVariantSeq'   : variantSeq,
            'description'      : description,
            'visualisation'    : visualisation,
            'errors'           : errors,
            'summary'          : summary,
            'messages'         : map(util.message_info,
                output.getMessages())
        }

        output.addMessage(__file__, -1, 'INFO',
            "Finished Description Extract request")

        return render.descriptionExtract(args)
    #descriptionExtract
#DescriptionExtract


class Bed:
    """
    Create BED track.
    """
    def GET(self):
        """
        Create a BED track for the given variant, listing the positioning of
        its raw variants. E.g. for use in the UCSC Genome Browser.

        Parameters:
        - name: Variant to create BED track for.

        This basically just runs the variant checker and extracts the raw
        variants with positions.
        """
        web.header('Content-Type', 'text/plain')

        i = web.input(name=None)
        name = i.name

        if not name:
            web.ctx.status = '404 Not Found'
            return 'Sorry, we have not BED track for this variant.'

        output = Output(__file__)

        variantchecker.check_variant(str(name), output)

        raw_variants = output.getIndexedOutput('rawVariantsChromosomal', 0)
        if not raw_variants:
            web.ctx.status = '404 Not Found'
            return 'Sorry, we have no BED track for this variant.'

        fields = {
            'name'       : 'Mutalyzer',
            'description': 'Mutalyzer track for ' + name,
            'visibility' : 'pack',
            'db'         : 'hg19',
            'url'        : web.ctx.homedomain + web.ctx.homepath +
                '/check?name=' + urllib.quote(name),
            'color':       '255,0,0'}

        bed = ' '.join(['track'] + [
            '%s="%s"' % field for field in fields.items()]) + '\n'

        for description, positions in raw_variants[2]:
            bed += '\t'.join([raw_variants[0], str(min(positions) - 1),
                str(max(positions)), description, '0', raw_variants[1]]) + '\n'

        return bed
    #GET
#Bed


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
        @todo: Actually, signaling 'OK' here only means the last entry was
            taken from the database queue. It might still be processing, in
            which case not all output is yet written to the result file.
            For the standard use case, this is no big deal, since any user
            will take more than a few milliseconds to actually click the
            download link.
            However, if we imagine some scripted batch uploader, it might get
            bitten by this bug. (This includes our unit tests, where we work
            around it by explicitely waiting a second.)
        """
        attr = {"percentage": 0}

        i = web.input(ajax = None)
        try:
            jobID = int(i.jobID)
            total = int(i.totalJobs)
        except ValueError:
            return

        D = Db.Batch()
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

        return render.progress(attr)
    #GET
#BatchProgress


class BatchChecker:
    """
    Run batch jobs.
    """
    def GET(self, batchType = None):
        """
        Render batch checker HTML form.

        @kwarg batchType: Type of batch job.
        """
        return self.batch(batchType=batchType)
    #GET

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
        i = web.input(batchEmail=None, batchFile={}, arg1='', batchType=None)

        return self.batch(email=i.batchEmail, inFile=i.batchFile, arg1=i.arg1,
                          batchType=i.batchType)
    #POST

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
        O = Output(__file__)

        maxUploadSize = config.get('batchInputMaxSize')

        attr = {
            "messages"      : [],
            "errors"        : [],
            "debug"         : [],
            "maxSize"       : float(maxUploadSize) / 1048576,
            "batchTypes"    : ["NameChecker", "SyntaxChecker",
                               "PositionConverter", "SnpConverter"],
            "hideTypes"     : batchType and 'none' or '',
            "selected"      : "0",
            "batchType"     : batchType or "",
            "avail_builds"  : config.get('dbNames')[::-1],
            "jobID"         : None,
            "totalJobs"     : None
        }

        # Make sure the correct page is displayed for an entrypoint
        if batchType:
            page = 'batch' + batchType
        else:
            page = 'batch'
            batchType = 'NameChecker'

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
                return 'Sorry, only files up to %s megabytes are accepted.' % (
                    float(maxUploadSize) / 1048576)

            D = Db.Batch()
            S = Scheduler.Scheduler(D)
            FileInstance = File.File(O)

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
                    " is scheduled, you will receive an email when the job is"
                    " finished.")

            attr["errors"].extend(map(util.message_info, O.getMessages()))

        return render.batch(attr, page=page)
    #batch
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
        called with result = '../../mutalyzer.conf' for example.

        The url routing currently makes sure to only call this with filenames
        of the form \d+.
        """
        filename = 'Results_%s.txt' % result
        handle = open(os.path.join(config.get('resultsDir'), filename))
        web.header('Content-Type', 'text/plain')
        web.header('Content-Disposition',
            'attachment; filename="%s"' % filename)

        return handle.read()
    #GET
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
    inpv = inpv.replace(',', '').replace('.', '').replace('-', '')

    try:
        return int(inpv)
    except ValueError:
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
        maxUploadSize = config.get('maxDldSize')
        available_assemblies = config.get('dbNames')[::-1]
        UD, errors = "", []
        args = {
            'UD'                   : UD,
            'available_assemblies' : available_assemblies,
            'maxSize'              : float(maxUploadSize) / 1048576,
            'errors'               : errors
        }
        return render.gbupload(args)
    #GET

    def POST(self):
        """
        Render reference sequence uploader form and handle a reference
        sequence upload or retrieval.

        This handler has four methods:
        1. The reference sequence file is a local file.
        2. The reference sequence file can be found at the following URL.
        3. Retrieve part of the reference genome for a (HGNC) gene symbol.
        4. Retrieve a range of a chromosome by accession number.
        5. Retrieve a range of a chromosome by name.

        Parameters:
        - invoermethode: Input method. One of 'file', 'url', 'gene', 'chr',
          'chrname'.

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

        Parameters (method 'chrname'):
        - chrnameassembly: Genome assembly (probably 'hg18' or 'hg19').
        - chrname: Chromosome name.
        - chrnamestart: Start position.
        - chrnamestop: Stop position.
        - chrnameorientation: Orientation.
        """
        maxUploadSize = config.get('maxDldSize')
        available_assemblies = config.get('dbNames')[::-1]

        O = Output(__file__)
        IP = web.ctx["ip"]
        D = Db.Cache()
        R = Retriever.GenBankRetriever(O, D)

        UD, errors = "", []

        i = web.input(invoermethode='', bestandsveld={}, urlveld='',
                      genesymbol='', organism='', fiveutr='', threeutr='',
                      chracc='', start='', stop='', orientation='',
                      chrnameassembly='', chrname='', chrnamestart='',
                      chrnamestop='', chrnameorientation='')

        O.addMessage(__file__, -1, 'INFO',
            'Received request'
            ' upload(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) from %s' % (
            i.invoermethode, i.urlveld, i.genesymbol, i.organism, i.fiveutr,
            i.threeutr, i.chracc, i.start, i.stop, i.orientation, i.chrnameassembly,
            i.chrname, i.chrnamestart, i.chrnamestop, i.chrnameorientation, IP))

        try:
            if i.invoermethode == "file":
                if not 'CONTENT_LENGTH' in web.ctx.environ.keys():
                    web.header('Content-Type', 'text/plain')
                    web.ctx.status = '411 Length required'
                    return 'Content length required.'

                if int(web.ctx.environ.get('CONTENT_LENGTH')) > maxUploadSize:
                    web.header('Content-Type', 'text/plain')
                    web.ctx.status = '413 Request entity too large'
                    return 'Upload limit exceeded.'

                # Non-conforming clients (read: LOVD) might send the form
                # request urlencoded (and not as the requested multipart/
                # form-data). We try to support this anyway.
                if web.ctx.env.get('CONTENT_TYPE', '') == \
                    'application/x-www-form-urlencoded' and \
                    isinstance(i.bestandsveld, str):
                    UD = R.uploadrecord(i.bestandsveld)
                elif not i.bestandsveld == None and i.bestandsveld.file:
                    # Todo: actually we should check if .file exists
                    UD = R.uploadrecord(i.bestandsveld.file.read())
                else:
                    raise web.badrequest()

            elif i.invoermethode == "url":
                UD = R.downloadrecord(i.urlveld)
            elif i.invoermethode == "gene":
                geneName = i.genesymbol
                organism = i.organism
                upStream = _checkInt(i.fiveutr, "5' flanking nucleotides")
                downStream = _checkInt(i.threeutr, "3' flanking nucleotides")
                UD = R.retrievegene(geneName, organism, upStream, downStream)
            elif i.invoermethode == "chr":
                accNo = i.chracc
                start = _checkInt(i.start, "Start position")
                stop = _checkInt(i.stop, "Stop position")
                orientation = int(i.orientation)
                UD = R.retrieveslice(accNo, start, stop, orientation)
            elif i.invoermethode == "chrname":
                build = i.chrnameassembly
                name = i.chrname
                start = _checkInt(i.chrnamestart, "Start position")
                stop = _checkInt(i.chrnamestop, "Stop position")
                orientation = int(i.chrnameorientation)

                if build not in available_assemblies:
                    raise InputException('Assembly not available: %s' % build)

                if not name.startswith('chr'):
                    name = 'chr%s' % name

                database = Db.Mapping(build)
                accession = database.chromAcc(name)

                if not accession:
                    raise InputException('Chromosome not available for build %s: %s' %
                                         (build, name))

                UD = R.retrieveslice(accession, start, stop, orientation)
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
            "UD"                   : UD,
            'available_assemblies' : available_assemblies,
            "maxSize"              : float(maxUploadSize) / 1048576,
            "errors"               : errors
        }

        O.addMessage(__file__, -1, 'INFO',
            'Finished request upload(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)' \
            % (i.invoermethode, i.urlveld, i.genesymbol, i.organism,
               i.fiveutr, i.threeutr, i.chracc, i.start, i.stop, i.orientation,
               i.chrnameassembly, i.chrname, i.chrnamestart, i.chrnamestop, i.chrnameorientation))

        return render.gbupload(args)
    #POST
#Uploader


class SoapApi:
    """
    SOAP web service documentation.
    """
    def GET(self):
        """
        HTML documentation for the web service.

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
               web service.
        @todo: Use configuration value for .xsl location.
        @todo: Cache this transformation.
        """
        url = web.ctx.homedomain + web.ctx.homepath + SERVICE_SOAP_LOCATION
        wsdl = Wsdl11(soap.application.interface)
        wsdl.build_interface_document(url)
        wsdl_handle = StringIO(wsdl.get_interface_document())
        xsl_handle = open(os.path.join(mutalyzer.package_root(), WSDL_VIEWER),
                          'r')
        wsdl_doc = etree.parse(wsdl_handle)
        xsl_doc = etree.parse(xsl_handle)
        transform = etree.XSLT(xsl_doc)

        web.header('Content-Type', 'text/html')
        return str(transform(wsdl_doc))
    #GET
#SoapApi


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
    #GET
#Static
