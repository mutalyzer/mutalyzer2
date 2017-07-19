"""
Mutalyzer website views.
"""


from __future__ import unicode_literals

import bz2
import os
import pkg_resources
import re
import StringIO
import urllib

from datetime import datetime
from flask import Blueprint
from flask import (abort, jsonify, make_response, redirect, render_template,
                   request, send_from_directory, url_for)
import jinja2
from lxml import etree
from spyne.server.http import HttpBase
from sqlalchemy.orm.exc import NoResultFound

import extractor

import mutalyzer
from mutalyzer import (announce, backtranslator, File, ncbi, Retriever,
                       Scheduler, stats, util, variantchecker)
from mutalyzer.config import settings
from mutalyzer.db.models import BATCH_JOB_TYPES
from mutalyzer.db.models import Assembly, BatchJob
from mutalyzer.grammar import Grammar
from mutalyzer.mapping import Converter
from mutalyzer.output import Output
from mutalyzer.services import soap


website = Blueprint('website', __name__)


def global_context():
    """
    Create a context of global template variables.
    """
    # Note that this cannot be a static module variable, since we want to call
    # `announce.get_announcement()` on each request, not once during startup.
    return {
        'mutalyzer_version'   : mutalyzer.__version__,
        'nomenclature_version': mutalyzer.NOMENCLATURE_VERSION,
        'release_date'        : mutalyzer.__date__,
        'release'             : mutalyzer.__version_info__[-1] != 'dev',
        'copyright_years'     : mutalyzer.COPYRIGHT_YEARS,
        'contact_email'       : settings.EMAIL,
        'soap_wsdl_url'       : settings.SOAP_WSDL_URL or 'https://mutalyzer.nl/services/?wsdl',
        'json_root_url'       : settings.JSON_ROOT_URL or 'https://mutalyzer.nl/json',
        'piwik'               : settings.PIWIK,
        'piwik_base_url'      : settings.PIWIK_BASE_URL,
        'piwik_site_id'       : settings.PIWIK_SITE_ID,
        'announcement'        : announce.get_announcement()
    }


def request_terms():
    """
    List of terms associated with the request (i.e., from the request path and
    query string).
    """
    terms = request.path.lstrip('/').split('/')
    terms += [s for item in request.args.iteritems() for s in item]
    return terms


@website.app_template_filter('short')
def short_seq_filter(s, max_len=21):
    if len(s) > max_len:
        return s[:max_len / 2 - 1] + '...' + s[-max_len / 2 + 2:]
    return s


@website.context_processor
def add_globals():
    return global_context()


@website.errorhandler(404)
def error_not_found(error):
    return render_template('not-found.html', terms=request_terms()), 404


@website.app_errorhandler(404)
def app_error_not_found(error):
    return render_template('not-found.html', terms=request_terms(),
                           **global_context()), 404


@website.route('/')
def homepage():
    """
    Website homepage.
    """
    return render_template('homepage.html')


@website.route('/about')
def about():
    """
    About page.
    """
    return render_template('about.html', counter_totals=stats.get_totals())


@website.route('/name-generator')
def name_generator():
    """
    Name generator page.
    """
    return render_template('name-generator.html')


@website.route('/webservices')
def webservices():
    """
    Webservices documentation page.
    """
    return render_template('webservices.html')


@website.route('/soap-api')
def soap_api():
    """
    SOAP web service documentation by rendering the WSDL document as HTML.

    Generate the documentation by a XSL transform of the WSDL document.
    The XSL transformation used is `WSDL viewer
    <http://tomi.vanek.sk/index.php?page=wsdl-viewer>`_ by Tomi Vanek.

    We apply a small patch to this transformation to show newlines in
    the SOAP method docstrings:

    1. Around line 1195, the description ``<div>``, replace
       ``<div class="value">`` by ``<div class="value documentation">``.

    2. In the style sheet, add:

           .documentation { white-space: pre-line; }
    """
    soap_server = HttpBase(soap.application)
    soap_server.doc.wsdl11.build_interface_document(
        settings.SOAP_WSDL_URL or 'https://mutalyzer.nl/services/?wsdl')
    wsdl_string = soap_server.doc.wsdl11.get_interface_document()

    xsl_file = os.path.join(
        pkg_resources.resource_filename('mutalyzer', 'website/templates'),
        'wsdl-viewer.xsl')
    wsdl_doc = etree.fromstring(wsdl_string)
    xsl_doc = etree.parse(xsl_file)
    transform = etree.XSLT(xsl_doc)

    return make_response(unicode(transform(wsdl_doc)))


@website.route('/downloads/<string:filename>')
def downloads(filename):
    """
    Files for download, such as example webservice client scripts and batch
    job input files.
    """
    try:
        response = make_response(render_template(os.path.join('downloads',
                                                              filename)))
    except jinja2.exceptions.TemplateNotFound:
        abort(404)

    response.headers['Content-Type'] = 'text/plain; charset=utf-8'
    response.headers['Content-Disposition'] = ('attachment; filename="%s"'
                                               % filename)
    return response


@website.route('/syntax-checker')
def syntax_checker():
    """
    Parse the given variant and render the syntax checker HTML form.
    """
    # Backwards compatibility.
    if 'variant' in request.args:
        return redirect(url_for('.syntax_checker',
                                description=request.args['variant']),
                        code=301)

    description = request.args.get('description')

    if not description:
        return render_template('syntax-checker.html')

    output = Output(__file__)
    output.addMessage(__file__, -1, 'INFO',
                      'Received request syntaxCheck(%s) from %s'
                      % (description, request.remote_addr))
    stats.increment_counter('syntax-checker/website')

    grammar = Grammar(output)
    grammar.parse(description)

    parse_error = output.getOutput('parseError')
    messages = map(util.message_info, output.getMessages())

    output.addMessage(__file__, -1, 'INFO',
                      'Finished request syntaxCheck(%s)' % description)

    return render_template('syntax-checker.html',
                           description=description,
                           messages=messages,
                           parse_error=parse_error)


@website.route('/name-checker')
def name_checker():
    """
    Name checker.
    """
    # For backwards compatibility with older LOVD versions, we support the
    # `mutationName` argument. If present, we redirect and add `standalone=1`.
    #
    # Also for backwards compatibility, we support the `name` argument as an
    # alias for `description`.
    if 'name' in request.args:
        return redirect(url_for('.name_checker',
                                description=request.args['name'],
                                standalone=request.args.get('standalone')),
                        code=301)
    if 'mutationName' in request.args:
        return redirect(url_for('.name_checker',
                                description=request.args['mutationName'],
                                standalone=1),
                        code=301)

    description = request.args.get('description')

    if not description:
        return render_template('name-checker.html')

    output = Output(__file__)
    output.addMessage(__file__, -1, 'INFO', 'Received variant %s from %s'
                      % (description, request.remote_addr))
    stats.increment_counter('name-checker/website')

    variantchecker.check_variant(description, output)

    errors, warnings, summary = output.Summary()
    parse_error = output.getOutput('parseError')

    record_type = output.getIndexedOutput('recordType', 0, '')
    reference = output.getIndexedOutput('reference', 0, '')
    if reference:
        if record_type == 'LRG':
            reference_filename = reference + '.xml'
        elif record_type == 'GB':
            reference_filename = reference + '.gb'
        else:
            reference_filename = None
    else:
        reference_filename = None

    genomic_dna = output.getIndexedOutput('molType', 0) != 'n'
    genomic_description = output.getIndexedOutput('genomicDescription', 0, '')

    # Create a link to the UCSC Genome Browser.
    browser_link = None
    raw_variants = output.getIndexedOutput('rawVariantsChromosomal', 0)
    if raw_variants:
        positions = [pos
                     for descr, (first, last) in raw_variants[2]
                     for pos in (first, last)]
        bed_url = url_for('.bed', description=description, _external=True)
        browser_link = ('http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&'
                        'position={chromosome}:{start}-{stop}&hgt.customText='
                        '{bed_file}'.format(chromosome=raw_variants[0],
                                            start=min(positions) - 10,
                                            stop=max(positions) + 10,
                                            bed_file=urllib.quote(bed_url)))

    # Experimental description extractor.
    if (output.getIndexedOutput('original', 0) and
        output.getIndexedOutput('mutated', 0)):
        allele = extractor.describe_dna(output.getIndexedOutput('original', 0),
                                        output.getIndexedOutput('mutated', 0))
        extracted = '(skipped)'
        if allele:
            extracted = unicode(allele)

    else:
        extracted = ''
    print("%s: End description extractor." % datetime.now())

    # Todo: Generate the fancy HTML views for the proteins here instead of in
    #   `mutalyzer.variantchecker`.
    arguments = {
        'description'         : description,
        'messages'            : map(util.message_info, output.getMessages()),
        'summary'             : summary,
        'parse_error'         : parse_error,
        'errors'              : errors,
        'genomicDescription'  : genomic_description,
        'chromDescription'    : output.getIndexedOutput(
                                  'genomicChromDescription', 0),
        'genomicDNA'          : genomic_dna,
        'visualisation'       : output.getOutput('visualisation'),
        'descriptions'        : output.getOutput('descriptions'),
        'protDescriptions'    : output.getOutput('protDescriptions'),
        'oldProtein'          : output.getOutput('oldProteinFancy'),
        'altStart'            : output.getIndexedOutput('altStart', 0),
        'altProtein'          : output.getOutput('altProteinFancy'),
        'newProtein'          : output.getOutput('newProteinFancy'),
        'transcriptInfo'      : output.getIndexedOutput('hasTranscriptInfo',
                                                        0, False),
        'transcriptCoding'    : output.getIndexedOutput('transcriptCoding', 0,
                                                        False),
        'exonInfo'            : output.getOutput('exonInfo'),
        'cdsStart_g'          : output.getIndexedOutput('cdsStart_g', 0),
        'cdsStart_c'          : output.getIndexedOutput('cdsStart_c', 0),
        'cdsStop_g'           : output.getIndexedOutput('cdsStop_g', 0),
        'cdsStop_c'           : output.getIndexedOutput('cdsStop_c', 0),
        'restrictionSites'    : output.getOutput('restrictionSites'),
        'legends'             : output.getOutput('legends'),
        'reference_filename'  : reference_filename,  # Todo: Download link is not shown...
        'browserLink'         : browser_link,
        'extractedDescription': extracted,
        'standalone'          : bool(request.args.get('standalone'))
    }

    output.addMessage(__file__, -1, 'INFO',
                      'Finished variant %s' % description)

    return render_template('name-checker.html', **arguments)


@website.route('/bed')
def bed():
    """
    Create a BED track for the given variant, listing the positions of its raw
    variants, e.g., for use in the UCSC Genome Browser.

    This basically just runs the variant checker and extracts the raw variants
    with positions.
    """
    # Backwards compatibility.
    if 'name' in request.args:
        return redirect(url_for('.bed',
                                description=request.args['name']),
                        code=301)

    description = request.args.get('description')

    if not description:
        abort(404)

    output = Output(__file__)

    variantchecker.check_variant(description, output)

    raw_variants = output.getIndexedOutput('rawVariantsChromosomal', 0)
    if not raw_variants:
        abort(404)

    # Todo: Hard-coded hg19.
    fields = {
        'name'       : 'Mutalyzer',
        'description': 'Mutalyzer track for ' + description,
        'visibility' : 'pack',
        'db'         : 'hg19',
        'url'        : url_for('.name_checker',
                               description=description,
                               _external=True),
        'color':       '255,0,0'}

    bed = ' '.join(['track'] +
                   ['%s="%s"' % field for field in fields.items()]) + '\n'

    for descr, positions in raw_variants[2]:
        bed += '\t'.join([raw_variants[0],
                          unicode(min(positions) - 1),
                          unicode(max(positions)),
                          descr,
                          '0',
                          raw_variants[1]]) + '\n'

    response = make_response(bed)
    response.headers['Content-Type'] = 'text/plain; charset=utf-8'
    return response


@website.route('/position-converter')
def position_converter():
    """
    Position converter.
    """
    # Backwards compatibility.
    if 'variant' in request.args:
        return redirect(url_for('.position_converter',
                                description=request.args['variant']),
                        code=301)

    assemblies = Assembly.query \
        .order_by(*Assembly.order_by_criteria) \
        .all()

    assembly_name_or_alias = request.args.get('assembly_name_or_alias',
                                              settings.DEFAULT_ASSEMBLY)
    description = request.args.get('description')

    if not description:
        return render_template('position-converter.html',
                               assemblies=assemblies,
                               assembly_name_or_alias=assembly_name_or_alias)

    output = Output(__file__)
    output.addMessage(__file__, -1, 'INFO',
                      'Received request positionConverter(%s, %s) from %s'
                      % (assembly_name_or_alias, description,
                         request.remote_addr))
    stats.increment_counter('position-converter/website')

    chromosomal_description = None
    transcript_descriptions = None

    try:
        assembly = Assembly.by_name_or_alias(assembly_name_or_alias)
    except NoResultFound:
        output.addMessage(__file__, 3, 'ENOASSEMBLY',
                          'Not a valid assembly.')
    else:
        converter = Converter(assembly, output)

        # Convert chromosome name to accession number.
        corrected_description = converter.correctChrVariant(description)

        if corrected_description:
            # Now we're ready to actually do position conversion.
            if not(':c.' in corrected_description or
                   ':n.' in corrected_description or
                   ':g.' in corrected_description or
                   ':m.' in corrected_description):
                grammar = Grammar(output)
                grammar.parse(corrected_description)

            if (':c.' in corrected_description or
                ':n.' in corrected_description):
                corrected_description = converter.c2chrom(
                        corrected_description)

            chromosomal_description = corrected_description

            if corrected_description and (':g.' in corrected_description or
                                          ':m.' in corrected_description):
                descriptions = converter.chrom2c(corrected_description, 'dict')
                if descriptions is None:
                    chromosomal_description = None
                elif descriptions:
                    transcript_descriptions = [
                        '%-10s:\t%s' % (key[:10], '\n\t\t'.join(value))
                        for key, value in descriptions.items()]

    messages = map(util.message_info, output.getMessages())

    output.addMessage(__file__, -1, 'INFO',
                      'Finished request positionConverter(%s, %s)'
                      % (assembly_name_or_alias, description))

    return render_template('position-converter.html',
                           assemblies=assemblies,
                           assembly_name_or_alias=assembly_name_or_alias,
                           description=description,
                           chromosomal_description=chromosomal_description,
                           transcript_descriptions=transcript_descriptions,
                           messages=messages)


@website.route('/snp-converter')
def snp_converter():
    """
    SNP converter.

    Convert a dbSNP rs number to HGVS description(s) of the SNP specified on
    the reference sequence(s) used by dbSNP.
    """
    # Backwards compatibility.
    if 'rsId' in request.args:
        return redirect(url_for('.snp_converter',
                                rs_id=request.args['rsId']),
                        code=301)

    rs_id = request.args.get('rs_id')

    if not rs_id:
        return render_template('snp-converter.html')

    output = Output(__file__)
    output.addMessage(__file__, -1, 'INFO',
                      'Received request snpConvert(%s) from %s'
                      % (rs_id, request.remote_addr))
    stats.increment_counter('snp-converter/website')

    descriptions = ncbi.rsid_to_descriptions(rs_id, output)

    messages = map(util.message_info, output.getMessages())

    output.addMessage(__file__, -1, 'INFO',
                      'Finished request snpConvert(%s)' % rs_id)

    return render_template('snp-converter.html',
                           rs_id=rs_id,
                           descriptions=descriptions,
                           messages=messages,
                           summary=output.Summary()[2])


@website.route('/reference-loader')
def reference_loader():
    """
    Reference sequence loader form.
    """
    assemblies = Assembly.query \
        .order_by(*Assembly.order_by_criteria) \
        .all()

    return render_template('reference-loader.html',
                           assemblies=assemblies,
                           assembly_name_or_alias=settings.DEFAULT_ASSEMBLY,
                           max_file_size=settings.MAX_FILE_SIZE // 1048576)


@website.route('/reference-loader', methods=['POST'])
def reference_loader_submit():
    """
    Reference sequence loader.

    There are five ways for the user to load a reference sequence,
    corresponding to values for the `method` field, each requiring some
    additional fields to be defined.:

    `method=upload_method`
      The reference sequence file is uploaded from a local file.

      - `file`: Reference sequence file to upload.

    `method=url_method`
      The reference sequence file can be found at the specified URL.

      - `url`: URL of reference sequence file to load.

    `method=slice_gene_method`
      Retrieve part of the reference genome for an HGNC gene symbol.

      - `genesymbol`: Gene symbol.
      - `organism`: Organism.
      - `upstream`: Number of 5' flanking nucleotides.
      - `downstream`: Number of 3' flanking nucleotides.

    `method=slice_accession_method`
      Retrieve a range of a chromosome by accession number.

      - `accession`: Chromosome Accession Number.
      - `accession_start`: Start position (one-based, inclusive, in reference
          orientation).
      - `accession_stop`: Stop position (one-based, inclusive, in reference
          orientation).
      - `accession_orientation`: Orientation.

    `method=slice_chromosome_method`
      Retrieve a range of a chromosome by name.

      - `assembly_name_or_alias`: Genome assembly by name or by alias.
      - `chromosome`: Chromosome name.
      - `chromosome_start`: Start position (one-based, inclusive, in reference
          orientation).
      - `chromosome_stop`: Stop position (one-based, inclusive, in reference
          orientation).
      - `chromosome_orientation`: Orientation.
    """
    method = request.form.get('method')

    output = Output(__file__)
    output.addMessage(__file__, -1, 'INFO',
                      'Received request upload(%s) with arguments %s from %s'
                      % (method, unicode(request.form), request.remote_addr))

    assemblies = Assembly.query \
        .order_by(*Assembly.order_by_criteria) \
        .all()

    retriever = Retriever.GenBankRetriever(output)
    ud, errors = '', []

    class InputException(Exception):
        pass

    def check_position(position, field):
        position = position.replace(',', '').replace('.', '').replace('-', '')
        try:
            return int(position)
        except AttributeError, ValueError:
            raise InputException('Expected an integer in field: %s' % field)

    try:
        if method == 'upload_method':
            # Todo: Non-conforming clients (read: LOVD) might send the form
            #   request urlencoded (and not as the requested multipart/
            #   form-data).
            file = request.files.get('file')
            if not file:
                raise InputException('Please select a local file for upload.')

            ud = retriever.uploadrecord(file.read())

        elif method == 'url_method':
            ud = retriever.downloadrecord(request.form.get('url'))

        elif method == 'slice_gene_method':
            genesymbol = request.form.get('genesymbol')
            organism = request.form.get('organism')
            upstream = check_position(request.form.get('upstream', ''),
                                      '5\' flanking nucleotides')
            downstream = check_position(request.form.get('downstream', ''),
                                        '3\' flanking nucleotides')
            ud = retriever.retrievegene(genesymbol, organism, upstream,
                                        downstream)

        elif method == 'slice_accession_method':
            accession = request.form.get('accession')
            start = check_position(request.form.get('accession_start', ''),
                                   'Start position')
            stop = check_position(request.form.get('accession_stop', ''),
                                  'Stop position')
            orientation = int(request.form.get('accession_orientation'))
            ud = retriever.retrieveslice(accession, start, stop, orientation)

        elif method == 'slice_chromosome_method':
            chromosome_name = request.form.get('chromosome')
            start = check_position(request.form.get('chromosome_start', ''),
                                   'Start position')
            stop = check_position(request.form.get('chromosome_stop', ''),
                                  'Stop position')
            orientation = int(request.form.get('chromosome_orientation'))

            assembly_name_or_alias = request.form.get('assembly_name_or_alias',
                                                      settings.DEFAULT_ASSEMBLY)
            try:
                assembly = Assembly.by_name_or_alias(assembly_name_or_alias)
            except NoResultFound:
                raise InputException('Invalid assembly')

            if not chromosome_name.startswith('chr'):
                chromosome_name = 'chr%s' % chromosome_name

            chromosome = assembly.chromosomes \
                .filter_by(name=chromosome_name) \
                .first()
            if not chromosome:
                raise InputException('Chromosome not available for assembly '
                                     '%s: %s' % (assembly.name, name))

            ud = retriever.retrieveslice(chromosome.accession, start, stop,
                                         orientation)

        else:
            raise InputException('Invalid method')

    except InputException as e:
        errors.append(e)

    if not ud:
        errors.append('The request could not be completed')
        errors.extend(unicode(m) for m in output.getMessages())

    output.addMessage(__file__, -1, 'INFO',
                      'Finished request upload(%s) with arguments %s from %s'
                      % (method, unicode(request.form), request.remote_addr))

    return render_template('reference-loader.html',
                           assemblies=assemblies,
                           assembly_name_or_alias=settings.DEFAULT_ASSEMBLY,
                           max_file_size=settings.MAX_FILE_SIZE // 1048576,
                           ud=ud,
                           errors=errors)


@website.route('/back-translator')
def back_translator():
    """
    Back translator.
    """
    output = Output(__file__)
    output.addMessage(
        __file__, -1, 'INFO',
        'Received Back Translate request from {}'.format(request.remote_addr))
    stats.increment_counter('back-translator/website')

    description = request.args.get('description')

    variants = []
    if description:
        variants = backtranslator.backtranslate(output, description)

    errors, warnings, summary = output.Summary()
    messages = map(util.message_info, output.getMessages())

    output.addMessage(__file__, -1, 'INFO', 'Finished Back Translate request')

    return render_template(
        'back-translator.html', errors=errors, summary=summary,
        description=description or '', messages=messages, variants=variants)


@website.route('/description-extractor')
def description_extractor():
    """
    Description extractor loader form.
    """
    return render_template('description-extractor.html',
                           extractor_max_input_length=settings.EXTRACTOR_MAX_INPUT_LENGTH)


@website.route('/description-extractor', methods=['POST'])
def description_extractor_submit():
    """
    The Variant Description Extractor (experimental service).

    There multiple ways for the user to provide two sequences, corresponding to
    the values for the `reference_method` and `sample_method` fields, each
    requiring some additional fields to be defined:

    `raw_method`
      The reference and sample sequences are pasted into the form fields.

      - `reference_sequence`: The reference sequence.
      - `sample_sequence`: The sample sequence.

    `file_method`
      The reference and sample sequences are uploaded.

      - `reference_file`: The reference file.
      - `sample_file`: The sample file.

    `refseq_method`
      The reference and sample sequences are given by RefSeq accession numbers.

      - `reference_accession_number`: RefSeq accession number for the reference
        sequence.
      - `sample_accession_number`: RefSeq accession number for the sample
        sequence.
    """
    output = Output(__file__)
    output.addMessage(__file__, -1, 'INFO',
                      'Received Description Extract request from %s'
                      % request.remote_addr)
    stats.increment_counter('description-extractor/website')

    r = s = ''
    reference_method = request.form.get('reference_method')
    sample_method = request.form.get('sample_method')
    reference_sequence = request.form.get('reference_sequence')
    sample_sequence = request.form.get('sample_sequence')
    reference_file = request.files.get('reference_file')
    sample_file = request.files.get('sample_file')
    reference_filename = ''
    sample_filename = ''
    reference_accession_number = request.form.get('reference_accession_number')
    sample_accession_number = request.form.get('sample_accession_number')

    if reference_method == 'refseq_method':
        if reference_accession_number:
            retriever = Retriever.GenBankRetriever(output)
            genbank_record = retriever.loadrecord(reference_accession_number)
            if genbank_record:
                r = unicode(genbank_record.seq)
        else:
            output.addMessage(__file__, 3, 'EEMPTYFIELD',
                'Reference accession number input fields is empty.')
    elif reference_method == 'file_method':
        if reference_file:
            reference_filename = reference_file.filename
            r = util.read_dna(reference_file)
        else:
            output.addMessage(__file__, 3, 'EEMPTYFIELD',
                'No reference file provided.')
    else: # raw_method
        if reference_sequence:
            r = util.read_dna(StringIO.StringIO(reference_sequence))
        else:
            output.addMessage(__file__, 3, 'EEMPTYFIELD',
                'Reference sequence number input fields is empty.')

    if sample_method == 'refseq_method':
        if sample_accession_number:
            retriever = Retriever.GenBankRetriever(output)
            genbank_record = retriever.loadrecord(sample_accession_number)
            if genbank_record:
                s = unicode(genbank_record.seq)
        else:
            output.addMessage(__file__, 3, 'EEMPTYFIELD',
                'Sample accession number input fields is empty.')
    elif sample_method == 'file_method':
        if sample_file:
            sample_filename = sample_file.filename
            s = util.read_dna(sample_file)
        else:
            output.addMessage(__file__, 3, 'EEMPTYFIELD',
                'No sample file provided.')
    else: # raw_method
        if sample_sequence:
            s = util.read_dna(StringIO.StringIO(sample_sequence))
        else:
            output.addMessage(__file__, 3, 'EEMPTYFIELD',
                'Sample sequence number input fields is empty.')

    # Todo: Move this to the describe module.
    if not r or not util.is_dna(r):
        output.addMessage(__file__, 3, 'ENODNA',
                          'Reference sequence is not DNA.')
    if not s or not util.is_dna(s):
        output.addMessage(__file__, 3, 'ENODNA',
                          'Sample sequence is not DNA.')

    raw_vars = None
    if r and s:
        if (len(r) > settings.EXTRACTOR_MAX_INPUT_LENGTH or
            len(s) > settings.EXTRACTOR_MAX_INPUT_LENGTH):
            output.addMessage(__file__, 3, 'EMAXSIZE',
                              'Input sequences are restricted to {:,} bp.'
                              .format(settings.EXTRACTOR_MAX_INPUT_LENGTH))
        else:
            raw_vars = extractor.describe_dna(r, s)

    errors, warnings, summary = output.Summary()
    messages = map(util.message_info, output.getMessages())

    output.addMessage(__file__, -1, 'INFO',
                      'Finished Description Extract request')

    return render_template('description-extractor.html',
        extractor_max_input_length=settings.EXTRACTOR_MAX_INPUT_LENGTH,
        reference_sequence=reference_sequence or '',
        sample_sequence=sample_sequence or '',
        reference_accession_number=reference_accession_number or '',
        sample_accession_number=sample_accession_number or '',
        reference_filename=reference_filename or '',
        sample_filename=sample_filename or '',
        raw_vars=raw_vars, errors=errors, summary=summary, messages=messages,
        reference_method=reference_method, sample_method=sample_method)


@website.route('/reference/<string:filename>')
def reference(filename):
    """
    Download reference file from cache.
    """
    file_path = os.path.join(settings.CACHE_DIR, '%s.bz2' % filename)

    if not os.path.isfile(file_path):
        abort(404)

    response = make_response(bz2.BZ2File(file_path, 'r').read())

    response.headers['Content-Type'] = 'text/plain; charset=utf-8'
    response.headers['Content-Disposition'] = ('attachment; filename="%s"'
                                               % filename)
    return response


@website.route('/batch-jobs')
def batch_jobs():
    """
    Batch jobs form.
    """
    job_type = request.args.get('job_type', 'name-checker')

    assemblies = Assembly.query \
        .order_by(*Assembly.order_by_criteria) \
        .all()
    assembly_name_or_alias = request.args.get('assembly_name_or_alias',
                                              settings.DEFAULT_ASSEMBLY)

    return render_template('batch-jobs.html',
                           assemblies=assemblies,
                           assembly_name_or_alias=assembly_name_or_alias,
                           job_type=job_type,
                           max_file_size=settings.MAX_FILE_SIZE // 1048576)


@website.route('/batch-jobs', methods=['POST'])
def batch_jobs_submit():
    """
    Run batch jobs and render batch checker HTML form. The batch jobs are
    added to the database by the scheduler and ran by the BatchChecker
    daemon.
    """
    job_type = request.form.get('job_type')
    email = request.form.get('email')

    # Note that this is always a seekable binary file object.
    batch_file = request.files.get('file')

    assemblies = Assembly.query \
        .order_by(*Assembly.order_by_criteria) \
        .all()
    assembly_name_or_alias = request.form.get('assembly_name_or_alias',
                                              settings.DEFAULT_ASSEMBLY)

    errors = []

    if not email:
        email = '{}@website.mutalyzer'.format(request.remote_addr)

    if job_type not in BATCH_JOB_TYPES:
        errors.append('Invalid batch job type.')

    if not file:
        errors.append('Please select a local file for upload.')

    if job_type == 'position-converter':
        try:
            Assembly.by_name_or_alias(assembly_name_or_alias)
        except NoResultFound:
            errors.append('Not a valid assembly.')
        argument = assembly_name_or_alias
    else:
        argument = None

    output = Output(__file__)

    if not errors:
        stats.increment_counter('batch-job/website')

        scheduler = Scheduler.Scheduler()
        file_instance = File.File(output)
        job, columns = file_instance.parseBatchFile(batch_file)

        if job is None:
            errors.append('Could not parse input file, please check your '
                          'file format.')
        else:
            result_id = scheduler.addJob(email, job, columns, job_type,
                                         argument=argument)

            # Todo: We now assume that the job was not scheduled if there are
            #   messages, which is probably not correct.
            if not output.getMessages():
                return redirect(url_for('.batch_job_progress',
                                        result_id=result_id))

    for error in errors:
        output.addMessage(__file__, 3, 'EBATCHJOB', error)

    messages = map(util.message_info, output.getMessages())

    return render_template('batch-jobs.html',
                           assemblies=assemblies,
                           assembly_name_or_alias=assembly_name_or_alias,
                           job_type=job_type,
                           max_file_size=settings.MAX_FILE_SIZE // 1048576,
                           messages=messages)


@website.route('/batch-job-progress')
def batch_job_progress():
    """
    Batch jobs progress viewer.
    """
    result_id = request.args.get('result_id')
    json = bool(request.args.get('json'))

    if not result_id:
        return render_template('batch-job-progress.html')

    batch_job = BatchJob.query.filter_by(result_id=result_id).first()

    if not batch_job:
        # Only now, the job can be complete. But since we don't keep completed
        # jobs in the database, we can only see if it ever existed by checking
        # the result file.
        path = os.path.join(settings.CACHE_DIR, 'batch-job-%s.txt' % result_id)
        if os.path.isfile(path):
            if json:
                return jsonify(items_left=1, complete=True)
            return render_template('batch-job-progress.html',
                                   result_id=result_id)
        else:
            return render_template('batch-job-progress.html')

    items_left = batch_job.batch_queue_items.count()

    if json:
        return jsonify(items_left=items_left, complete=False)
    return render_template('batch-job-progress.html',
                           result_id=result_id,
                           items_left=items_left)


@website.route('/batch-job-result/batch-job-<string:result_id>.txt')
def batch_job_result(result_id):
    """
    Batch job result file download.
    """
    if not result_id:
        abort(404)

    batch_job = BatchJob.query.filter_by(result_id=result_id).first()
    if batch_job:
        # If the batch job exists, it is not done yet.
        abort(404)

    return send_from_directory(settings.CACHE_DIR,
                               'batch-job-%s.txt' % result_id,
                               mimetype='text/plain; charset=utf-8',
                               as_attachment=True)


# Todo: Is this obsolete?
@website.route('/getGS')
def lovd_get_gs():
    """
    LOVD bypass to get the correct GeneSymbol incl Transcript variant.

    Used by LOVD to get the correct transcript variant out of a genomic
    record. LOVD uses a genomic reference (``NC_``?) in combination with a
    gene symbol to pass variant info to mutalyzer. Mutalyzer 1.0 was only
    using the first transcript. LOVD supplies the NM of the transcript needed
    but this was ignored. This helper allows LOVD to get the requested
    transcript variant from a genomic reference.

    Parameters:

    mutationName
      The mutationname without gene symbol.
    variantRecord
      The NM reference of the variant.
    forward
      If set this forwards the request to the name checker.

    Returns: Output of name checker if `forward` is set, otherwise the
    gene symbol with the variant notation as string.
    """
    mutation_name = request.args['mutationName']
    variant_record = request.args['variantRecord']
    forward = request.args.get('forward')

    output = Output(__file__)
    output.addMessage(__file__, -1, 'INFO',
                      'Received request getGS(%s, %s, %s) from %s'
                      % (mutation_name, variant_record, forward,
                         request.remote_addr))

    variantchecker.check_variant(mutation_name, output)

    output.addMessage(__file__, -1, 'INFO',
                      'Finished request getGS(%s, %s, %s)'
                      % (mutation_name, variant_record, forward))

    legends = output.getOutput('legends')

    # Filter the transcript from the legend.
    legends = [l for l in legends if '_v' in l[0]]
    for l in legends:
        if l[1] == variant_record:
            if forward:
                p, a = mutation_name.split(':')
                return redirect(url_for('.name_checker',
                                        description='%s(%s):%s' % (p, l[0], a),
                                        standalone=1))
            else:
                response = make_response(l[0])
                response.headers['Content-Type'] = 'text/plain; charset=utf-8'
                return response

    response = make_response('Transcript not found')
    response.headers['Content-Type'] = 'text/plain; charset=utf-8'
    return response


# Todo: Is this obsolete?
@website.route('/Variant_info')
def lovd_variant_info():
    """
    The chromosomal to coding and vice versa conversion interface for LOVD.

    Search for an NM number in the database, if the version number matches,
    get the start and end positions in a variant and translate these positions
    to chromosomal notation if the variant is in coding notation and vice
    versa.

    - If no end position is present, the start position is assumed to be the
      end position.
    - If the version number is not found in the database, an error message is
      generated and a suggestion for an other version is given.
    - If the reference sequence is not found at all, an error is returned.
    - If no variant is present, the transcription start and end and CDS end
      in coding notation is returned.
    - If the variant is not accepted by the nomenclature parser, a parse error
      will be printed.

    Get variant info and return the result as plain text.

    Parameters:

    LOVD_ver
      The version of the calling LOVD.
    build
      The human genome build (hg19 assumed).
    acc
      The accession number (NM number).
    var
      A description of the variant.

    Returns:

    start_main
      The main coordinate of the start position in I{c.} (non-star) notation.
    start_offset
      The offset coordinate of the start position in I{c.} notation (intronic
      position).
    end_main
      The main coordinate of the end position in I{c.} (non-star) notation.
    end_offset
      The offset coordinate of the end position in I{c.} notation (intronic
      position).
    start_g
      The I{g.} notation of the start position.
    end_g
      The I{g.} notation of the end position.
    type
      The mutation type.

    Returns (alternative):

    trans_start
      Transcription start in I{c.} notation.
    trans_stop
      Transcription stop in I{c.} notation.
    CDS_stop
      CDS stop in I{c.} notation.
    """
    lovd_version = request.args['LOVD_ver']
    build = request.args['build']
    accession = request.args['acc']
    description = request.args.get('var')

    output = Output(__file__)
    output.addMessage(__file__, -1, 'INFO',
                      'Received request variantInfo(%s:%s (LOVD_ver %s, '
                      'build %s)) from %s'
                      % (accession, description, lovd_version, build,
                         request.remote_addr))

    try:
        assembly = Assembly.by_name_or_alias(build)
    except NoResultFound:
        response = make_response('invalid build')
        response.headers['Content-Type'] = 'text/plain; charset=utf-8'
        return response

    converter = Converter(assembly, output)

    result = ''

    # If no variant is given, return transcription start, transcription
    # end and CDS stop in c. notation.
    if description:
        ret = converter.mainMapping(accession, description)
    else:
        ret = converter.giveInfo(accession)
        if ret:
            result = '%i\n%i\n%i' % ret

    if not result and not getattr(ret, 'startmain', None):
        out = output.getOutput('LOVDERR')
        if out:
            result = out[0]
        else:
            result = 'Unknown error occured'

    output.addMessage(__file__, -1, 'INFO',
                      'Finished request variantInfo(%s:%s (LOVD_ver %s, '
                      'build %s))'
                      % (accession, description, lovd_version, build))

    if not result and getattr(ret, 'startmain', None):
        result = '%i\n%i\n%i\n%i\n%i\n%i\n%s' % (
            ret.startmain, ret.startoffset, ret.endmain, ret.endoffset,
            ret.start_g, ret.end_g, ret.mutationType)

    # Todo: Obsoleted error messages, remove soon.
    if lovd_version == '2.0-23':
        response = re.sub('^Error \(.*\):', 'Error:', result)

    response = make_response(result)
    response.headers['Content-Type'] = 'text/plain; charset=utf-8'
    return response


# Register redirects for backwards compatibility.
website.add_url_rule('/index', view_func=homepage, alias=True)
website.add_url_rule('/nameGenerator', view_func=name_generator, alias=True)
website.add_url_rule('/syntaxCheck', view_func=syntax_checker, alias=True)
website.add_url_rule('/check', view_func=name_checker, alias=True)
website.add_url_rule('/checkForward', view_func=name_checker, alias=True)
website.add_url_rule('/positionConverter', view_func=position_converter,
                     alias=True)
website.add_url_rule('/snp', view_func=snp_converter, alias=True)
website.add_url_rule('/upload', view_func=reference_loader, alias=True)
website.add_url_rule('/descriptionExtract', view_func=description_extractor,
                     alias=True)
website.add_url_rule('/Reference/<string:filename>', view_func=reference,
                     alias=True)
website.add_url_rule('/batch', view_func=batch_jobs, alias=True)
website.add_url_rule('/batchNameChecker', view_func=batch_jobs, alias=True)
website.add_url_rule('/batchSyntaxChecker', view_func=batch_jobs, alias=True)
website.add_url_rule('/batchPositionConverter', view_func=batch_jobs,
                     alias=True)
website.add_url_rule('/batchSnpConverter', view_func=batch_jobs, alias=True)
