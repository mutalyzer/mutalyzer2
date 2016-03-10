"""
Tests for the mutalyzer.website module.
"""


from __future__ import unicode_literals

import bz2
from mock import patch
import os
from io import BytesIO
import re
import urlparse

from Bio import Entrez
import lxml.html
import pytest

from mutalyzer import announce, Scheduler
from mutalyzer.db.models import BatchJob
from mutalyzer.website import create_app

from fixtures import with_references


# TODO: Tests for /upload.


def get_links(data, path=None):
    """
    Extract all link targets, or only those targeting the specified path, from
    the page and parse them.
    """
    def parse_link(link):
        splitted = urlparse.urlsplit(link)
        return splitted.path, urlparse.parse_qs(splitted.query)

    links = [parse_link(link) for link in re.findall('href="([^"]+)"', data)]
    return [(p, q) for p, q in links if path is None or p == path]


@pytest.fixture
def website():
    return create_app().test_client()


def test_homepage(website):
    """
    Expect the index HTML page.
    """
    r = website.get('/')
    assert r.status_code == 200
    assert 'Welcome to the Mutalyzer website' in r.data


def test_about(website):
    """
    See if people get proper credit.
    """
    r = website.get('/about')
    assert r.status == '200 OK'
    assert 'Jonathan Vis' in r.data


def test_non_existing(website):
    """
    Expect a 404 response.
    """
    r = website.get('/this/doesnotexist')
    assert r.status_code == 404


@pytest.mark.usefixtures('db')
def test_menu_links(website):
    """
    Test all links in the main menu.
    """
    # This could contain relative links we want to skip.
    ignore = []
    r = website.get('/')

    dom = lxml.html.fromstring(r.data)

    for link in dom.cssselect('nav a'):
        href = link.get('href')
        if (href.startswith('http://') or
            href.startswith('https://') or
            href.startswith('mailto:') or
            href.startswith('#') or
            href in ignore):
            continue
        if not href.startswith('/'):
            href = '/' + href

        r = website.get(href)
        assert r.status_code == 200


def test_announcement(website):
    """
    We should always see the current announcement.
    """
    announce.set_announcement('Test announcement')
    r = website.get('/syntax-checker')
    assert r.status == '200 OK'
    assert 'Test announcement' in r.data

    announce.set_announcement('New announcement')
    r = website.get('/syntax-checker')
    assert r.status == '200 OK'
    assert 'New announcement' in r.data

    announce.unset_announcement()
    r = website.get('/syntax-checker')
    assert r.status == '200 OK'
    assert 'nnouncement' not in r.data


def test_description_extractor_raw(website):
    """
    Submit two sequences to the variant description extractor.
    """
    r = website.post('/description-extractor', data={
        'reference_method': 'raw_method',
        'sample_method': 'raw_method',
        'reference_sequence': 'ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA',
        'sample_sequence': 'ATGATTTGATCAGATACATGTGATACCGGTAGTTAGGACAA'})
    assert '[5_6insTT;17del;26A&gt;C;35dup]' in r.data


def test_description_extractor_raw_fastq(website):
    """
    Submit two sequences to the variant description extractor.
    """
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'data', 'extractor_input.fq')
    r = website.post('/description-extractor', data={
        'reference_method': 'raw_method',
        'sample_method': 'raw_method',
        'reference_sequence': 'ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA',
        'sample_sequence': open(path).read()})
    assert '[5_6insTT;17del;26A&gt;C;35dup]' in r.data


@with_references('NM_004006.1', 'NM_004006.2')
def test_description_extractor_refseq(website):
    """
    Submit two accession numbers to the variant description extractor.
    """
    r = website.post('/description-extractor', data={
        'reference_method': 'refseq_method',
        'sample_method': 'refseq_method',
        'reference_accession_number': 'NM_004006.1',
        'sample_accession_number': 'NM_004006.2'})
    assert '[12749G&gt;A;13729G&gt;A]' in r.data


def test_description_extractor_file_fasta(website):
    """
    Submit a sequence and a FASTA file to the variant description
    extractor.
    """
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'data', 'extractor_input.fa')
    r = website.post('/description-extractor', data={
        'reference_method': 'raw_method',
        'sample_method': 'file_method',
        'reference_sequence': 'ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA',
        'sample_file': (open(path), 'extractor_input.fa')})
    assert '[5_6insTT;17del;26A&gt;C;35dup]' in r.data


def test_description_extractor_file_fastq(website):
    """
    Submit a sequence and a FASTQ file to the variant description
    extractor.
    """
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'data', 'extractor_input.fq')
    r = website.post('/description-extractor', data={
        'reference_method': 'raw_method',
        'sample_method': 'file_method',
        'reference_sequence': 'ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA',
        'sample_file': (open(path), 'extractor_input.fq')})
    assert '[5_6insTT;17del;26A&gt;C;35dup]' in r.data


def test_description_extractor_file_text(website):
    """
    Submit a sequence and a text file to the variant description
    extractor.
    """
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'data', 'extractor_input.txt')
    r = website.post('/description-extractor', data={
        'reference_method': 'raw_method',
        'sample_method': 'file_method',
        'reference_sequence': 'ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA',
        'sample_file': (open(path), 'extractor_input.txt')})
    assert '[5_6insTT;17del;26A&gt;C;35dup]' in r.data


def test_description_extractor_ref_too_long(settings, website):
    """
    Submit a reference sequence exceeding the maximum length to the variant
    description extractor.
    """
    r = website.post('/description-extractor', data={
        'reference_method': 'raw_method',
        'sample_method': 'raw_method',
        'reference_sequence': 'A' * (settings.EXTRACTOR_MAX_INPUT_LENGTH + 1),
        'sample_sequence': 'A'})
    assert '2_{}del'.format(settings.EXTRACTOR_MAX_INPUT_LENGTH + 1) not in r.data
    assert 'Input sequences are restricted to ' in r.data
    assert '1 Error, 0 Warnings.' in r.data


def test_description_extractor_sample_too_long(settings, website):
    """
    Submit a sample sequence exceeding the maximum length to the variant
    description extractor.
    """
    r = website.post('/description-extractor', data={
        'reference_method': 'raw_method',
        'sample_method': 'raw_method',
        'reference_sequence': 'A' * (settings.EXTRACTOR_MAX_INPUT_LENGTH),
        'sample_sequence': 'A' * (settings.EXTRACTOR_MAX_INPUT_LENGTH + 1)})
    assert '{}dup'.format(settings.EXTRACTOR_MAX_INPUT_LENGTH) not in r.data
    assert 'Input sequences are restricted to ' in r.data
    assert '1 Error, 0 Warnings.' in r.data


def test_description_extractor_lowercase(website):
    """
    Submit a sample sequence with a base in lowercase to the variant
    description extractor.
    """
    r = website.post('/description-extractor', data={
        'reference_method': 'raw_method',
        'sample_method': 'raw_method',
        'reference_sequence': 'TTT',
        'sample_sequence': 'TaT'})
    assert '<pre class="description">2T&gt;A</pre>' in r.data


def test_checksyntax_valid(website):
    """
    Submit the check syntax form with a valid variant.
    """
    r = website.get('/syntax-checker',
                    query_string={'description': 'AB026906.1:c.274G>T'})
    assert 'The syntax of this variant description is OK!' in r.data


def test_checksyntax_invalid(website):
    """
    Submit the check syntax form with an invalid variant.
    """
    r = website.get('/syntax-checker',
                    query_string={'description': 'AB026906.1:c.27'})
    assert 'Fatal' in r.data
    assert 'The &quot;^&quot; indicates the position where the error occurred' in r.data


@with_references('NM_002001.2')
def test_check_valid(website):
    """
    Submit the name checker form with a valid variant.
    Should include form and main HTML layout.
    """
    r = website.get('/name-checker',
                    query_string={'description': 'NM_002001.2:g.1del'})
    assert '0 Errors' in r.data
    assert '0 Warnings' in r.data
    assert 'Raw variant 1: deletion of 1' in r.data
    assert 'value="NM_002001.2:g.1del"' in r.data


def test_check_invalid(website):
    """
    Submit the name checker form with an invalid variant.
    """
    r = website.get('/name-checker',
                    query_string={'description': 'NM_002001.2'})
    assert '1 Error' in r.data
    assert '0 Warnings' in r.data
    assert 'The &quot;^&quot; indicates the position where the error occurred' in r.data


@with_references('NP_064445.1')
def test_check_protein_reference(website):
    """
    Submit the name checker form with a protein reference sequence (not
    supported).
    """
    r = website.get('/name-checker',
                    query_string={'description': 'NP_064445.1:c.274G>T'})
    assert '1 Error' in r.data
    assert '0 Warnings' in r.data
    assert 'Protein reference sequences are not supported' in r.data


@with_references('NM_002001.2')
def test_check_noninteractive(website):
    """
    Submit the name checker form non-interactively.
    Should not include form and main layout HTML.
    """
    r = website.get('/name-checker',
                    query_string={'description': 'NM_002001.2:g.1del',
                                  'standalone': '1'})
    assert '<a href="#bottom" class="hornav">go to bottom</a>' not in r.data
    assert '<input type="text" name="description" value="NM_002001.2:g.1del" style="width:100%">' not in r.data
    assert '0 Errors' in r.data
    assert '0 Warnings' in r.data
    assert 'Raw variant 1: deletion of 1' in r.data


@with_references('NG_012772.1')
def test_check_noninteractive_links(website):
    """
    Submitting non-interactively should have links to transcripts also
    non-interactive.
    """
    r = website.get('/name-checker',
                    query_string={'description': 'NG_012772.1:g.128del',
                                  'standalone': '1'})
    assert '0 Errors' in r.data

    links = get_links(r.data, path='/name-checker')
    assert len(links) >= 2
    assert all(q['standalone'] == ['1'] for _, q in links)


@with_references('NG_012772.1')
def test_check_interactive_links(website):
    """
    Submitting interactively should have links to transcripts also
    interactive.
    """
    r = website.get('/name-checker',
                    query_string={'description': 'NG_012772.1:g.128del'})
    assert '0 Errors' in r.data

    links = get_links(r.data, path='/name-checker')
    assert len(links) >= 2
    assert all('standalone' not in q for _, q in links)


def test_snp_converter_valid(website):
    """
    Submit the SNP converter form with a valid SNP.
    """
    # Patch Retriever.snpConvert to return rs9919552.
    def mock_efetch(*args, **kwargs):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data',
                            'rs9919552.xml.bz2')
        return bz2.BZ2File(path)

    with patch.object(Entrez, 'efetch', mock_efetch):
        r = website.get('/snp-converter',
                        query_string={'rs_id': 'rs9919552'})
    assert '0 Errors' in r.data
    assert '0 Warnings' in r.data
    assert 'NC_000011.9:g.111959625C&gt;T' in r.data
    assert 'NG_012337.2:g.7055C&gt;T' in r.data
    assert 'NM_003002.3:c.204C&gt;T' in r.data
    assert 'NP_002993.1:p.Ser68=' in r.data


def test_snp_converter_invalid(website):
    """
    Submit the SNP converter form with an invalid SNP.
    """
    r = website.get('/snp-converter',
                    query_string={'rs_id': 'r9919552'})

    assert '1 Error' in r.data
    assert '0 Warnings' in r.data
    assert 'Fatal' in r.data
    assert 'This is not a valid dbSNP id' in r.data


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_position_converter_c2g(website):
    """
    Submit the position converter form with a valid variant.
    """
    r = website.get('/position-converter',
                    query_string={'assembly_name_or_alias': 'hg19',
                                  'description': 'NM_003002.2:c.204C>T'})
    assert 'NC_000011.9:g.111959625C&gt;T' in r.data


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_position_converter_g2c(website):
    """
    Submit the position converter form with a valid variant.
    """
    r = website.get('/position-converter',
                    query_string={'assembly_name_or_alias': 'hg19',
                                  'description': 'NC_000011.9:g.111959625C>T'})
    assert 'NM_003002.2:c.204C&gt;T' in r.data


def _batch(website, job_type='name-checker', assembly_name_or_alias=None,
           file='', size=0, header='', lines=None):
    """
    Submit a batch form.

    @kwarg batch_type: Type of batch job to test. One of name-checker,
                       syntax-checker, position-converter.
    @kwarg argument: Optional extra argument for the batch job.
    @kwarg file: String with variants to use as input for the batch job.
    @kwarg size: Number of variants in input.
    @kwarg header: Message that must be found in the batch job result.
    @kwarg lines: Number of result rows expected.

    @return: The batch result document.
    @rtype: string
    """
    data = {'job_type': job_type,
            'email': 'test@test.test',
            'file': (BytesIO(file.encode('utf-8')), 'test.txt')}
    if assembly_name_or_alias is not None:
        data['assembly_name_or_alias'] = assembly_name_or_alias

    r = website.post('/batch-jobs',
                     data=data)
    progress_url = '/' + r.location.split('/')[-1]

    r = website.get(progress_url)
    assert '<div id="if_items_left">' in r.data
    assert '<div id="ifnot_items_left" style="display:none">' in r.data
    assert ('<span id="items_left">%d</span>' % size) in r.data

    scheduler = Scheduler.Scheduler()
    scheduler.process()

    r = website.get(progress_url)
    assert '<div id="if_items_left" style="display:none">' in r.data
    assert '<div id="ifnot_items_left">' in r.data

    dom = lxml.html.fromstring(r.data)
    result_url = dom.cssselect('#ifnot_items_left a')[0].attrib['href']

    if not lines:
        lines = size

    r = website.get(result_url)
    assert 'text/plain' in r.headers['Content-Type']
    assert header in r.data
    assert len(r.data.strip().split('\n')) - 1 == lines

    return r.data


@with_references('AB026906.1', 'NM_003002.2', 'AL449423.14')
def test_batch_namechecker(website):
    """
    Submit the batch name checker form.
    """
    variants = ['AB026906.1(SDHD):g.7872G>T',
                'NM_003002.2:c.3_4insG',
                'AL449423.14(CDKN2A_v002):c.5_400del']
    _batch(website,
           'name-checker',
           file='\n'.join(variants),
           size=len(variants),
           header='Input\tErrors and warnings')


@pytest.mark.usefixtures('db')
def test_batch_namechecker_extra_tab(website):
    """
    Submit the batch syntax checker form with lines ending with tab
    characters.
    """
    variants = ['AB026906.1(SDHD):g.7872G>T\t',
                'AB026906.1(SDHD):g.7872G>T\t',
                'AB026906.1(SDHD):g.7872G>T\t']
    _batch(website,
           'syntax-checker',
           file='\n'.join(variants),
           size=len(variants) * 2,
           lines=len(variants),
           header='Input\tStatus')


@pytest.mark.usefixtures('db')
def test_batch_syntaxchecker(website):
    """
    Submit the batch syntax checker form.
    """
    variants = ['AB026906.1(SDHD):g.7872G>T',
                'NM_003002.1:c.3_4insG',
                'AL449423.14(CDKN2A_v002):c.5_400del']
    _batch(website,
           'syntax-checker',
           file='\n'.join(variants),
           size=len(variants),
           header='Input\tStatus')


@pytest.mark.usefixtures('hg19')
def test_batch_positionconverter(website):
    """
    Submit the batch position converter form.
    """
    variants = ['NM_003002.2:c.204C>T',
                'NC_000011.9:g.111959625C>T']
    _batch(website,
           'position-converter',
           assembly_name_or_alias='hg19',
           file='\n'.join(variants),
           size=len(variants),
           header='Input Variant')


@pytest.mark.usefixtures('db')
def test_batch_syntaxchecker_newlines_unix(website):
    """
    Submit batch syntax checker job with Unix line endings.
    """
    variants = ['AB026906.1(SDHD):g.7872G>T',
                'NM_003002.1:c.3_4insG',
                'AL449423.14(CDKN2A_v002):c.5_400del']
    _batch(website,
           'syntax-checker',
           file='\n'.join(variants),
           size=len(variants),
           header='Input\tStatus')


@pytest.mark.usefixtures('db')
def test_batch_syntaxchecker_newlines_mac(website):
    """
    Submit batch syntax checker job with Mac line endings.
    """
    variants = ['AB026906.1(SDHD):g.7872G>T',
                'NM_003002.1:c.3_4insG',
                'AL449423.14(CDKN2A_v002):c.5_400del']
    _batch(website,
           'syntax-checker',
           file='\r'.join(variants),
           size=len(variants),
           header='Input\tStatus')


@pytest.mark.usefixtures('db')
def test_batch_syntaxchecker_newlines_windows(website):
    """
    Submit batch syntax checker job with Windows line endings.
    """
    variants = ['AB026906.1(SDHD):g.7872G>T',
                'NM_003002.1:c.3_4insG',
                'AL449423.14(CDKN2A_v002):c.5_400del']
    _batch(website,
           'syntax-checker',
           file='\r\n'.join(variants),
           size=len(variants),
           header='Input\tStatus')


@pytest.mark.usefixtures('db')
def test_batch_syntaxchecker_newlines_big_unix(website):
    """
    Submit big batch syntax checker job with Unix line endings.
    """
    samples = ['AB026906.1(SDHD):g.7872G>T',
               'NM_003002.1:c.3_4insG',
               'AL449423.14(CDKN2A_v002):c.5_400del']
    variants = []
    # Create 240 variants out of 3 samples
    for i in range(80):
        variants.extend(samples)
    _batch(website,
           'syntax-checker',
           file='\n'.join(variants),
           size=len(variants),
           header='Input\tStatus')


@pytest.mark.usefixtures('db')
def test_batch_syntaxchecker_newlines_big_mac(website):
    """
    Submit big batch syntax checker job with Mac line endings.
    """
    samples = ['AB026906.1(SDHD):g.7872G>T',
               'NM_003002.1:c.3_4insG',
               'AL449423.14(CDKN2A_v002):c.5_400del']
    variants = []
    # Create 240 variants out of 3 samples
    for i in range(80):
        variants.extend(samples)
    _batch(website,
           'syntax-checker',
           file='\r'.join(variants),
           size=len(variants),
           header='Input\tStatus')


@pytest.mark.usefixtures('db')
def test_batch_syntaxchecker_newlines_big_windows(website):
    """
    Submit big batch syntax checker job with Windows line endings.
    """
    samples = ['AB026906.1(SDHD):g.7872G>T',
               'NM_003002.1:c.3_4insG',
               'AL449423.14(CDKN2A_v002):c.5_400del']
    variants = []
    # Create 240 variants out of 3 samples
    for i in range(80):
        variants.extend(samples)
    _batch(website,
           'syntax-checker',
           file='\r\n'.join(variants),
           size=len(variants),
           header='Input\tStatus')


@pytest.mark.usefixtures('db')
def test_batch_syntaxchecker_oldstyle(website):
    """
    Submit the batch syntax checker form with old style input file.
    """
    variants = ['AccNo\tGenesymbol\tMutation',
                'AB026906.1\tSDHD\tg.7872G>T',
                'NM_003002.1\t\tc.3_4insG',
                'AL449423.14\tCDKN2A_v002\tc.5_400del']
    _batch(website,
           'syntax-checker',
           file='\n'.join(variants),
           size=len(variants)-1,
           header='Input\tStatus')


@with_references('AB026906.1')
def test_batch_namechecker_restriction_sites(website):
    """
    Submit the batch name checker form and see if restriction site effects
    are added.
    """
    variants = ['AB026906.1:c.274G>T',
                'AB026906.1:c.[274G>T;143A>G;15G>T]']
    results = _batch(website,
                     'name-checker',
                     file='\n'.join(variants),
                     size=len(variants),
                     header='Input\tErrors and warnings').strip().split('\n')
    assert 'Restriction Sites Created\tRestriction Sites Deleted' in results[0]
    assert 'CviQI,RsaI\tBccI' in results[1]
    assert 'CviQI,RsaI;HhaI,HinP1I;SfcI\tBccI;;BpmI,BsaXI (2),LpnPI,MnlI' in results[2]


@pytest.mark.usefixtures('db')
def test_batch_multicolumn(website):
    """
    Submit the batch syntax checker with a multiple-colums input file.

    This by the way also tests for the correct order of batch results.
    """
    variants = [('AB026906.1(SDHD):g.7872G>T', 'NM_003002.1:c.3_4insG'),
                ('NM_003002.1:c.3_4insG', 'AB026906.1(SDHD):g.7872G>T'),
                ('AL449423.14(CDKN2A_v002):c.5_400del', 'AL449423.14(CDKN2A_v002):c.5_400del')]
    result = _batch(website,
                    'syntax-checker',
                    file='\n'.join(['\t'.join(r) for r in variants]),
                    size=len(variants) * 2,
                    header='Input\tStatus',
                    lines=len(variants))
    for line in result.splitlines()[1:]:
        assert len(line.split('\t')) == len(variants[0]) * 2


def test_download_py(website):
    """
    Download a Python example client for the web service.
    """
    r = website.get('/downloads/client-suds.py')
    assert 'text/plain' in r.headers['Content-Type']
    assert '#!/usr/bin/env python' in r.data


def test_download_rb(website):
    """
    Download a Ruby example client for the web service.
    """
    r = website.get('/downloads/client-savon.rb')
    assert 'text/plain' in r.headers['Content-Type']
    assert '#!/usr/bin/env ruby' in r.data


def test_download_cs(website):
    """
    Download a C# example client for the web service.
    """
    r = website.get('/downloads/client-mono.cs')
    assert 'text/plain' in r.headers['Content-Type']
    assert 'public static void Main(String [] args) {' in r.data


def test_download_php(website):
    """
    Download a PHP example client for the web service.
    """
    r = website.get('/downloads/client-php.php')
    assert 'text/plain' in r.headers['Content-Type']
    assert '<?php' in r.data


def test_downloads_batchtest(website):
    """
    Download the batch test example file.
    """
    r = website.get('/downloads/batchtestnew.txt')
    assert 'text/plain' in r.headers['Content-Type']
    assert 'NM_003002.1:c.3_4insG' in r.data


def test_annotated_soap_api(website):
    """
    Test the SOAP documentation generated from the WSDL.
    """
    r = website.get('/soap-api')
    assert 'text/html' in r.headers['Content-Type']
    assert 'Web Service: Mutalyzer' in r.data


@with_references('NG_012337.1')
def test_getgs(website):
    """
    Test the /getGS interface used by LOVD2.
    """
    r = website.get('/getGS',
                    query_string={'variantRecord': 'NM_003002.2',
                                  'forward': '1',
                                  'mutationName': 'NG_012337.1:g.7055C>T'},
                    follow_redirects=True)
    assert '0 Errors' in r.data
    assert '0 Warnings' in r.data
    assert 'Raw variant 1: substitution at 7055' in r.data
    assert 'go to bottom' not in r.data
    assert '<input' not in r.data


@with_references('NG_012337.1')
def test_getgs_coding_multiple_transcripts(website):
    """
    Test the /getGS interface on a coding description and genomic
    reference with multiple transcripts.
    """
    r = website.get('/getGS',
                    query_string={'variantRecord': 'NM_003002.2',
                                  'forward': '1',
                                  'mutationName': 'NG_012337.1:c.45A>T'},
                    follow_redirects=False)
    assert '/name-checker?' in r.location
    assert 'description=NG_012337.1' in r.location


@with_references('NG_008939.1')
def test_getgs_variant_error(website):
    """
    Test the /getGS interface on a variant description with an error.
    """
    # The error is that position c.45 is a C, not an A.
    r = website.get('/getGS',
                    query_string={'variantRecord': 'NM_000532.4',
                                  'forward': '1',
                                  'mutationName': 'NG_008939.1:c.45A>T'},
                    follow_redirects=False)
    assert '/name-checker?' in r.location
    assert 'description=NG_008939.1' in r.location


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_variantinfo_g2c(website):
    """
    Test the /Variant_info interface used by LOVD2 (g to c).
    """
    r = website.get('/Variant_info',
                    query_string={'LOVD_ver': '2.0-29',
                                  'build': 'hg19',
                                  'acc': 'NM_203473.1',
                                  'var': 'g.48374289_48374389del'})
    assert 'text/plain' in r.headers['Content-Type']
    expected = '\n'.join(['1020', '0', '1072', '48', '48374289', '48374389', 'del'])
    assert r.data == expected


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_variantinfo_c2g(website):
    """
    Test the /Variant_info interface used by LOVD2 (c to g).
    """
    r = website.get('/Variant_info',
                    query_string={'LOVD_ver': '2.0-29',
                                  'build': 'hg19',
                                  'acc': 'NM_203473.1',
                                  'var': 'c.1020_1072+48del'})
    assert 'text/plain' in r.headers['Content-Type']
    expected = '\n'.join(['1020', '0', '1072', '48', '48374289', '48374389', 'del'])
    assert r.data == expected


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_variantinfo_c2g_downstream(website):
    """
    Test the /Variant_info interface used by LOVD2 (c variant downstream
    notation to g).
    """
    r = website.get('/Variant_info',
                    query_string={'LOVD_ver': '2.0-29',
                                  'build': 'hg19',
                                  'acc': 'NM_203473.1',
                                  'var': 'c.1709+d187del'})
    assert 'text/plain' in r.headers['Content-Type']
    expected = '\n'.join(['1709', '187', '1709', '187', '48379389', '48379389', 'del'])
    assert r.data == expected


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_variantinfo_no_variant(website):
    """
    Test the /Variant_info interface used by LOVD2 (without variant).
    """
    r = website.get('/Variant_info',
                    query_string={'LOVD_ver': '2.0-29',
                                  'build': 'hg19',
                                  'acc': 'NM_203473.1'})
    assert 'text/plain' in r.headers['Content-Type']
    assert 'text/plain' in r.content_type
    expected = '\n'.join(['-158', '1709', '1371'])
    assert r.data == expected


@pytest.mark.usefixtures('hg19_transcript_mappings')
def test_variantinfo_ivs(website):
    """
    Test the /Variant_info interface used by LOVD2 (with IVS positioning).
    """
    r = website.get('/Variant_info',
                    query_string={'LOVD_ver': '2.0-33',
                                  'build': 'hg19',
                                  'acc': 'NM_000249.3',
                                  'var': 'c.IVS10+3A>G'})
    assert 'text/plain' in r.headers['Content-Type']
    expected = '\n'.join(['884', '3', '884', '3', '37059093', '37059093', 'subst'])
    assert r.data == expected


@pytest.mark.usefixtures('db')
def test_upload_local_file(website):
    """
    Test the genbank uploader.
    """
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'data',
                        'AB026906.1.gb.bz2')
    r = website.post('/reference-loader',
                     data={'method': 'upload_method',
                           'file': (bz2.BZ2File(path), 'AB026906.1.gb')})
    assert 'Your reference sequence was loaded successfully.' in r.data

    dom = lxml.html.fromstring(r.data)
    reference_url = dom.cssselect('#reference_download')[0].attrib['href']

    r = website.get(reference_url)
    assert r.data == bz2.BZ2File(path).read()


@pytest.mark.usefixtures('db')
def test_upload_local_file_invalid(website):
    """
    Test the genbank uploader with a non-genbank file.
    """
    r = website.post('/reference-loader',
                     data={'method': 'upload_method',
                           'file': (BytesIO('this is not a genbank file'.encode('utf-8')), 'AB026906.1.gb')})
    assert 'Your reference sequence was loaded successfully.' not in r.data
    assert 'The file could not be parsed.' in r.data


@with_references('NM_002001.2')
def test_reference(website):
    """
    Test if reference files are cached.
    """
    r = website.get('/name-checker',
                    query_string={'description': 'NM_002001.2:g.1del'})
    assert '0 Errors' in r.data

    r = website.get('/reference/NM_002001.2.gb')
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'data',
                        'NM_002001.2.gb.bz2')
    assert r.data == bz2.BZ2File(path).read()


@with_references('NM_002001.2')
def test_reference_head(website):
    """
    Test if reference files are cached, by issuing a HEAD request.
    """
    r = website.get('/name-checker',
                    query_string={'description': 'NM_002001.2:g.1del'})
    assert '0 Errors' in r.data

    r = website.head('/reference/NM_002001.2.gb')
    assert r.status_code == 200


@pytest.mark.usefixtures('db')
def test_reference_head_none(website):
    """
    Test if non-existing reference files gives a 404 on a HEAD request.
    """
    r = website.head('/reference/NM_002001.2.gb')
    assert r.status_code == 404


@pytest.mark.usefixtures('hg19_transcript_mappings')
@with_references('NM_003002.2')
def test_bed(website):
    """
    BED track for variant.
    """
    r = website.get('/bed',
                    query_string={'description': 'NM_003002.2:c.274G>T'})
    assert 'text/plain' in r.headers['Content-Type']
    assert '\t'.join(['chr11', '111959694', '111959695', '274G>T', '0', '+']) in r.data


@pytest.mark.usefixtures('hg19_transcript_mappings')
@with_references('NM_000132.3')
def test_bed_reverse(website):
    """
    BED track for variant on reverse strand.
    """
    r = website.get('/bed',
                    query_string={'description': 'NM_000132.3:c.[4374A>T;4380_4381del]'})
    assert 'text/plain' in r.headers['Content-Type']
    assert '\t'.join(['chrX', '154157690', '154157691', '4374A>T', '0', '-']) in r.data
    assert '\t'.join(['chrX', '154157683', '154157685', '4380_4381del', '0', '-']) in r.data


def test_checksyntax_unicode(website):
    """
    Run check syntax form with an invalid variant description containing
    non-ASCII unicode characters.
    """
    r = website.get('/syntax-checker',
                    query_string={'description': 'La Pe\xf1a'})
    body = r.get_data(as_text=True)
    assert 'Fatal' in body
    assert 'The &quot;^&quot; indicates the position where the error occurred' in body
    assert 'Expected W:(0123...) (at char 2), (line:1, col:3)' in body


@pytest.mark.usefixtures('db')
def test_batch_unicode(website):
    """
    Submit a batch form with non-ASCII unicode characters in the input
    file.
    """
    file = '\n'.join(['\u2026AB026906.1:c.274G>T',
                      '\u2026AL449423.14(CDKN2A_v002):c.5_400del'])
    expected = [['\u2026AB026906.1:c.274G>T',
                 '(grammar): Expected W:(0123...) (at char 0), (line:1, col:1)'],
                ['\u2026AL449423.14(CDKN2A_v002):c.5_400del',
                 '(grammar): Expected W:(0123...) (at char 0), (line:1, col:1)']]

    data = {'job_type': 'syntax-checker',
            'email': 'test@test.test',
            'file': (BytesIO(file.encode('utf-8')), 'test.txt')}

    r = website.post('/batch-jobs',
                     data=data)
    progress_url = '/' + r.location.split('/')[-1]

    assert BatchJob.query.first().email == 'test@test.test'

    scheduler = Scheduler.Scheduler()
    scheduler.process()

    r = website.get(progress_url)

    dom = lxml.html.fromstring(r.data)
    result_url = dom.cssselect('#ifnot_items_left a')[0].attrib['href']

    r = website.get(result_url)
    assert 'text/plain' in r.headers['Content-Type']

    result = r.get_data(as_text=True).strip().split('\n')[1:]
    assert expected == [line.split('\t') for line in result]


@pytest.mark.usefixtures('db')
def test_batch_unicode_email(website):
    """
    Submit a batch form with non-ASCII unicode characters in the email
    address.
    """
    file = '\n'.join(['AB026906.1:c.274G>T',
                      'AL449423.14(CDKN2A_v002):c.5_400del'])
    expected = [['AB026906.1:c.274G>T',
                 'OK'],
                ['AL449423.14(CDKN2A_v002):c.5_400del',
                 'OK']]

    data = {'job_type': 'syntax-checker',
            'email': 'pe\xf1a@test.test',
            'file': (BytesIO(file.encode('utf-8')), 'test.txt')}

    r = website.post('/batch-jobs',
                     data=data)
    progress_url = '/' + r.location.split('/')[-1]

    assert BatchJob.query.first().email == 'pe\xf1a@test.test'

    scheduler = Scheduler.Scheduler()
    scheduler.process()

    r = website.get(progress_url)

    dom = lxml.html.fromstring(r.data)
    result_url = dom.cssselect('#ifnot_items_left a')[0].attrib['href']

    r = website.get(result_url)
    assert 'text/plain' in r.headers['Content-Type']

    result = r.get_data(as_text=True).strip().split('\n')[1:]
    assert expected == [line.split('\t') for line in result]
