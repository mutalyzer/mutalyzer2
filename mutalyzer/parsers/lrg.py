"""
This module contains one public function create_record, which returns a
mutalyzer GenRecord.Record populated with data from an LRG file.

An LRG file is an XML formatted file and consists of a fixed and
updatable section. The fixed section contains a DNA sequence
and for that sequence a number of transcripts.

The updatable region could contain all sorts of annotation for the
sequence and transcripts. It can also contain additional (partial)
transcripts and mapping information.

More information on LRG files:
    https://www.lrg-sequence.org/faq/
    http://ftp.ebi.ac.uk/pub/databases/lrgex/docs/LRG_XML_schema_documentation_1_9.pdf
    http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG.rnc
    http://ftp.ebi.ac.uk/pub/databases/lrgex/docs/LRG.pdf

This module is based on the result of the minidom xml parser.

NOTE: A strong alternative to the minidom parser would be lxml, which is
already employed by Mutalyzer in other circumstances.
"""


from __future__ import unicode_literals

import xml.dom.minidom
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from mutalyzer import GenRecord


def _get_content(data, refname):
    """
    Return string-content of an XML textnode.

    @arg data:     a minidom object
    @type data:    object
    @arg refname:  the name of a member of the minidom object
    @type refname: unicode

    @return: The content of the textnode or an emtpy string
    @rtype: string
    """
    temp = data.getElementsByTagName(refname)
    if temp:
        return temp[0].lastChild.data
    return ""
#_get_content


def _attr2dict(attr):
    """
    Create a dictionary from the attributes of an XML node

    @arg attr: a minidom node
    @type attr: object

    @return: A dictionary with pairing of node-attribute names and values.
    Integer string values are converted to integers.
    @rtype: dictionary
    """
    ret = {}
    for key, value in attr.items():
        if value.isdigit():
            value = int(value)
        ret[key] = value
    return ret
#_attr2dict


def _get_coordinates(data, system=None):
    """
    Get attributes from descendent <coordinates> element as a dictionary. If
    more than one <coordinates> element is found, we have a preference for the
    one with 'coord_system' attribute equal to the `system` argument, if
    defined.
    """
    result = None
    coordinates = data.getElementsByTagName('coordinates')
    for coordinate in coordinates:
        attributes = _attr2dict(coordinate.attributes)
        if result and system and attributes.get('coord_system') != system:
            continue
        result = attributes
    return result
#_get_coordinates


def _get_gene_name(section):
    """
    Extract the gene name from the LRG record updatable section.

    NOTE: It is necessary to use the updatable section since there is no
    other way to identify the main gene directly from the LRG file.
    A possibility would be to use the HGNC id with some external service.
    Another way would be to make use of the special file with genes to LRG:
    http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_transcripts_GRCh38.txt

    :param section: (updatable) section of the LRG file
    :return: gene name present under the lrg annotation set
    """
    gene_name = ""
    annotation_nodes = section.getElementsByTagName("annotation_set")
    for anno in annotation_nodes:
        if anno.getAttribute("type") == "lrg":
            gene_name = _get_content(anno, "lrg_locus")
    return gene_name
#_get_gene_name


def _get_transcripts(section):
    """
    Extracts the transcripts present in the (fixed) section of the LRG file.

    :param section: (fixed) section of the LRG file
    :return: list of transcripts (GenRecord.Locus)
    """
    lrg_id = _get_content(section, 'id')

    transcripts = []
    for tdata in section.getElementsByTagName("transcript"):
        # iterate over the transcripts in the fixed section.
        # get the transcript from the updatable section and combine results
        transcript_name = tdata.getAttribute("name")[1:]
        transcription = GenRecord.Locus(transcript_name)

        coordinates = tdata.getElementsByTagName('coordinates')[0]

        # Set the locusTag, linkMethod (used in the output) and the location
        # LRG file transcripts can (for now) always be linked via the locustag
        transcription.locusTag = transcript_name and "t" + transcript_name
        transcription.linkMethod = "Locus Tag"
        transcription.location = [int(coordinates.getAttribute("start")),
                                  int(coordinates.getAttribute("end"))]

        # Get the transcript exons and store them in a position list.
        exonPList = GenRecord.PList()
        for exon in tdata.getElementsByTagName("exon"):
            coordinates = _get_coordinates(exon, lrg_id)
            exonPList.positionList.extend([int(coordinates["start"]),
                                           int(coordinates["end"])])
        exonPList.positionList.sort()

        # Get the CDS of the transcript and store them in a position list.
        # NOTE: up until now all CDSlists only consisted of a starting end
        # ending position, keep the possibility in mind that multiple CDS
        # regions are given
        CDSPList = GenRecord.PList()
        for cds_id, CDS in enumerate(tdata.getElementsByTagName("coding_region")):
            if cds_id > 0:
                # Todo: For now, we only support one CDS per transcript and
                #   ignore all others.
                #   By the way, I don't think the loop and sorting of CDS
                #   positions makes any sense here, but I leave it in place
                #   and just ignore everything except the first iteration.
                #translationName = CDS.getElementsByTagName("translation")[0].getAttribute("name").encode("utf8")[1:]
                #print 'Ignoring transcript %s translation %s' % (transcript_name, translationName)
                continue
            coordinates = _get_coordinates(CDS, lrg_id)
            CDSPList.positionList.extend([int(coordinates["start"]),
                                          int(coordinates["end"])])
        CDSPList.positionList.sort()

        # If there is a CDS position List set the transcriptflag to True
        if CDSPList.positionList:
            transcription.molType = 'c'
            CDSPList.location = [CDSPList.positionList[0],
                                 CDSPList.positionList[-1]]
            # If we only got the flanking CDS positions, we clear it and let
            # GenRecord.checkRecord reconstruct the correct CDS list
            # from the mRNA list later on
            if len(CDSPList.positionList) == 2:
                CDSPList.positionList = []
            transcription.translate = True
        else:
            transcription.molType = 'n'

        # Note: Not all the transcripts contain a coding_region.
        if tdata.getElementsByTagName('coding_region'):
            transcription.transcribe = True
            # Store CDS position lists in the transcription
            transcription.CDS = CDSPList

        # Store exon position list in the transcription
        transcription.exon = exonPList

        # Add the transcription to the transcripts list
        transcripts.append(transcription)
    return transcripts
#_get_transcripts


def create_record(data):
    """
    Create a GenRecord.Record of a LRG <xml> formatted string.

    @arg data: Content of LRG file
    @type data: byte string

    @return: GenRecord.Record instance
    @rtype: object
    """
    # Initiate the GenRecord.Record
    record = GenRecord.Record()
    record._sourcetype = "LRG"

    # Parse the xml data string and extract the fixed and updatable section
    data = xml.dom.minidom.parseString(data)
    fixed = data.getElementsByTagName("fixed_annotation")[0]
    updatable = data.getElementsByTagName("updatable_annotation")[0]

    # Get the organism
    record.organism = _get_content(data, 'organism')

    # Get the sequence from the fixed section
    record.molType = 'g'
    record.seq = Seq(_get_content(fixed, "sequence"), IUPAC.unambiguous_dna)

    # Get the gene name
    gene = GenRecord.Gene(_get_gene_name(updatable))
    # Add transcripts information from the fixed section to the main gene.
    gene.transcriptList = _get_transcripts(fixed)
    record.geneList = [gene]

    return record
#create_record
