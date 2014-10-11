"""
Module contains one public function create_record which returns a
mutalyzer GenRecord.Record populated with data from a LRG file.

A LRG file is an XML formatted file and consists of a fixed and
updatable section. The fixed section contains a DNA sequence
and for that sequence a number of transcripts.

The updatable region could contain all sorts of annotation for the
sequence and transcripts. It can also contain additional (partial)
transcripts and mapping information.

This module is based on the result of the minidom xml parser.

Todo: Check document to Relax NG LRG schema.
    ftp://ftp.ebi.ac.uk/pub/databases/lrgex/

NOTE: A strong alternative to the minidom parser would be ElementTree which is
added in python2.5. Its main strengths are speed and readability [pythonesque].
(http://docs.python.org/library/xml.etree.elementtree.html)
"""


from __future__ import unicode_literals

import xml.dom.minidom
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from mutalyzer import GenRecord


def _debug_parsed_data(title, data):
    """
    Output additional data to stdout. Used for debugging the
    intermediate format used while parsing a LRG file.

    @requires: pprint

    @arg title:
    @type title: string
    @arg data: minidom object
    @type data: object
    """
    import pprint       #Only imported when the debug flag is set
    print "#"*79+"\nDEBUG: Start of "+title+"\n"+"#"*79
    pprint.pprint(data)
    print "#"*79+"\nDEBUG: End of "+title+"\n"+"#"*79
#_debug_parsed_data


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
    else:
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
    cos = data.getElementsByTagName('coordinates')
    for co in cos:
        attributes = _attr2dict(co.attributes)
        if result and system and attributes.get('coord_system') != system:
            continue
        result = attributes
    return result
#_get_coordinates


def create_record(data):
    """
    Create a GenRecord.Record of a LRG <xml> formatted string.

    @arg data: Content of LRG file
    @type data: string

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

    # Get LRG id
    lrg_id = _get_content(data, 'id')

    # Get organism
    record.organism = _get_content(data, 'organism')

    # Get all raw information from updatable section into a nested dict
    updParsed = parseUpdatable(updatable, lrg_id)

    # NOTE: To get insight in the structure of the intermediate
    # nested dictionary format please comment out the following line
    #_debug_parsed_data('Updatable Section', updParsed)

    # Get the genomic mapping from the Updatable Section -> LRG
    # NOTE: The mapping is not yet used in the mutalyzer program
    record.mapping = getMapping(updParsed["LRG"]["mapping"])

    # Get the geneList from the Updatable Section
    record.geneList = genesFromUpdatable(updParsed)

    # At this point the GenRecord.Record instance contains all information
    # from the updatable section.

    # get sequence from Fixed Section
    #assert(_get_content(fixed, "mol_type") == "dna")
    record.molType = 'g'
    record.seq = Seq(_get_content(fixed, "sequence"), IUPAC.unambiguous_dna)

    # Get the genename of the fixed gene in the LRG
    # and put that gene on top of the geneList.
    geneName = updParsed["LRG"]["genename"]
    gene = [gene for gene in record.geneList if gene.name == geneName][0]
    record.geneList.remove(gene)
    record.geneList.insert(0, gene)
    # NOTE: This is necessary because there is no way to otherwise
    # identify the main gene in the LRG file.

    # Update the Gene object in the Record
    # Add transcripts information from the fixed section to the main gene.
    for tData in fixed.getElementsByTagName("transcript"):
        # iterate over the transcripts in the fixed section.
        # get the transcript from the updatable section and combine results
        transcriptName = tData.getAttribute("name")[1:]
        transcription = [t for t in gene.transcriptList if t.name ==
                transcriptName][0]  #TODO?: swap with gene.findLocus

        coordinates = tData.getElementsByTagName('coordinates')[0]

        # Set the locusTag, linkMethod (used in the output) and the location
        # LRG file transcripts can (for now) always be linked via the locustag
        transcription.locusTag = transcriptName and "t"+transcriptName
        transcription.linkMethod = "Locus Tag"
        transcription.location = \
          [int(coordinates.getAttribute("start")), int(coordinates.getAttribute("end"))]

        #get the Exons of transcript and store them in a position list.
        exonPList = GenRecord.PList()
        for exon in tData.getElementsByTagName("exon"):
            coordinates = _get_coordinates(exon, lrg_id)
            exonPList.positionList.extend(
              [int(coordinates["start"]), int(coordinates["end"])])
        exonPList.positionList.sort()

        # get the CDS of the transcript and store them in a position list.
        # NOTE: up until now all CDSlists only consisted of a starting end
        # ending position, keep the possibility in mind that multiple CDS
        # regions are given
        CDSPList = GenRecord.PList()
        for cds_id, CDS in enumerate(tData.getElementsByTagName("coding_region")):
            if cds_id > 0:
                # Todo: For now, we only support one CDS per transcript and
                #   ignore all others.
                #   By the way, I don't think the loop and sorting of CDS
                #   positions makes any sense here, but I leave it in place
                #   and just ignore everything except the first iteration.
                #translationName = CDS.getElementsByTagName("translation")[0].getAttribute("name").encode("utf8")[1:]
                #print 'Ignoring transcript %s translation %s' % (transcriptName, translationName)
                continue
            coordinates = _get_coordinates(CDS, lrg_id)
            CDSPList.positionList.extend(\
            [int(coordinates["start"]), int(coordinates["end"])])
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

        # all transcripts in the fixed section are transcribable
        transcription.transcribe = True

        # Store the exon and CDS position lists in the transcription
        transcription.exon = exonPList
        transcription.CDS = CDSPList
    #for
    return record
#create_record


def genesFromUpdatable(updParsed):
    """
    Populate GenRecord.Gene instances with updatable LRG node data.

    @arg updParsed: Intermediate nested dict of updatable section
    @type updParsed: dictionary

    @return: genes ; List of GenRecord.Gene instances, populated with the
    content of the updatable section
    @rtype: list
    """
    genes = []
    for geneSymbol, geneData in updParsed["NCBI"].items():
        gene = GenRecord.Gene(geneSymbol)
        gene.location = [geneData["coordinates"]["start"],
                         geneData["coordinates"]["end"]]
        gene.longName = geneData["geneLongName"]
        gene.orientation = int(geneData["coordinates"]["strand"])
        # Get transcripts
        transcripts = transcriptsFromParsed(geneData["transcripts"])
        if not transcripts:
            transcripts = _emptyTranscripts(geneData)
            # Todo: For now we skip genes for which we have no transcripts
            #   since we don't know how to name them anyway.
            continue
        gene.transcriptList = transcripts
        genes.append(gene)
    #for
    for geneSymbol, geneData in updParsed["Ensembl"].items():
        #TODO: add annotation from Ensembl section, what should we add?
        pass
    return genes
#genesFromUpdatable


def transcriptsFromParsed(parsedData):
    """
    Populate GenRecord.Locus instances with updatable LRG node data

    @arg parsedData: Dictionary of transcript data. See getFeaturesAnnotation
    @type parsedData: dictionary

    @return: transcripts ; List of GenRecord.Locus instances, populated with the
    content of the parsed Data
    @rtype: list
    """
    transcripts = []

    # Temporary remove transcripts without a FixedID from the data
    nofixed = parsedData.pop("noFixedId")

    # First add the transcripts linked to the transcripts in the fixed section
    for trName, trData in parsedData.items():
        transcripts.append(_transcriptPopulator(trName, trData))

    # Second add the transcripts not linked
    # FIXME: How to name these transcripts?
    # Todo: For now we skip these transcripts since we don't know how to name
    #   them.
    #for trData in nofixed:
    #    transcripts.append(_transcriptPopulator("", trData))

    return transcripts
#transcriptsFromParsed


def _emptyTranscripts(data):
    """
    Populate a GenRecord.Locus instance with minimal data to make the
    gene compatible with mutalyzer. Data abstracted from the gene.

    @arg data: Data from the gene which is used to populate the create a minimal
    GenRecord.Locus instance
    @type data: dictionary

    @return: List with a single bogus GenRecord.Locus instance, in which
    location and mRNA are copied from the gene
    @rtype: list

    @todo: This function can be moved to the GenRecord.checkRecord method.
    """
    transcript = GenRecord.Locus('')
    transcript.molType = 'n'
    mRNA = GenRecord.PList()
    location = [data["coordinates"]["start"],
            data["coordinates"]["end"]]
    mRNA.location = location
    transcript.mRNA = mRNA
    return [transcript,]
#_emptyTranscripts


def _transcriptPopulator(trName, trData):
    """
    Populate GenRecord.Locus instance with updatable LRG node data.

    @arg trName: Name of the transcript
    @type trName: string
    @arg trData: Data associated with the transcript
    @type trData: dictionary

    @return: transcript  ; GenRecord.Locus instance, populated with the content
    of the parsed Data
    @rtype: object
    """
    transcript = GenRecord.Locus(trName)
    transcript.transcriptProduct = trData.get("transLongName")
    if trData.has_key("transAttr"):
        tA = trData["transAttr"]
        transcript.transcriptID = tA.get("accession")
        transcript.location = [tA.get("start"), tA.get("end")]

    if trData.has_key("proteinAttr"):
        transcript.proteinProduct = trData.get("proteinLongName")

        pA = trData["proteinAttr"]
        transcript.proteinID = pA.get("accession")
        CDS = GenRecord.PList()
        CDS.location = [pA.get("start"), pA.get("end")]
        # If the CDS list is empty set it to None. This will be fixed by the
        # GenRecord.checkRecord method
        if not any(CDS.location): CDS = None
        transcript.CDS = CDS
    else:
        transcript.molType = 'n'

    return transcript
#_transcriptPopulator


def getMapping(rawMapData):
    """
    Collect all necessary info to map the current LRG sequence to the
    genomic reference supplied by the file.

    @arg rawMapData: A list of dictionaries with the raw mapping info
    @type rawMapData: list

    @return: dictionary with the mapping info
    @rtype: dictionary
    """

    mapp, span, diffs = rawMapData
    ret = { "assembly":     mapp.get("coord_system"),      # Assembly Reference
            "chr_name":     mapp.get("other_name"),        # Chromosome name
            "chr_id":       mapp.get("other_id"),          # Chr. Reference NC
            "chr_location": [int(span.get("other_start")), # Start & End Position
                             int(span.get("other_end"))],  #    of seq on the ref
            "strand":       int(span.get("strand")),       # Forward:1 Reverse:-1
            "lrg_location": [int(span.get("lrg_start")),   # Start & End Position
                             int(span.get("lrg_end"))],    #    of seq on the LRG
            "diffs":        diffs}                         # Unformatted Diffs
    # NOTE: diffs contains the differences between the LRG and reference seq.
    return ret
#getMapping


def parseUpdatable(data, lrg_id):
    """
    Mediator function which transforms the minidom object to a nested dict
    and filters information needed to construct the GenRecord.Record.

    NOTE: an xml node has attributes and elements, this function squashes this
    ambiguity and collects only the attributes and elements of interest

    @arg data: The LRG file's updatable section node
    @type data: dictionary

    @return: Contains the fields of interest of the LRG NCBI and Ensembl
    sections of the updatable node
    @rtype: dictionary
    """
    ret = {"LRG":{}, "NCBI":{}, "Ensembl":{}}
    annotation_nodes = data.getElementsByTagName("annotation_set")
    for anno in annotation_nodes:
        name = _get_content(anno, "name")
        if name == "LRG":
            ret["LRG"] = getLrgAnnotation(anno)
        elif name == "NCBI RefSeqGene":
            ret["NCBI"] = getFeaturesAnnotation(anno, lrg_id)
        elif name == "Ensembl":
            ret["Ensembl"] = getFeaturesAnnotation(anno, lrg_id)
        else:
            #This annotation node is not yet implemented
            pass
    #for
    return ret
#parseUpdatable


def getLrgAnnotation(data):
    """
    Retrieves three parts of the LRG annotation:
        - the mapping of this LRG file to a genomic reference
        - a diference list between the LRG sequence and the ref seq
        - the genename of the main gene annotated by this LRG file

    @arg data: updatable section -> Annotations -> LRG node
    @type data:  dictionary

    @return: Contains the mapping [+ opt. diffs] and the genename
    @rtype: dictionary
    """
    ret = {"mapping": (), "genename":""}
    # Get the mapping
    for mapp in data.getElementsByTagName("mapping"):
        mapattr = _attr2dict(mapp.attributes)
        # only the most recent mapping
        if ret['mapping'] and not mapattr.has_key("most_recent"):
            continue
        # check if span exists
        for span in mapp.getElementsByTagName("mapping_span"):
            spanattr = _attr2dict(span.attributes)
            break
        else: continue
        diffs = []
        for diff in span.getElementsByTagName("diff"):
            diffs.append(_attr2dict(diff.attributes))
        ret["mapping"] = (mapattr,spanattr,diffs)
    #for
    # Get the LRG Gene Name, this is the main gene in this LRG
    ret["genename"] = _get_content(data, "lrg_locus")
    return ret
#getLrgAnnotation


def getFeaturesAnnotation(data, lrg_id):
    """
    Retrieves feature annotations from NCBI & Ensembl nodes.
        - List of genes
        - List of transcripts per gene
        - Potential Product of a transcript

    If a transcript can not be linked to a transcript from the fixed section it
    is stored in the noFixedId list.

    NOTE: an xml node has attributes and elements, this function squashes this
    ambiguity and collects only the attributes and elements of interest

    @todo: check documentation

    @arg data: updatable section -> Annotations -> NCBI | Ensembl
    @type data: dictionary

    @return: nested dict ; toplevel contains the genesymbols as keys e.g:
        - COL1A1 :
            - coordinates     : {}
            - geneLongName    : ""
            - transcripts     : {}
        - coordinates contains the start, end and strand info
        - transcripts contains a list of transcripts that could not be linked to
        the fixed section AND it contains each linked transcript with the
        locustag as key e.g:
            - 1 :
                - transAttr       : {}
                - transLongName   : ""
                - proteinAttr     : {}
                - proteinLongName : ""
            - transAttr & proteinAttr contain reference, start and end info
    @rtype: dictionary
    """
    ret = {} # Get annotation per gene symbol {"COL1A1":{}}
    #Check if features exists
    if not data.getElementsByTagName("features"): return ret
    feature = data.getElementsByTagName("features")[0]
    for gene in feature.getElementsByTagName("gene"):
        symbol = _get_content(gene, 'symbol')
        coordinates = _get_coordinates(gene, lrg_id)
        geneLongName = _get_content(gene, "long_name")
        transcripts = {"noFixedId": []}
        for transcript in gene.getElementsByTagName("transcript"):
            transAttr = _attr2dict(transcript.attributes)
            transAttr.update(_get_coordinates(transcript, lrg_id))
            transLongName = _get_content(transcript, "long_name")
            # Check if the transcript has a protein product
            proteinProduct =\
                    transcript.getElementsByTagName("protein_product")
            if proteinProduct:
                protein = proteinProduct[0]
                proteinAttr = _attr2dict(protein.attributes)
                proteinAttr.update(_get_coordinates(protein, lrg_id))
                proteinLongName = _get_content(protein, "long_name")
            else:
                proteinAttr = {}
                proteinLongName = ""

            #TODO: Add CCDS

            #transRet contains the fields to return for a transcript
            transRet = {\
                    "transAttr":        transAttr,
                    "transLongName":    transLongName,
                    "proteinAttr":      proteinAttr,
                    "proteinLongName":  proteinLongName}

            # Check if the transcript is linked to the fixed section
            if transAttr.has_key("fixed_id"):
                transcripts[transAttr["fixed_id"][1:]] = transRet
            else:
                transcripts["noFixedId"].append(transRet)

        #for transcript
        ret[symbol] = {\
                    "coordinates":      coordinates,
                    "geneLongName":     geneLongName,
                    "transcripts":      transcripts}
    #for gene
    return ret
#getFeaturesAnnotation
