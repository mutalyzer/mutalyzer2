from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Modules import GenRecord
import xml.dom.minidom
DEBUG = False

def __debugParsedData(title, data):
    if not DEBUG: return
    import pprint
    print "#"*79+"\nDEBUG: Start of "+title+"\n"+"#"*79
    pprint.pprint(data)
    print "#"*79+"\nDEBUG: End of "+title+"\n"+"#"*79

def _getContent(data, refname):
    """Return string-content of an XML textnode"""
    temp = data.getElementsByTagName(refname)
    if temp:
        return temp[0].lastChild.data.encode("utf8")
    else:
        return ""

def _attr2dict(attr):
    """Create a dictionary from the attributes of an XML node"""
    ret = {}
    for key, value in attr.items():
        value = value.isdigit() and int(value) or value.encode("utf-8")
        ret[key.encode("utf-8")] = value
    return ret



def createLrgRecord(data):
    """
        Create a GenRecord.Record of a LRG <xml> formatted string

        Input:
            data    :   String content of LRG file

        Output
            record  :   GenRecord.Record instance


    """
    record = GenRecord.Record()
    record._sourcetype = "LRG"

    data = xml.dom.minidom.parseString(data)
    fixed = data.getElementsByTagName("fixed_annotation")[0]
    updatable = data.getElementsByTagName("updatable_annotation")[0]

    # Get all raw information from updatable section into a nested dict
    updParsed = parseUpdatable(updatable)
    __debugParsedData("Updatable Section",updParsed)

    # Updatable Section get mapping and geneList
    record.mapping = getMapping(updParsed["LRG"]["mapping"])
    record.geneList = genesFromUpdatable(updParsed)

    # Fixed Section get sequence
    assert(_getContent(fixed, "mol_type") == "dna")
    record.molType = 'g'
    record.seq = Seq(_getContent(fixed, "sequence"), IUPAC.unambiguous_dna)

    # Get the genename of the fixed gene in the LRG and put that gene on top of
    # the geneList. This is necessary because there is no way to otherwise
    # identify the main gene in the LRG file.
    geneName = updParsed["LRG"]["genename"]
    gene = [gene for gene in record.geneList if gene.name == geneName][0]
    record.geneList.remove(gene)
    record.geneList.insert(0, gene)

    # Update the Gene object with transcripts from fixed section
    for tData in fixed.getElementsByTagName("transcript"):
        transcriptName = tData.getAttribute("name").encode("utf8")[1:]
        transcription = [t for t in gene.transcriptList if t.name ==
                transcriptName][0]

        transcription.location = \
          [int(tData.getAttribute("start")), int(tData.getAttribute("end"))]

        #get Exons
        exonPList = GenRecord.PList()
        for exon in tData.getElementsByTagName("exon"):
            co = exon.getElementsByTagName("lrg_coords")[0]
            exonPList.positionList.extend(\
              [int(co.getAttribute("start")), int(co.getAttribute("end"))])
        exonPList.positionList.sort()

        # get CDS, keep the possibility in mind 
        # that multiple CDS regions are given
        CDSPList = GenRecord.PList()
        for CDS in tData.getElementsByTagName("coding_region"):
            CDSPList.positionList.extend(\
            [int(CDS.getAttribute("start")), int(CDS.getAttribute("end"))])
        CDSPList.positionList.sort()

        if CDSPList.positionList:
            transcription.molType = 'c'
            CDSPList.location = [CDSPList.positionList[0],
                                 CDSPList.positionList[-1]]
            # if only 1 CDS is found clear the positionList 
            # GenRecord.checkRecord will reconstruct it later on
            if len(CDSPList.positionList) == 2:
                CDSPList.positionList = []
            transcription.translate = True
        else:
            transcription.molType = 'n'

        transcription.transcribe = True

        transcription.exon = exonPList
        transcription.CDS = CDSPList


    #for
    return record
#createLrgRecord

def genesFromUpdatable(updParsed):
    """Get the genes from the updatable section"""
    genes = []
    for geneSymbol, geneData in updParsed["NCBI"].items():
        gene = GenRecord.Gene(geneSymbol)
        gene.location = [geneData["geneAttr"]["start"],
                         geneData["geneAttr"]["end"]]
        gene.longName = geneData["geneLongName"]
        gene.orientation = int(geneData["geneAttr"]["strand"])
        # Get transcripts
        transcripts = transcriptsFromParsed(geneData["transcripts"])
        if not transcripts:
            transcripts = _emptyTranscripts(geneSymbol, geneData)
        gene.transcriptList = transcripts
        genes.append(gene)
    #for
    for geneSymbol, geneData in updParsed["Ensembl"].items():
        #add notation from gathered Ensembl annotation TODO
        pass
    return genes
#genesFromUpdatable

def transcriptsFromParsed(parsedData):
    """Get the transcripts for each gene from parsed Data"""
    #parsedData contains the pre-parsed transcripts in dict format

    transcripts = []

    # Store transcripts without a FixedID
    nofixed = parsedData.pop("noFixedId")

    # First the fixed transcripts
    for trName, trData in parsedData.items():
        transcripts.append(_transcriptPopulator(trName, trData))

    # Second the nonfixed transcripts
    # FIXME: How to name these transcripts? 
    for trData in nofixed:
        transcripts.append(_transcriptPopulator("", trData))

    return transcripts
#transcriptsFromParsed

def _emptyTranscripts(symbol, data):
    transcript = GenRecord.Locus('')
    transcript.molType = 'n'
    mRNA = GenRecord.PList()
    location = [data["geneAttr"]["start"],
            data["geneAttr"]["end"]]
    mRNA.location = location
    transcript.mRNA = mRNA
    return [transcript,]

def _transcriptPopulator(trName, trData):
    transcript = GenRecord.Locus(trName)
    transcript.longName = trData.get("transLongName")
    if trData.has_key("transAttr"):
        tA = trData["transAttr"]
        transcript.transcriptID = tA.get("transcript_id")
        transcript.location = [tA.get("start"), tA.get("end")]

    if trData.has_key("proteinAttr"):
        transcript.protLongName = trData.get("proteinLongName")

        pA = trData["proteinAttr"]
        transcript.proteinID = pA.get("accession")
        CDS = GenRecord.PList()
        CDS.location = [pA.get("cds_start"), pA.get("cds_end")]
        if not any(CDS.location): CDS = None
        transcript.CDS = CDS
    else:
        transcript.molType = 'n'

    # Check if the transcript has a name, if not; use transcriptid
    # This is needed for the transcripts without a fixed ID
    #if trName == "":
    #    transcript.name = '1' #transcript.transcriptID

    return transcript
#_transcriptPopulator

def getMapping(rawMapData):
    mapp, span, diffs = rawMapData
    ret = { "assembly":     mapp.get("assembly"),
            "chr_name":     mapp.get("chr_name"),
            "chr_id":       mapp.get("chr_id"),
            "chr_location": [int(span.get("start")),
                             int(span.get("end"))],
            "strand":       int(span.get("strand")),
            "lrg_location": [int(span.get("lrg_start")),
                             int(span.get("lrg_end"))],
            "diffs":        diffs}
    return ret
#getMapping

def parseUpdatable(data):
    ret = {"LRG":{}, "NCBI":{}, "Ensembl":{}}
    annotation_nodes = data.getElementsByTagName("annotation_set")
    for anno in annotation_nodes:
        name = _getContent(anno, "name")
        if name == "LRG":
            ret["LRG"] = getLrgAnnotation(anno)
        elif name == "NCBI RefSeqGene":
            ret["NCBI"] = getFeaturesAnnotation(anno)
        elif name == "Ensembl":
            ret["Ensembl"] = getFeaturesAnnotation(anno)
        else:
            #This annotation node is not yet implemented
            pass
    #for
    return ret
#parseUpdatable

def getLrgAnnotation(data):
    ret = {"mapping": (), "genename":""}

    # Get the mapping
    for mapp in data.getElementsByTagName("mapping"):
        mapattr = _attr2dict(mapp.attributes)

        # only the most recent mapping
        if not(mapattr.has_key("most_recent")): continue

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
    ret["genename"] = _getContent(data, "lrg_gene_name")

    return ret
#getLrgAnnotation

def getFeaturesAnnotation(data):
    ret = {} # Get annotation per gene symbol {"COL1A1":{}}

    #Check if features exists
    if not data.getElementsByTagName("features"): return ret
    feature = data.getElementsByTagName("features")[0]

    for gene in feature.getElementsByTagName("gene"):
        geneAttr = _attr2dict(gene.attributes)
        geneLongName = _getContent(gene, "long_name")
        transcripts = {"noFixedId": []}
        for transcript in gene.getElementsByTagName("transcript"):
            transAttr = _attr2dict(transcript.attributes)
            transLongName = _getContent(transcript, "long_name")

            # Check if the transcript has a protein product
            proteinProduct =\
                    transcript.getElementsByTagName("protein_product")
            if proteinProduct:
                protein = proteinProduct[0]
                proteinAttr = _attr2dict(protein.attributes)
                proteinLongName = _getContent(protein, "long_name")
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
        ret[geneAttr["symbol"]] = {\
                    "geneAttr":         geneAttr,
                    "geneLongName":     geneLongName,
                    "transcripts":      transcripts}
    #for gene
    return ret
#getFeaturesAnnotation

if __name__ == "__main__":
    print "Use the unit tests to test this Module"
