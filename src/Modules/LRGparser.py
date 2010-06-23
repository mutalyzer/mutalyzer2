from Modules import GenRecord
import xml.dom.minidom

#TODO: Add verification: e.g. fixed annotation section, transcript etc.
#TODO: Add additional genes from updatable section

def _getContent(data, refname):
    return data.getElementsByTagName(refname)[0].lastChild.data.encode("utf8")

def createLrgRecord(data):
    record = GenRecord.Record()
    data = xml.dom.minidom.parseString(data)
    fixed = data.getElementsByTagName("fixed_annotation")[0]

    #LRG file should have dna as mol_type
    assert(_getContent(fixed, "mol_type") == "dna")
    record.mol_type = 'g'

    genename = _getContent(data, "lrg_gene_name")
    gene = GenRecord.Gene(genename)
    refseq = _getContent(fixed, "sequence")

    for tData in fixed.getElementsByTagName("transcript"):
        transcription = GenRecord.Locus(tData.getAttribute("name").encode("utf8"))
        #NOTE: Location is a property of the transcription, gene has no location
        transcription.location = \
            [int(tData.getAttribute("start")), int(tData.getAttribute("end"))]

        #get Exons
        exonPList = GenRecord.PList()
        for exon in tData.getElementsByTagName("exon"):
            co = exon.getElementsByTagName("lrg_coords")[0]
            exonPList.positionList.extend(\
                [int(co.getAttribute("start")), int(co.getAttribute("end"))])

        #get CDS
        CDSPList = GenRecord.PList()
        for CDS in tData.getElementsByTagName("coding_region"):
            CDSPList.positionList.extend(\
            [int(CDS.getAttribute("start")), int(CDS.getAttribute("end"))]

        if CDSPList.positionList:
            transcription.molType = 'c'
            CDSPList.location = [CDSPList.positionList[0], CDSPList.positionList[-1]]
            if len(CDSPList.positionList) == 2: # if only 1 CDS is found remove and
                CDSPList.positionList = []      # construnct positionList from exons
        else:
            transcription.molType = 'n'

        transcription.exon = exonPList
        transcription.CDS = CDSPList
        gene.transcriptList.append(transcription)

    #only 1 gene
    record.geneList.append(gene)

    #FIXME add a BOGUS id as GenBank ID
    record.id = "LRG_TMP"
    #FIXME add the ref_seq
    record.seq = refseq

    return record

#interesting fields in fixed annotation:
# id
# organism % taxon
# source => name
# mol_type
# sequence
# tanscript % name %start %end
#  cdna
#   sequence
#  coding_region %start %end
#   translation
#    sequence
#  exon
#   lrg_coords % start %end
#  intron %phase


if __name__ == "__main__":
    filename = "LRG_11.xml"
    f = open(filename)
    record = createLrgRecord(f.read())
