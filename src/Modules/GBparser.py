from Bio import SeqIO  # read()
import bz2             # BZ2Compressor(), BZ2File()
from GenRecord import PList, Locus, Gene, Record, GenRecord

"""
    Module contains one public function createGBRecord which returns a
    mutalyzer GenRecord.Record populated with data from a GenBank file.
"""

def __location2pos(location) :
    """
        Convert a location object to a tuple of integers.

        Arguments:
            location ; A location object (see the BioPython documentation).

        Returns:
            List ; A tuple of integers.
    """

    ret = []

    if not str(location.start).isdigit() or \
       not str(location.end).isdigit() :
        return None
    #if

    ret.append(location.start.position + 1)
    ret.append(location.end.position)

    return ret
#__location2pos

def __locationList2posList(locationList) :
    """
        Convert a list of locations to a list of integers.

        Arguments:
            locationList ; A list of locations
                            (see the BioPython documentation).

        Returns:
            List ; A list (of even length) of integers.
    """

    ret = []

    if not str(locationList.location.start).isdigit() or \
       not str(locationList.location.end).isdigit() :
        return None
    #if

    for i in locationList.sub_features :
        if i.ref : # This is a workaround for a bug in BioPython.
            ret = None
            break
        #if
        temp = __location2pos(i.location)
        ret.append(temp[0])
        ret.append(temp[1])
    #for

    return ret
#__locationList2posList

def createGBRecord(filename):
    """
        Create a GenRecord.Record from a GenBank file

        Arguments:
            filename    ; The full path to the compresed GenBank file

        Returns:
            record      ; A GenRecord.Record instance
    """

    # first create an intermediate genbank record with BioPython
    file_handle = bz2.BZ2File(filename, "r")
    biorecord = SeqIO.read(file_handle, "genbank")
    file_handle.close()

    record = Record()
    record.seq = biorecord.seq

    mRNAProducts = []
    CDSProducts = []
    for i in biorecord.features :
        if i.qualifiers :
            if i.qualifiers.has_key("gene") :
                if i.type == "mRNA" :
                    if i.qualifiers.has_key("product") :
                        mRNAProducts.append(i.qualifiers["product"][0])
                if i.type == "CDS" :
                    if i.qualifiers.has_key("product") :
                        CDSProducts.append(i.qualifiers["product"][0])
            #if

    for i in biorecord.features :
        if i.qualifiers :
            if i.type == "source" :
                if i.qualifiers.has_key("mol_type") :
                    if i.qualifiers["mol_type"][0] == "mRNA" :
                        record.molType = 'n'
                    else :
                        record.molType = 'g'
                #if
                if i.qualifiers.has_key("organelle") :
                    record.organelle = i.qualifiers["organelle"][0]
                    if record.organelle == "mitochondrion" :
                        record.molType = 'm'
                #if

                fakeGene = Locus("001")
                record.source.transcriptList.append(fakeGene)
                fakeGene.CDS = PList()
                fakeGene.CDS.location = __location2pos(i.location)
            #if

            if i.qualifiers.has_key("gene") :
                gene = i.qualifiers["gene"][0]

                GeneInstance = record.findGene(gene)
                if not GeneInstance :
                    GeneInstance = Gene(gene)
                    record.geneList.append(GeneInstance)
                #if

                if i.type == "gene" :
                    if i.strand :
                        GeneInstance.orientation = i.strand
                    GeneInstance.location = __location2pos(i.location)
                #if

                # Look if there is a locus tag present, if not, give it the
                # default tag `001'.
                locusName = "001"
                locusTag = None
                if i.qualifiers.has_key("locus_tag") :
                    locusTag = i.qualifiers["locus_tag"][0]
                    locusName = locusTag[-3:]
                #if

                LocusInstance = GeneInstance.findLocus(locusName)
                if not LocusInstance :
                    LocusInstance = Locus(locusName)
                    GeneInstance.transcriptList.append(LocusInstance)
                #if

                if i.type == "mRNA" :
                    PListInstance = PList()
                    LocusInstance.mRNA = PListInstance

                    posList = __locationList2posList(i)
                    if posList != None :
                        PListInstance.location = \
                            __location2pos(i.location)
                        PListInstance.positionList = posList
                    #if
                    if i.qualifiers.has_key("transcript_id") :
                        LocusInstance.transcriptID = \
                            i.qualifiers["transcript_id"][0]
                    LocusInstance.locusTag = locusTag
                #if
                if i.type == "CDS" :
                    PListInstance = PList()
                    LocusInstance.CDS = PListInstance

                    PListInstance.location = __location2pos(i.location)
                    PListInstance.positionList = \
                        __locationList2posList(i)

                    if i.qualifiers.has_key("transl_table") :
                        LocusInstance.txTable = \
                            int(i.qualifiers["transl_table"][0])
                    if i.qualifiers.has_key("protein_id") :
                        LocusInstance.proteinID = \
                            i.qualifiers["protein_id"][0]
                    LocusInstance.locusTag = locusTag
                #if
                if i.type == "exon" :
                    if not LocusInstance.exon :
                        LocusInstance.exon = PList()
                    LocusInstance.exon.positionList.extend(
                        __location2pos(i.location))
                #if
            #if                            
        #if
    #for
    return record
#parseRecord

