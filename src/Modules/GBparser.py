from Bio import SeqIO, Entrez  # read()
import bz2                     # BZ2Compressor(), BZ2File()
from GenRecord import PList, Locus, Gene, Record, GenRecord
import Db

"""
    Module contains one public function createGBRecord which returns a
    mutalyzer GenRecord.Record populated with data from a GenBank file.
"""

class GBparser() :
    """
    """

    def __init__(self) :
        """
        """

        import Config
        config = Config.Config()
        Entrez.email = config.Retriever.email
        self.__database = Db.Cache(config.Db)
    #__init__

    def __location2pos(self, location) :
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
    
    def __locationList2posList(self, locationList) :
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
            temp = self.__location2pos(i.location)
            ret.append(temp[0])
            ret.append(temp[1])
        #for
    
        return ret
    #__locationList2posList
    
    def __transcriptToProtein(self, transcriptAcc) :
        """
        """
    
        proteinAcc = self.__database.getProtAcc(transcriptAcc)
        if not proteinAcc :
            handle = Entrez.esearch(db = "nucleotide", term = transcriptAcc)
            result = Entrez.read(handle)
            handle.close()
            
            transcriptGI = result["IdList"][0]
            
            handle = Entrez.elink(dbfrom = "nucleotide", db = "protein", 
                                  id = transcriptGI)
            result = Entrez.read(handle)
            handle.close()
            
            proteinGI = result[0]["LinkSetDb"][0]["Link"][0]["Id"]
            
            handle = Entrez.efetch(db = "protein", id = proteinGI, 
                                   rettype = "acc")

            proteinAcc = handle.read().split('.')[0]
            handle.close()

            self.__database.insertLink(transcriptAcc, proteinAcc)
        #if
    
        return proteinAcc
    #__transcriptToProtein
    
    def createGBRecord(self, filename):
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
    
        #mRNAProducts = []
        #CDSProducts = []
        #for i in biorecord.features :
        #    if i.qualifiers :
        #        if i.qualifiers.has_key("gene") :
        #            if i.type == "mRNA" :
        #                if i.qualifiers.has_key("product") :
        #                    mRNAProducts.append(i.qualifiers["product"][0])
        #            if i.type == "CDS" :
        #                if i.qualifiers.has_key("product") :
        #                    CDSProducts.append(i.qualifiers["product"][0])
        #        #if
        #print mRNAProducts
        #print CDSProducts
    
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
                    fakeGene.CDS.location = self.__location2pos(i.location)
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
                        GeneInstance.location = self.__location2pos(i.location)
                        if not GeneInstance.location :
                            GeneInstance.transcribe = False
                    #if

                    # RESOLV
                    LocusInstance = None
                    locusTag = None
                    if i.qualifiers.has_key("locus_tag") :
                        locusTag = i.qualifiers["locus_tag"][0]
                        locusName = locusTag[-3:]
                        LocusInstance = GeneInstance.findLocus(locusTag[-3:])
                    #if
                    else :
                        if i.qualifiers.has_key("transcript_id") :
                            LocusInstance = GeneInstance.findLink(
                                self.__transcriptToProtein(
                                i.qualifiers["transcript_id"][0].split('.')[0]))
                        if i.qualifiers.has_key("protein_id") :
                            LocusInstance = GeneInstance.findLink(
                                i.qualifiers["protein_id"][0].split('.')[0])
                    #else                                
                    if not LocusInstance and (i.type == "mRNA" or i.type == "CDS") :
                        if record.molType != 'n' :
                            LocusInstance = Locus(GeneInstance.newLocusTag())
                            GeneInstance.transcriptList.append(LocusInstance)
                        else :
                            LocusInstance = GeneInstance.transcriptList[0]

                    if not LocusInstance and i.type == "exon" :
                        if GeneInstance.transcriptList :
                            LocusInstance = GeneInstance.transcriptList[0]
                        else :
                            LocusInstance = Locus(GeneInstance.newLocusTag())
                            GeneInstance.transcriptList.append(LocusInstance)
                    # /RESOLV

                    if i.type == "mRNA" :
                        PListInstance = PList()
                        LocusInstance.mRNA = PListInstance
    
                        posList = self.__locationList2posList(i)
                        if posList != None :
                            PListInstance.location = \
                                self.__location2pos(i.location)
                            PListInstance.positionList = posList
                        #if
                        if i.qualifiers.has_key("transcript_id") :
                            LocusInstance.transcriptID = \
                                i.qualifiers["transcript_id"][0]
                            LocusInstance.link = self.__transcriptToProtein(
                                LocusInstance.transcriptID.split('.')[0])

                        LocusInstance.locusTag = locusTag
                    #if
                    if i.type == "CDS" :
                        PListInstance = PList()
                        LocusInstance.CDS = PListInstance
    
                        PListInstance.location = self.__location2pos(i.location)
                        PListInstance.positionList = \
                            self.__locationList2posList(i)
    
                        if i.qualifiers.has_key("transl_table") :
                            LocusInstance.txTable = \
                                int(i.qualifiers["transl_table"][0])
                        if i.qualifiers.has_key("protein_id") :
                            LocusInstance.proteinID = \
                                i.qualifiers["protein_id"][0]
                            LocusInstance.link = \
                                LocusInstance.proteinID.split('.')[0]

                        LocusInstance.locusTag = locusTag
                    #if
                    if i.type == "exon" :
                        if not LocusInstance.exon :
                            LocusInstance.exon = PList()
                        LocusInstance.exon.positionList.extend(
                            self.__location2pos(i.location))
                    #if
                #if                            
            #if
        #for
        return record
    #parseRecord
#GBparser
