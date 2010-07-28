from Bio import SeqIO, Entrez  # read()
import bz2                     # BZ2Compressor(), BZ2File()
from GenRecord import PList, Locus, Gene, Record, GenRecord
import Db

"""
    Module contains one public function createGBRecord which returns a
    mutalyzer GenRecord.Record populated with data from a GenBank file.
"""

class tempGene() :
    """
    """

    def __init__(self, name) :
        """
        """

        self.name = name
        self.rnaList = []
        self.cdsList = []
    #__init__
#tempGene    

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

            if not result[0]["LinkSetDb"] :
                return None

            proteinGI = result[0]["LinkSetDb"][0]["Link"][0]["Id"]

            handle = Entrez.efetch(db = "protein", id = proteinGI,
                                   rettype = "acc")

            proteinAcc = handle.read().split('.')[0]
            handle.close()

            self.__database.insertLink(transcriptAcc, proteinAcc)
        #if

        return proteinAcc
    #__transcriptToProtein

    def __findMismatch(self, productList, direction) :
        """
        """

        i = 0
        while i < productList[0].count(' ') + 1 :
            for j in range(1, len(productList)) :
                if productList[0][::direction].split(' ')[i] != \
                   productList[j][::direction].split(' ')[i] :
                    return i
            i += 1
        #while
        return 0
    #__findMismatch

    def __tagByDict(self, locus, key) :
        """
        """

        if locus.qualifiers.has_key(key) :
            setattr(locus, key, locus.qualifiers[key][0])
        else :
            setattr(locus, key, None)
    #__tagByDict

    def __tagLocus(self, locusList) :
        """
        """

        productList = []
        for i in locusList :
            productList.extend(i.qualifiers["product"])
            self.__tagByDict(i, "locus_tag")
            self.__tagByDict(i, "transcript_id")
            self.__tagByDict(i, "protein_id")
            self.__tagByDict(i, "gene")
            self.__tagByDict(i, "product")
            i.proteinLink = None
            i.linked = False
            if not i.transcript_id :
                if i.protein_id :
                    i.proteinLink = i.protein_id.split('.')[0]
            #if                    
            else :                    
                i.proteinLink = \
                    self.__transcriptToProtein(i.transcript_id.split('.')[0])
            i.positionList = self.__locationList2posList(i)
            i.location = self.__location2pos(i.location) #FIXME
            #if not i.positionList : # FIXME ???
            #    i.positionList = i.location
            if i.positionList :
                i.usable = True
            else :
                i.usable = False
        #for

        if productList :
            a = self.__findMismatch(productList, 1)
            b = productList[0].count(' ') - \
                self.__findMismatch(productList, -1) + 1

            for i in range(len(locusList)) :
                locusList[i].productTag = \
                    ' '.join(productList[i].split(' ')[a:b])
        #if
    #__tagLocus


    def __checkTags(self, locusList, tagName) :
        """
            Check whether all tags in a locus list are unique.

            Arguments:
                locusList ; A list of loci.
                tagName   ; Name of the tag to be checked.

            Returns:
                boolean ; True if the tags are unique, False otherwise.
        """

        tags = []
        for i in locusList :
            tags.append(getattr(i, tagName))
        badTags = []
        for i in locusList :
            myTag = getattr(i, tagName)
            numberOfTags = tags.count(myTag)
            if numberOfTags > 1 :
                badTags.append(myTag)
        #for
        for i in locusList :
            if getattr(i, tagName) in badTags :
                setattr(i, tagName, None)
        #for
    #__checkTags

    def __matchByRange(self, mrna, cds) :
        """
            Match the mRNA list to the CDS list.

            Arguments:
                mrnaList ; List of splice sites.
                cdsList  ; CDS list (including internal splice sites).

            Returns:
                integer ; -1 : False.
                           0 : Don't know.
                           1 : Maybe true.
                           2 : Probably true.
        """

        if not cds or not mrna :
            return 0          # No information -> Don't know.

        mrnaList = mrna.positionList
        if not mrnaList :
            mrnaList = mrna.location
        cdsList = cds.positionList
        if not cdsList :
            cdsList = cds.location

        if not cdsList or not mrnaList :
            return 0          # No information -> Don't know.
        if cdsList[0] < mrnaList[0] or cdsList[-1] > mrnaList[-1] :
            return -1         # CDS is outside transcript range -> False.
        if len(cdsList) > 2 : # The CDS spans more than one exon.
            if not cdsList[1] in mrnaList :
                return -1     # At least one splice site doesn't match -> False.
            x = mrnaList.index(cdsList[1])
            y = x + len(cdsList) - 2
            if mrnaList[x:y] == cdsList[1:-1] :
                return 2 # All splice sites match -> Probably true.
            return -1    # At least one splice site doesn't match -> False.
        #if
        return 1         # Everything matches, but there is little information.
    #__matchByRange

    def link(self, rnaList, cdsList) :
        """
        """

        self.__tagLocus(rnaList)
        self.__tagLocus(cdsList)

        self.__checkTags(rnaList, "locus_tag")
        self.__checkTags(cdsList, "locus_tag")
        self.__checkTags(rnaList, "proteinLink")
        self.__checkTags(cdsList, "proteinLink")
        self.__checkTags(rnaList, "productTag")
        self.__checkTags(cdsList, "productTag")

        for i in rnaList :
            i.link = None
            i.linkMethod = None
            for j in cdsList :
                if self.__matchByRange(i, j) > 0 :
                    if i.locus_tag and i.locus_tag == j.locus_tag :
                        i.link = j
                        i.linkMethod = "locus"
                        j.linked = True
                        print "Linked:", j.locus_tag
                        break
                    #if
                    if i.proteinLink and i.proteinLink == j.proteinLink :
                        i.link = j
                        i.linkMethod = "protein"
                        j.linked = True
                        break
                    #if
                    if i.productTag and i.productTag == j.productTag :
                        i.link = j
                        i.linkMethod = "product"
                        j.linked = True
                        break
                    #if
                #if
            #for

        # One *could* also do exhaustion per matched range...
        for i in rnaList :
            if not i.link :
                leftOverCount = 0
                leftOverTranscript = None
                leftOverProtein = None
                for j in cdsList :
                    if self.__matchByRange(i, j) > 0 and not j.linked :
                        leftOverCount += 1
                        leftOverTranscript = i
                        leftOverProtein = j
                    #if
                #for
                if leftOverCount == 1 :
                    leftOverTranscript.link = leftOverProtein
                    leftOverTranscript.linkMethod = "exhaustion"
                    leftOverProtein.linked = True
                #if
            #if
        #for
    #link

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

        record.version = biorecord.id.split('.')[1]

        exonList = []
        geneDict = {}

        accInfo = biorecord.annotations['accessions']
        if len(accInfo) >= 3 :
            region = accInfo[2]
            if "complement" in region :
                record.orientation = -1
                record.chromOffset = int(region.split('.')[2][:-1])
            #if
            else :
                record.chromOffset = int(accInfo[2].split('.')[0])
        #if
        record.recordId = biorecord.id
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
                    if i.type == "gene" :
                        geneName = i.qualifiers["gene"][0]
                        if not geneDict.has_key(geneName) :
                            myGene = Gene(geneName)
                            record.geneList.append(myGene)
                            if i.strand :
                                myGene.orientation = i.strand
                            myGene.location = self.__location2pos(i.location)
                            geneDict[geneName] = tempGene(geneName)
                    #if

                    if i.type == "mRNA" :
                        geneDict[geneName].rnaList.append(i)
                    if i.type == "CDS" :
                        geneDict[geneName].cdsList.append(i)
                    if i.type == "exon" :
                        exonList.extend(self.__location2pos(i.location))
                #if
            #if
        #for
        if record.molType == 'g' :
            for j in geneDict.keys() :
                myGene = geneDict[j]
                self.link(myGene.rnaList, myGene.cdsList)
                for i in myGene.rnaList :
                    if i.usable :
                        myRealGene = record.findGene(i.gene)
                        if i.locus_tag :
                            myTranscript = Locus(i.locus_tag[-3:])
                        else :
                            myTranscript = Locus(myRealGene.newLocusTag())
                        myTranscript.mRNA = PList()
                        myTranscript.mRNA.positionList = i.positionList
                        myTranscript.mRNA.location = i.location
                        myTranscript.transcribe = True
                        myTranscript.transcriptID = i.transcript_id
                        myTranscript.transcriptProduct = i.product
                        myTranscript.locusTag = i.locus_tag
                        if i.link :
                            myTranscript.CDS = PList()
                            myTranscript.CDS.positionList = i.link.positionList
                            myTranscript.CDS.location = i.link.location
                            myTranscript.translate = True
                            myTranscript.proteinID = i.link.protein_id
                            myTranscript.linkMethod = i.linkMethod
                            myTranscript.proteinProduct = i.link.product
                            if i.link.qualifiers.has_key("transl_table") :
                                myTranscript.txTable = \
                                    int(i.qualifiers["transl_table"][0]) 
                        #if
                        myRealGene.transcriptList.append(myTranscript)
                    #if
                #for
                for i in myGene.cdsList :
                    if not i.linked and i.usable :
                        myRealGene = record.findGene(i.gene)
                        if i.locus_tag :
                            myTranscript = Locus(i.locus_tag[-3:])
                        else :
                            myTranscript = Locus(myRealGene.newLocusTag())
                        myTranscript.CDS = PList()
                        myTranscript.CDS.positionList = i.positionList
                        myTranscript.CDS.location = i.location
                        myTranscript.proteinID = i.protein_id
                        myTranscript.proteinProduct = i.product
                        if i.qualifiers.has_key("transl_table") :
                            myTranscript.txTable = \
                                int(i.qualifiers["transl_table"][0]) 
                        myRealGene.transcriptList.append(myTranscript)
                        #if
                    #if
                #for
            #for
        #if
        else :
            myGene = geneDict[geneDict.keys()[0]]
            myRealGene = record.geneList[0]
            myCDS = myGene.cdsList[0]
            self.__tagByDict(myCDS, "protein_id")
            self.__tagByDict(myCDS, "product")
            myTranscript = Locus("001")
            myTranscript.exon = PList()
            if exonList :
                myTranscript.exon.positionList = exonList
            else :
                myTranscript.exon.location = myRealGene.location
            myTranscript.CDS = PList()
            myTranscript.CDS.location = self.__location2pos(myCDS.location)
            myTranscript.transcriptID = biorecord.id
            myTranscript.proteinID = myCDS.protein_id
            myTranscript.proteinProduct = myCDS.product
            myTranscript.linkMethod = "exhaustion"
            myTranscript.transcribe = True
            if myCDS.qualifiers.has_key("transl_table") :
                myTranscript.txTable = \
                    int(i.qualifiers["transl_table"][0]) 
            myRealGene.transcriptList.append(myTranscript)
        #else
        for i in record.geneList :
            if not i.transcriptList :
                record.geneList.remove(i)

        return record
    #parseRecord
#GBparser
