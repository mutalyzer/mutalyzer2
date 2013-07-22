"""
Module contains one public function create_record which returns a
mutalyzer GenRecord. Record populated with data from a GenBank file.
"""


import re
import bz2
from itertools import izip_longest

from Bio import SeqIO, Entrez
from Bio.Alphabet import ProteinAlphabet

from mutalyzer import config
from mutalyzer import Db
from mutalyzer.GenRecord import PList, Locus, Gene, Record, GenRecord


# Regular expression used to find version number in locus tag
LOCUS_TAG_VERSION = re.compile('\d{1,3}$')


class tempGene():
    """
    Container class for a given gene name.

    Special methods:
        - __init__(name) ; Initialise the class.

    Public variables:
        - rnaList ; List of splice sites.
        - cdsList ; CDS list (including internal splice sites).
    """

    def __init__(self, name):
        """
        Initialise the class for a given gene name.

        Public variables:
            - rnaList ; List of splice sites.
            - cdsList ; CDS list (including internal splice sites).

        @arg name: Gene name
        @type name: string
        """

        self.name = name
        self.rnaList = []
        self.cdsList = []
    #__init__
#tempGene


class GBparser():
    """
    @todo: documentation
    """
    def __init__(self):
        """
        Initialise the class

        Private variables:
            - __database ; Db.Cache object
        """
        Entrez.email = config.get('email')
        self.__database = Db.Cache()
    #__init__

    def __location2pos(self, location):
        """
        Convert a location object to a tuple of integers.

        @arg location: A location object (see the BioPython documentation)
        @type location: location object

        @return: A tuple of integers
        @rtype: list
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

    def __locationList2posList(self, locationList):
        """
        Convert a list of locations to a list of integers.

        @arg locationList: A list of locations (see the BioPython documentation)
        @type locationList: list (location objects)

        @return: A list (of even length) of integers
        @rtype: list (integers)
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
            if temp :
                ret.append(temp[0])
                ret.append(temp[1])
            #if
        #for

        return ret
    #__locationList2posList

    def __transcriptToProtein(self, transcriptAcc):
        """
        Try to find the protein linked to a transcript id.

        First look in our database, if a link can not be found, try to
        retrieve it via the NCBI. Store the result in our database.

        @arg transcriptAcc: Accession number of the transcript for which we
                            want to find the protein
        @type transcriptAcc: string

        @return: Accession number of a protein or None if nothing can be found
        @rtype: string
        """
        try:
            return self.__database.getProtAcc(transcriptAcc)
        except IndexError:
            pass

        handle = Entrez.esearch(db = "nucleotide", term = transcriptAcc)
        result = Entrez.read(handle)
        handle.close()

        transcriptGI = result["IdList"][0]

        handle = Entrez.elink(dbfrom = "nucleotide", db = "protein",
                              id = transcriptGI)
        result = Entrez.read(handle)
        handle.close()

        if not result[0]["LinkSetDb"] :
            self.__database.insertLink(transcriptAcc, None)
            return None

        proteinGI = result[0]["LinkSetDb"][0]["Link"][0]["Id"]

        handle = Entrez.efetch(db='protein', id=proteinGI, rettype='acc', retmode='text')

        proteinAcc = handle.read().split('.')[0]
        handle.close()

        self.__database.insertLink(transcriptAcc, proteinAcc)
        return proteinAcc
    #__transcriptToProtein

    def _find_mismatch(self, sentences):
        """
        Find the indices of the first and last words that distinguishes one
        sentence from another. The index of the last word is counted backwards.

        @arg sentences: A list of sentences.
        @type sentences: list of strings

        @return: The indices of the words where sentences start to differ,
            both are -1 when no mismatches are found.
        @rtype: tuple(int, int)

        Example usage:

            >>> parser._find_mismatch(['a b c d e' , 'a B c d e', 'a b C d e'])
            (1, 2)
            >>> parser._find_mismatch(['a b c d e' , 'a b c d e', 'a b C D e'])
            (2, 1)
            >>> parser._find_mismatch(['a b c' , 'a b c', 'a b c'])
            (-1, -1)

        Note: The result can be used to slice the mismatching part from the
            sentences where you take the negative value of the second returned
            index. For the second example above:

                >>> 'a b c d e'.split()[1:-2]
                ['b', 'c']

            But be careful since the second index may be 0, but slicing syntax
            does not permit taking the -0 index from the right:

                >>> 'a b c d e'.split()[2:0] == ['c', 'd', 'e']
                False

            Although less elegant, first check the second index for 0 and in
            that case leave it out:

                >>> 'a b c d e'.split()[2:] == ['c', 'd', 'e']
                True

            The case where no mismatch is found just works, since slicing with
            [-1:1] yields the empty list.
        """
        # Create lists of words
        lists = map(str.split, sentences)

        try:
            forward, reverse = [next(i for i, v in
                                     enumerate(izip_longest(*lists))
                                     if not len(set(v)) <= 1)
                                for lists in (lists, map(reversed, lists))]
        except StopIteration:
            # No mismatch found
            forward = reverse = -1

        return forward, reverse
    #_find_mismatch

    def __tagByDict(self, locus, key):
        """
        Transfer a variable in the qualifiers dictionary to the locus
        object. If the variable does not exist, set it to the empty string.

        @arg locus: The locus object on which the transfer should be performed
        @type locus: locus object
        @arg key: The name of the variable that should be transferred
        @type key: string
        """

        if locus.qualifiers.has_key(key) :
            setattr(locus, key, locus.qualifiers[key][0])
        else :
            setattr(locus, key, "")
    #__tagByDict

    def __tagLocus(self, locusList, output):
        """
        Enrich a list of locus objects (mRNA or CDS) with information used
        for linking (locus_tag, proteinLink and productTag). Also
        transfer the variables transcript_id, protein_id, gene and product
        to each of the locus objects. If these variables do not exist, set
        them to the empty string.

        @arg locusList: A list of locus objects
        @type locusList: list
        """

        productList = []
        for i in locusList :
            # Transfer some variables from the dictionary to the locus object.
            self.__tagByDict(i, "locus_tag")
            self.__tagByDict(i, "transcript_id")
            self.__tagByDict(i, "protein_id")
            self.__tagByDict(i, "gene")
            self.__tagByDict(i, "product")

            # Gather the product tags.
            productList.append(i.product)

            i.proteinLink = None
            i.linked = False
            if not i.transcript_id :
                if i.protein_id : # Tag a CDS with the protein id.
                    i.proteinLink = i.protein_id.split('.')[0]
            #if
            else :                # Tag an mRNA with the protein id too.
                i.proteinLink = \
                    self.__transcriptToProtein(i.transcript_id.split('.')[0])
            i.positionList = self.__locationList2posList(i)
            i.location = self.__location2pos(i.location) #FIXME
            #if not i.positionList : # FIXME ???
            #    i.positionList = i.location
            if i.positionList or i.location :
                i.usable = True
            else :
                i.usable = False
                
        #for

        if productList :
            # Find the defining words in the product list.
            a, b = self._find_mismatch(productList)

            # Add the defining words to the locus.
            for i in range(len(locusList)):
                if b == 0:
                    locusList[i].productTag = \
                        ' '.join(productList[i].split()[a:])
                else:
                    locusList[i].productTag = \
                        ' '.join(productList[i].split()[a:-b])
        #if
    #__tagLocus


    def __checkTags(self, locusList, tagName):
        """
        Check whether all tags in a locus list are unique. Prune all the
        non unique tags.

        @arg locusList: A list of loci
        @type locusList: list
        @arg tagName: Name of the tag to be checked
        @type tagName: string
        """

        tags = []
        for i in locusList : # Gather all the tags.
            tags.append(getattr(i, tagName))

        badTags = []
        for i in locusList : # Collect the tags that can not be used.
            myTag = getattr(i, tagName)
            numberOfTags = tags.count(myTag)
            if numberOfTags > 1 :
                badTags.append(myTag)
        #for

        for i in locusList : # Remove unusable tags.
            if getattr(i, tagName) in badTags :
                setattr(i, tagName, None)
        #for
    #__checkTags

    def __matchByRange(self, mrna, cds):
        """
        Match the mRNA list to the CDS list.

        @arg mrna: List of splice sites
        @type mrna: list
        @arg cds: CDS list (including internal splice sites)
        @type cds: list

        @return:
            - E{-}1 : False
            - 0 : Don't know
            - 1 : Maybe true
            - 2 : Probably true
        @rtype: integer
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

    def link(self, rnaList, cdsList, output):
        """
        Link mRNA loci to CDS loci (all belonging to one gene).

        First of all, the range of the CDS must be a subrange of that of
        the mRNA. If this is true, then we try to link both loci. The first
        method is by looking at the locus_tag, if this fails, we try to
        match the proteinLink tags, if this also fails, we try the
        productTag.

        If no link could be found, but there is only one possibility left,
        the loci are linked too.

        The method that was used to link the loci, is put in the linkmethod
        variable of the transcript locus. The link variable of the
        transcript locus is a pointer to the CDS locus. Furthermore, the
        linked variable of the CDS locus is set to indicate that this locus
        is no longer available for linking.

        Available link methods are: locus, protein, product and exhaustion.

        @arg rnaList: A list of mRNA loci
        @type rnaList: list
        @arg cdsList: A list of CDS loci
        @type cdsList: list
        """

        # Enrich the lists with as much information we can find.
        self.__tagLocus(rnaList, output)
        
        self.__tagLocus(cdsList, output)

        # Prune the tags based upon uniqueness.
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
                    # Try to link via the locus tag first.
                    if i.locus_tag and i.locus_tag == j.locus_tag :
                        i.link = j
                        i.linkMethod = "locus"
                        j.linked = True
                        #print "Linked:", j.locus_tag
                        break
                    #if
                    # Try the proteinLink tag.
                    if i.proteinLink and i.proteinLink == j.proteinLink :
                        i.link = j
                        i.linkMethod = "protein"
                        j.linked = True
                        break
                    #if
                    # Try the productTag.
                    if i.productTag and i.productTag == j.productTag :
                        i.link = j
                        i.linkMethod = "product"
                        j.linked = True
                        break
                    #if
                #if
            #for

        # Now look if there is only one possibility left.
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

    def create_record(self, filename, output):
        """
        Create a GenRecord.Record from a GenBank file

        @arg filename: The full path to the compressed GenBank file
        @type filename: string

        @return: A GenRecord.Record instance
        @rtype: object (record)
        """
        # first create an intermediate genbank record with BioPython
        file_handle = bz2.BZ2File(filename, "r")
        biorecord = SeqIO.read(file_handle, "genbank")
        file_handle.close()

        record = Record()
        record.seq = biorecord.seq

        # Note: The .source_* values may be different from the values we are
        #     working with, e.g. for UD slices where these values (taken from
        #     the genbank file) are from the original NC reference. We try to
        #     set the .id field to the working value in the caller.
        record.source_id = biorecord.id
        record.source_accession, record.source_version = biorecord.id.split('.')[:2]
        record.source_gi = biorecord.annotations['gi']
        record.organism = biorecord.annotations['organism']

        # Todo: This will change once we support protein references
        if isinstance(biorecord.seq.alphabet, ProteinAlphabet):
            return record

        exonList = []
        geneDict = {}

        accInfo = biorecord.annotations['accessions']
        if len(accInfo) >= 3 and accInfo[1] == "REGION:":
            # Todo: This information is present in the genbank file if it is a
            #     UD sliced from a chromosome. We can get the same information
            #     for NM references from our mapping database and that way
            #     also provide chromosomal variant descriptions for those.
            region = accInfo[2]
            if "complement" in region :
                record.orientation = -1
                record.chromOffset = int(region.split('.')[2][:-1])
            #if
            else :
                record.chromOffset = int(accInfo[2].split('.')[0])
        #if
        for i in biorecord.features :
            if i.qualifiers :
                if i.type == "source" :
                    if i.qualifiers.has_key("mol_type") :
                        if i.qualifiers["mol_type"][0] in ["mRNA", \
                           "transcribed RNA"] :
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
                    geneName = i.qualifiers["gene"][0]
                    if i.type == "gene" :
                        if not geneDict.has_key(geneName) :
                            myGene = Gene(geneName)
                            record.geneList.append(myGene)
                            if i.strand :
                                myGene.orientation = i.strand
                            myGene.location = self.__location2pos(i.location)
                            geneDict[geneName] = tempGene(geneName)
                        #if
                    #if

                    if i.type in ["mRNA", "misc_RNA", "ncRNA", "rRNA", "tRNA",
                       "tmRNA"] :
                        geneDict[geneName].rnaList.append(i)
                    if i.type == "CDS" :
                        geneDict[geneName].cdsList.append(i)
                    if i.type == "exon" :
                        exonLocation = self.__location2pos(i.location)
                        if exonLocation :
                            exonList.extend(exonLocation)
                    #if
                #if
            #if
        #for
        if record.molType in ['g', 'm'] :
            for j in geneDict.keys() :
                myGene = geneDict[j]
                self.link(myGene.rnaList, myGene.cdsList, output)
                for i in myGene.rnaList :
                    if i.usable :

                        myRealGene = record.findGene(i.gene)
                        if i.locus_tag :
			    
                            # Note: We use the last three characters of the
                            # locus_tag as a unique transcript version id.
                            # This is also used to for the protein-transcript
                            # link table.
                            # Normally, locus_tag ends with three digits, but
                            # for some (e.g. mobA on NC_011228, a plasmid) it
                            # ends with two digits prepended with an
                            # underscore. Or prepended with a letter. We
                            # really want a number, so 'fix' this by only
                            # looking for a numeric part.
                            try:
                                version = LOCUS_TAG_VERSION.findall(
                                    i.locus_tag)[0].zfill(3)
                            except IndexError:
                                version = '000'

                            myTranscript = Locus(version)
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
                            if "transl_except" in i.link.qualifiers:
                                myTranscript.transl_except=self.create_exception(i.link)
                               
                        #if
                        myRealGene.transcriptList.append(myTranscript)
                    #if
                    else:
                        output.addMessage(__file__, 2, 'WPOSITION',
                              "The gene's %s coordinates extend beyound transcript" % i.gene)
                #for
                for i in myGene.cdsList :
                    if not i.linked and \
                       (i.usable or not geneDict[myGene.name].rnaList) :
              
                        myRealGene = record.findGene(i.gene)
                        if i.locus_tag :
                            # Note: We use the last three characters of the
                            # locus_tag as a unique transcript version id.
                            # This is also used to for the protein-transcript
                            # link table.
                            # Normally, locus_tag ends with three digits, but
                            # for some (e.g. mobA on NC_011228, a plasmid) it
                            # ends with two digits prepended with an
                            # underscore. Or prepended with a letter. We
                            # really want a number, so 'fix' this by only
                            # looking for a numeric part.
                            try:
                                version = LOCUS_TAG_VERSION.findall(
                                    i.locus_tag)[0].zfill(3)
                            except IndexError:
                                version = '000'
                            myTranscript = Locus(version)
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
                        if "transl_except" in i.qualifiers:
                            myTranscript.transl_except=self.create_exception(i)
                           
                        myRealGene.transcriptList.append(myTranscript)
                        #if
                    #if
                #for
            #for
        #if
        else :
            if geneDict :
                myGene = geneDict[geneDict.keys()[0]]
                myRealGene = record.geneList[0]
                if myGene.cdsList :
                    myCDS = myGene.cdsList[0]
                    self.__tagByDict(myCDS, "protein_id")
                    self.__tagByDict(myCDS, "product")
                #if
                else :
                    myCDS = None
                myTranscript = Locus("001")
                myTranscript.exon = PList()
                if exonList :
                    myTranscript.exon.positionList = exonList
                else :
                    myTranscript.exon.location = myRealGene.location
                if myCDS :
                    myTranscript.CDS = PList()
                    myTranscript.CDS.location = \
                        self.__location2pos(myCDS.location)
                #if
                if exonList or myRealGene.location or \
                   myTranscript.CDS.location :
                    myTranscript.transcriptID = biorecord.id
                    if myCDS :
                        myTranscript.proteinID = myCDS.protein_id
                        myTranscript.proteinProduct = myCDS.product
                        myTranscript.linkMethod = "exhaustion"
                        myTranscript.transcribe = True
                        if myCDS.qualifiers.has_key("transl_table") :
                            myTranscript.txTable = \
                                int(i.qualifiers["transl_table"][0])
                        if "transl_except" in myCDS.qualifiers:
                            myTranscript.transl_except=self.create_exception(myCDS)
                           
                    #if
                    myRealGene.transcriptList.append(myTranscript)
                #if
            #if
        #else
        for i in record.geneList :
            if not i.transcriptList :
                record.geneList.remove(i)

        return record
    #create_record

    def create_exception(self, SeqFeature):
	''' create-exception return a list with tuples of coordinates of start and end positions of uncanonical amino acids.
	If your entry does not contain it, this function return an empty list.'''

	sec_coord_list=[]
	# We are looking for translational exception: single codon 
	# the translation of which does not conform to genetic code
	# indicated by transl_table. For more information see 
	# http://www.ddbj.nig.ac.jp/sub/ref6-e.html#transl_except'''
	for transl_except in SeqFeature.qualifiers["transl_except"]:
		intermediate=re.split("[,:.]", transl_except.strip("()"))
		triplet_dict={"Ala":"A", "Gly":"G", "Val":"V", "Leu":"L", "Ile":"I",
			      "Met":"M", "Phe":"F", "Asn":"N", "Gln":"Q", "Asp":"D",
                  "Glu":"E", "His":"H", "Lys":"K", "Arg":"R", "Ser":"S",
			      "Thr":"T", "Tyr":"Y", "Trp":"W", "Cys":"C", "Pro":"P", 
			      "Sec":"U", "Pyl":"O", "TERM":"Stop", "OTHER": "X", "Asx" : "B", 
                  "Glx" : "Z", "Xle" : "J"}
		sec_coord_list.append((int(intermediate[1]), triplet_dict[intermediate[-1]], "g."))			
	#for
        return sec_coord_list
    #create_exception
   

#GBparser
