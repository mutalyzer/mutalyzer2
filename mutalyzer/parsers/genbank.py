"""
Module contains one public function create_record which returns a
mutalyzer GenRecord. Record populated with data from a GenBank file.
"""


from __future__ import unicode_literals

import codecs
import re
import bz2
from itertools import izip_longest

from Bio import SeqIO
from Bio.Alphabet import ProteinAlphabet

from .. import ncbi
from ..GenRecord import PList, Locus, Gene, Record


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
        @type name: unicode
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
    def __location2pos(self, location, require_exact=True):
        """
        Convert a location object to a tuple of integers.

        @arg location: A location object (see the BioPython documentation)
        @type location: location object
        @arg require_exact: Require exact positions.
        @type require_exact: bool

        @return: A tuple of integers
        @rtype: list
        """

        ret = []

        if require_exact:
            if not unicode(location.start).isdigit() or \
               not unicode(location.end).isdigit() :
                return None

        ret.append(location.start.position + 1)
        ret.append(location.end.position)

        return ret
    #__location2pos

    def __location2posList(self, location, require_exact=True):
        """
        Convert a location object to a list of integers.

        @arg location: A location object (see the BioPython documentation)
        @type location: location object
        @arg require_exact: Require exact positions.
        @type require_exact: bool

        @return: A list (of even length) of integers
        @rtype: list (integers)
        """

        ret = []

        if require_exact:
            if not unicode(location.start).isdigit() or \
               not unicode(location.end).isdigit() :
                return None
        #if

        for part in location.parts[::location.strand]:
            pos = self.__location2pos(part, require_exact=require_exact)
            if not pos:
                return None

            ret.append(pos[0])
            ret.append(pos[1])
            #if
        #for

        if not ret:
            # No subfeatures found, in that case just use the feature itself
            # as if it were its only subfeature.
            ret = self.__location2pos(location, require_exact=require_exact)

        return ret
    #__location2posList

    def _find_mismatch(self, sentences):
        """
        Find the indices of the first and last words that distinguishes one
        sentence from another. The index of the last word is counted backwards.

        @arg sentences: A list of sentences.
        @type sentences: list of unicode strings

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
        lists = [s.split() for s in sentences]

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
        @type key: unicode
        """

        if locus.qualifiers.has_key(key) :
            setattr(locus, key, locus.qualifiers[key][0])
        else :
            setattr(locus, key, "")
    #__tagByDict

    def __tagLocus(self, locusList):
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
                accession, version = i.transcript_id.split('.')
                try:
                    # We ignore the version.
                    i.proteinLink = ncbi.transcript_to_protein(
                        accession, int(version), match_version=False)[0]
                except ncbi.NoLinkError:
                    pass
            i.original_location = i.location
            if i.ref:
                # This is a workaround for a bug in BioPython.
                # But seriously I have no idea for which bug and couldn't find
                # any hints in the commit history. So I just copied it over
                # with the last changes to this code, but it can probably be
                # removed.
                i.positionList = None
            else:
                i.positionList = self.__location2posList(i.location)
            i.location = self.__location2pos(i.location) #FIXME
            #if not i.positionList : # FIXME ???
            #    i.positionList = i.location
            if i.positionList :
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
        @type tagName: unicode
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
        if not mrnaList :
            # If the mRNA doesn't have exact positions (e.g., it's annotated
            # at `join(<1..11,214..548,851..4143)`), we still want to use the
            # part that is in this reference for matching.
            mrnaList = self.__location2posList(mrna.original_location,
                                               require_exact=False)

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

    def link(self, rnaList, cdsList):
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
        self.__tagLocus(rnaList)
        self.__tagLocus(cdsList)

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

    def create_record(self, filename):
        """
        Create a GenRecord.Record from a GenBank file

        @arg filename: The full path to the compressed GenBank file
        @type filename: unicode

        @return: A GenRecord.Record instance
        @rtype: object (record)
        """
        # first create an intermediate genbank record with BioPython
        file_handle = bz2.BZ2File(filename, "r")
        file_handle = codecs.getreader('utf-8')(file_handle)
        biorecord = SeqIO.read(file_handle, "genbank")
        file_handle.close()

        record = Record()
        record.seq = biorecord.seq

        # Note: The .source_* values may be different from the values we are
        #     working with, e.g. for UD slices where these values (taken from
        #     the genbank file) are from the original NC reference. We try to
        #     set the .id field to the working value in the caller.
        record.source_id = biorecord.id
        try:
            record.source_accession, record.source_version = biorecord.id.split('.')[:2]
        except ValueError:
            record.source_accession = biorecord.id
            record.source_version = '1'
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
                self.link(myGene.rnaList, myGene.cdsList)
                for i in myGene.rnaList :
                    myRealGene = record.findGene(i.gene)
                    version = myRealGene.newLocusTag()
                    # TODO: Here we discard transcripts that are not complete
                    # in this reference, but it might be nicer to still keep
                    # them so that we can (for example) show them in the
                    # legend. Of course they should still not be allowed to be
                    # selected in the variant description.
                    # (Same for leftover CDS features below.)
                    if i.usable :
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
                            # (Same for leftover CDS features below.)
                            try:
                                version = LOCUS_TAG_VERSION.findall(
                                    i.locus_tag)[0].zfill(3)
                            except IndexError:
                                pass
                        myTranscript = Locus(version)
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

                # We now look for leftover CDS entries that were not linked to
                # any transcript. We add them and the RNA will be constructed
                # for them later.
                # This does mean that these transcripts always come last (and
                # are shown last in for example the legend).
                for i in myGene.cdsList :
                    if not i.linked:
                        myRealGene = record.findGene(i.gene)
                        version = myRealGene.newLocusTag()
                        if i.usable:
                            if i.locus_tag :
                                try:
                                    version = LOCUS_TAG_VERSION.findall(
                                        i.locus_tag)[0].zfill(3)
                                except IndexError:
                                    pass
                            myTranscript = Locus(version)
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
                    #if
                    myRealGene.transcriptList.append(myTranscript)
                #if
            #if
        #else

        # Discard genes for which we haven't constructed any transcripts.
        record.geneList = [gene for gene in record.geneList
                           if gene.transcriptList]

        return record
    #create_record
#GBparser
