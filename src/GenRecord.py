#!/usr/bin/python

class Plist(object) :
    """
        A position list object, to store a general location and a list of 
        specific splice sites (if available).

        These objects are used to describe either a list of mRNA splice sites
        or a list of CDS splice sites. These splice sites are stored in the
        list element. The location element is a fallback in case the splice
        sites are not available.

        Special methods:
            __init__() ; Initialise the class.

        Public variables:
            location ; A tuple of integers between which the object resides.
            list     ; A list (with an even amount of entries) of splice sites.
    """

    def __init__(self) :
        """
            Initialise the class.

            Public variables (altered):
                location ; A tuple of integers between which the object
                           resides. 
                list     ; A list (with an even amount of entries) of splice 
                           sites.
        """

        self.location = []
        self.list = []
    #__init__
#plist

class Locus(object) :
    """
        A Locus object, to store data about the mRNA and CDS splice sites.

        Special methods:
            __init__() ; Initialise the class.

        Public variables:
            mRNA ; A position list object.
            CDS  ; A position list object.
            exon ; A position list object.
    """

    def __init__(self) :
        """
            Initialise the class.

            Public variables (altered):
                mRNA ; A position list object.
                CDS  ; A position list object.
                exon ; A position list object.
                CM   ; A Crossmap object.
        """

        self.mRNA = None
        self.CDS = None
        self.exon = None
        self.CM = None

    #__init__
#locus

class Gene(object) :
    """
        A Gene object, to store a list of Locus objects and the orientation of 
        the gene.
        
        Special methods:
            __init__() ; Initialise the class.

        Public variables:
            orientation ; The orientation of the gene: 1 = forward, 
                                                      -1 = reverse.
            list        ; A list of Locus objects.
    """

    def __init__(self) :
        """
            Initialise the class.

            Public variables (altered):
                orientation ; The orientation of the gene.
                list        ; A list of Locus objects.
        """

        self.orientation = 0
        self.list = {}
    #__init__
#gene

class GenRecord() :
    """
        Hmmmm.

        Private methods:
            __location2pos(location)             ;
            __locationList2posList(locationList) ;

        Public methods:
    """

    def __location2pos(self, location) :
        """
            Convert a location object to a tuple of integers.
    
            Arguments:
                location ; A location object (see the BioPython documentation).
    
            Returns:
                List ; A tuple of integers.
        """

        ret = []
    
        ret.append(location.start.position + 1)
        ret.append(location.end.position)
    
        return ret
    #__location2pos
    
    def __locationList2posList(self, locationList) :
        """
            Convert a list of locations to a list of integers.
    
            Arguments:
                locationList ; A list of locations (see the BioPython 
                               documentation).
    
            Returns:
                List ; A list (of even length) of integers.
        """

        ret = []
    
        for i in locationList.sub_features :
            temp = self.__location2pos(i.location)
            ret.append(temp[0])
            ret.append(temp[1])
        #for
        
        return ret
    #__locationList2posList

    """
    def __sortins(self, position, posList) :
        last = 0

        for i in range(0, len(posList), 2) :
            if position[0] == posList[i] :
                return posList
            if position[0] > last and position[0] < posList[i] :
                return posList[:i] + position + posList[i:]
            last = posList[i]
        #for        
        return posList + position
    #__sortins
    """
    
    def record2dict(self, record) :
        recordDict = {}
        for i in  record.features :
            if i.qualifiers and i.qualifiers.has_key("gene") :
                gene = i.qualifiers["gene"][0]
                if not recordDict.has_key(gene) :
                    recordDict[gene] = Gene()
                if i.type == "gene" :
                    recordDict[gene].orientation = i.strand
    
                # Look if there is a locus tag present, if not, give it the
                # default tag `001'.
                locus_tag = "001"
                if i.qualifiers.has_key("locus_tag") :
                    locus_tag = i.qualifiers["locus_tag"][0][-3:]
                if not recordDict[gene].list.has_key(locus_tag) :
                    recordDict[gene].list[locus_tag] = Locus()
    
                if i.type == "mRNA" :
                    recordDict[gene].list[locus_tag].mRNA = Plist()
                    recordDict[gene].list[locus_tag].mRNA.location = \
                        self.__location2pos(i.location)
                    recordDict[gene].list[locus_tag].mRNA.list = \
                        self.__locationList2posList(i)
                #if
                if i.type == "CDS" :
                    recordDict[gene].list[locus_tag].CDS = Plist()
                    recordDict[gene].list[locus_tag].CDS.location = \
                        self.__location2pos(i.location)
                    recordDict[gene].list[locus_tag].CDS.list = \
                        self.__locationList2posList(i)
                #if
                if i.type == "exon" :
                    if not recordDict[gene].list[locus_tag].exon :
                        recordDict[gene].list[locus_tag].exon = Plist()
                    recordDict[gene].list[locus_tag].exon.list.extend(
                        self.__location2pos(i.location))
            #if
        #for

        return recordDict
    #record2dict
    
    def printRecordDict(self, d, record) :
        for i in d :
            print i
            print "  Orientation: " + str(d[i].orientation)
            for j in d[i].list :
                print "  Locus: " + str(j)
                if d[i].list[j].mRNA :
                    print "    mRNA: "
                    print "      " + str(d[i].list[j].mRNA.location)
                    if d[i].list[j].mRNA.list :
                        print "      " + str(d[i].list[j].mRNA.list)
                        print splice(record, d[i].list[j].mRNA.list)
                    #if
                    else :
                        print splice(record, d[i].list[j].mRNA.location)
                #if
                if d[i].list[j].CDS :
                    print "    CDS: "
                    print "      " + str(d[i].list[j].CDS.location)
                    if d[i].list[j].CDS.list :
                        print "      " + str(d[i].list[j].CDS.list)
                        print splice(record, d[i].list[j].CDS.list)
                    #if
                    else :
                        print splice(record, d[i].list[j].CDS.location)
                #if
            #for
        #for
    #printRecordDict
#GenRecord

if __name__ == "__main__" :
    R = GenRecord()
    bla = R._GenRecord__sortins([10, 20], [4, 5])
    print R._GenRecord__sortins([1, 2], bla)
    print R._GenRecord__sortins([8, 9], bla)
#if
