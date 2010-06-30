#!/usr/bin/python

import Crossmap
import Bio

"""
    Module to convert a GenBank record to a nested dictionary consisting of
    a list of genes, which itself consists of a list of loci. This structure
    makes it possible to iterate over genes and transcripts without having to
    search for them each time.

    Public classes:
        PList     ; Store a general location and a list of splice sites.
        Locus     ; Store data about the mRNA and CDS splice sites.
        Gene      ; Store a list of Locus objects and the orientation.
        Record    ; Store a geneList and other additional information.
        GenRecord ; Convert a GenBank record to a nested dictionary.
"""

class PList(object) :
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
        self.positionList = []
    #__init__
#PList

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

    def __init__(self, name) :
        """
            Initialise the class.

            Public variables (altered):
                mRNA     ; A position list object.
                CDS      ; A position list object.
                location ;
                exon     ; A position list object.
                txTable  ; The translation table.
                CM       ; A Crossmap object.
        """

        self.name = name
        self.mRNA = None
        self.CDS = None
        self.location = []
        self.exon = None
        self.txTable = 1
        self.CM = None
        self.transcriptID = None
        self.proteinID = None
        self.molType = 'c'
        self.description = ""
        self.proteinDescription = "?"
        self.locusTag = None
        self.transLongName = ""
        self.protLongName = ""
    #__init__

    def addToDescription(self, rawVariant) :
        """
        """

        if self.description :
            self.description = "%s;%s" % (self.description, rawVariant)
        else :
            self.description = rawVariant
    #addToDescription
#Locus

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

    def __init__(self, name) :
        """
            Initialise the class.

            Public variables (altered):
                orientation ; The orientation of the gene.
                list        ; A list of Locus objects.
        """

        self.name = name
        self.orientation = 1
        self.transcriptList = []
        self.location = []
        self.longName = ""
    #__init__

    def findLocus(self, name) :
        """
        """

        for i in self.transcriptList :
            if i.name == name :
                return i
        return None
    #findLocus
#Gene

class Record(object) :
    """
        A Record object, to store a geneList and other additional 
        information.

        Special methods:
            __init__() ; Initialise the class.

        Public variables:
            geneList  ; List of Gene objects.
            mol_type  ; Variable to indicate the sequence type (DNA, RNA, ...)
            organelle ; Variable to indicate whether the sequence is from the
                        nucleus or from an onganelle (if so, also from which
                        one).
            source    ; A fake gene that can be used when no gene information
                        is present.
    """

    def __init__(self) :
        """
            Initialise the class.


            Public variables (altered):
                geneList  ; List of Gene objects.
                molType   ; Variable to indicate the sequence type (DNA, RNA,
                            ...)
                seq       ; The reference sequence
                mapping   ; The mapping of the reference sequence to the genome
                            include a list of differences between the sequences
                organelle ; Variable to indicate whether the sequence is from
                            the nucleus or from an onganelle (if so, also from
                            which one).
                source    ; A fake gene that can be used when no gene
                            information is present.
        """

        self.geneList = []
        self.molType = ''
        self.seq = ""
        self.mapping = []
        self.organelle = None
        self.source = Gene(None)
        self.description = ""
        self._sourcetype = None           #LRG or GB
    #__init__

    def findGene(self, name) :
        """
        """

        for i in self.geneList :
            if i.name == name :
                return i
        return None
    #findGene

    def addToDescription(self, rawVariant) :
        """
        """

        if self.description :
            self.description = "%s;%s" % (self.description, rawVariant)
        else :
            self.description = rawVariant
    #addToDescription
#Record

class GenRecord() :
    """
        Convert a GenBank record to a nested dictionary.

        Public methods:
            checkRecord()   ;   Check and repair self.record
    """

    def __init__(self, config, output) :
        """
        """

        self.__config = config
        self.__output = output
        self.record = None
    #__init__

    def __constructCDS(self, mRNA, CDSpos) :
        """
        """
    
        i = 1
        ret = [CDSpos[0]]
    
        while CDSpos[0] > mRNA[i] :
            i += 2
    
        j = i
        while CDSpos[1] > mRNA[j] :
            j += 2
    
        ret.extend(mRNA[i:j])
        ret.append(CDSpos[1])
    
        return ret
    #__constructCDS

    def __maybeInvert(self, gene, string) :
        """
        """
    
        if gene.orientation == -1 :
            return Bio.Seq.reverse_complement(string)
        return string
    #__maybeInvert

    def checkRecord(self) :
        """
            Check if the record in self.record is compatible with mutalyzer

            update the mRNA PList with the exon and CDS data
        """
        #TODO:  This function should really check 
        #       the record for minimal requirements.
        for i in self.record.geneList :
            for j in i.transcriptList :
                if not j.mRNA :
                    if not j.exon:
                        self.__output.addMessage(__file__, 2, "WNOMRNA",
                            "No mRNA field found for gene %s, transcript " \
                            "variant %s in record, constructing " \
                            "it from CDS." % (i.name, j.name))
                        if j.CDS :
                            if not j.CDS.positionList :
                                self.__output.addMessage(__file__, 2,
                                    "WNOCDSLIST", "No CDS list found for " \
                                    "gene %s, transcript variant %s in " \
                                    "record, constructing it from " \
                                    "CDS location." % (i.name, j.name))
                                j.mRNA = j.CDS
                                j.mRNA.positionList = j.CDS.location
                            #if
                            else :
                                j.mRNA = j.CDS
                        #if
                        else :
                            self.__output.addMessage(__file__, 2, "WNOCDS",
                                "No CDS found for gene %s, transcript " \
                                "variant %s in record, " \
                                "constructing it from genelocation." % (
                                i.name, j.name))
                            j.mRNA = PList()
                            j.mRNA.positionList = i.location
                            j.molType = 'n'
                        #else
                    #if
                    else :
                        self.__output.addMessage(__file__, 2, "WNOMRNA",
                            "No mRNA field found for gene %s, transcript " \
                            "variant %s in record, constructing " \
                            "it from gathered exon information." % (
                            i.name, j.name))
                        j.mRNA = j.exon
                    #else
                #if
                if not j.mRNA.positionList :
                    j.mRNA.positionList = j.mRNA.location
                if j.CDS :
                    if not j.CDS.positionList :
                        self.__output.addMessage(__file__, 2, "WNOCDS",
                            "No CDS list found for gene %s, transcript " \
                            "variant %s in record, constructing " \
                            "it from mRNA list and CDS location." % (i.name, 
                            j.name))
                        if j.mRNA.positionList :
                            j.CDS.positionList = self.__constructCDS(
                                j.mRNA.positionList, j.CDS.location)
                        else :
                            j.CDS.positionList = self.__constructCDS(
                                j.mRNA.location, j.CDS.location)
                    #if
                    j.CM = Crossmap.Crossmap(j.mRNA.positionList, 
                                             j.CDS.location, i.orientation)
                #if
                else :
                    j.molType = 'n'
                    if j.mRNA.positionList :
                        j.CM = Crossmap.Crossmap(j.mRNA.positionList, 
                                                 [], i.orientation)
                    else :
                        j.description = '?'
                #else                                                 
            #for
        #for
    #checkRecord

    def name(self, start_g, stop_g, varType, arg1, arg2, roll) :
        """
        """

        forwardStart = start_g
        forwardStop = stop_g
        reverseStart = stop_g
        reverseStop = start_g
        if roll :
            forwardStart += roll[1]
            forwardStop += roll[1]
            reverseStart -= roll[0]
            reverseStop -= roll[0]
        #if

        if varType != "subst" :
            if forwardStart != forwardStop :
                self.record.addToDescription("%s_%s%s%s" % (forwardStart, 
                    forwardStop, varType, arg1))
            else :
                self.record.addToDescription("%s%s%s" % (forwardStart, varType,
                                                         arg1))
        #if
        else :
            self.record.addToDescription("%s%c>%c" % (forwardStart, arg1, arg2))

        for i in self.record.geneList :
            for j in i.transcriptList :
                if j.CM :
                    orientedStart = forwardStart
                    orientedStop = forwardStop
                    if i.orientation == -1 :
                        orientedStart = reverseStart
                        orientedStop = reverseStop
                    #if

                    # Check whether the variant hits CDS start.
                    if j.molType == 'c' and forwardStop >= j.CM.x2g(1, 0) \
                       and forwardStart <= j.CM.x2g(3, 0) :
                        self.__output.addMessage(__file__, 2, "WSTART", 
                            "Mutation in start codon of gene %s transcript " \
                            "%s." % (i.name, j.name))

                    # FIXME Check whether the variant hits a splice site.

                    if varType != "subst" :
                        if orientedStart != orientedStop :
                            j.addToDescription("%s_%s%s%s" % (
                                j.CM.g2c(orientedStart), j.CM.g2c(orientedStop),
                                varType, self.__maybeInvert(i, arg1)))
                        else :
                            j.addToDescription("%s%s%s" % (
                                j.CM.g2c(orientedStart), varType, 
                                self.__maybeInvert(i, arg1)))
                    #if
                    else :
                        j.addToDescription("%s%c>%c" % (j.CM.g2c(orientedStart),
                            self.__maybeInvert(i, arg1), 
                            self.__maybeInvert(i, arg2)))
                #if
            #for
        #for
    #name
#GenRecord

if __name__ == "__main__" :
    R = GenRecord()
    del R
#if
