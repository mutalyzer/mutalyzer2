#!/usr/bin/python

"""
    Search for an NM number in the MySQL database, if the version number 
    matches, get the start and end positions in a variant and translate these
    positions to g. notation if the variant is in c. notation and vice versa.

    - If no end position is present, the start position is assumed to be the
      end position. 
    - If the version number is not found in the database, an error message is
      generated and a suggestion for an other version is given.
    - If the reference sequence is not found at all, an error is returned.
    - If no variant is present, the transcription start and end and CDS end
      in c. notation is returned.
    - If the variant is not accepted by the nomenclature parser, a parse error
      will be printed.

"""
import sys                     # argv
from Modules import Config     # Config()
from Modules import Db         # Db(), get_NM_version(), get_NM_info()
from Modules import Crossmap   # Crossmap(), g2x(), x2g(), main2int(), 
                               # offset2int(), info()
from Modules import Parser     # Nomenclatureparser(), parse()
from Modules import Output     # Output(), LogMsg()
from ZSI.fault import Fault    # Fault()
from ZSI import TC             # Struct()

from soaplib.serializers.primitive import String, Integer
from soaplib.serializers.clazz import ClassSerializer

class Mapping(ClassSerializer) :
    """
        Extended ClassSerializer object with mixed types of attributes
        
        Attributes:
            startmain ; Define the type of startmain.
            startoffset ; Define the type of startoffset.
            endmain ; Define the type of endmain value.
            endoffset ; Define the type of endoffset value.
            start_g ; Define the type of start_g value.
            end_g ; Define the type of end_g value.
            mutationType ; Define the type of mutation type
    """

    class types() :
        """
            Types are defined here for the soaplib module.
        """

        startmain = Integer
        startoffset = Integer
        endmain = Integer
        endoffset = Integer
        start_g = Integer
        end_g = Integer
        mutationType = String
    #types

    def __init__(self) :
        """
            Types are defined here for the TC module.
        """

        self.typecode = TC.Struct(Mapping, [ 
            TC.Integer('startmain'),
            TC.Integer('startoffset'),
            TC.Integer('endmain'),
            TC.Integer('endoffset'),
            TC.Integer('start_g'),
            TC.Integer('end_g'),
            TC.String('mutationType') 
            ], 'Mapping')
    #__init__                                                            
#Mapping

class Transcript(ClassSerializer) :
    """
        Extended ClassSerializer object with mixed types of attributes
        
        Attributes:
            trans_start ; Define the type of trans_start
            trans_stop  ; Define the type of trans_stop
            CDS_stop    ; Define the type of CDS_stop
    """

    class types() :
        """
        """

        trans_start = Integer
        trans_stop = Integer
        CDS_stop = Integer
    #types

    def __init__(self) :
        """
        """

        self.typecode = TC.Struct(Transcript, [ 
            TC.Integer('trans_start'),
            TC.Integer('trans_stop'),
            TC.Integer('CDS_stop') 
            ], 'Transcript')
    #__init__                                                            
#Transcript

def __sl2il(l) :
    """
        Convert a list of strings to a list of integers.

        Arguments: l ; A list of strings.

        Returns: list ; A list of integers.
    """

    for i in range(len(l)) :
        l[i] = int(l[i])
    return l
#__sl2il

def __getcoords(C, Loc, Type) :
    """
        Return main, offset and g positions given either a position in
        c. or in g. notation.

        Arguments:
            C    ; A crossmapper.
            Loc  ; Either a location in g. or c. notation.
            Type ; The reference type.
        Returns:
            triple:
                0 ; Main coordinate in c. notation.
                1 ; Offset coordinate in c. notation.
                2 ; Position in g. notation.
    """

    if Type == 'c' :
        main = C.main2int(Loc.MainSgn +  Loc.Main)
        offset = C.offset2int(Loc.OffSgn +  Loc.Offset)
        g = C.x2g(main, offset)
        main, offset = C.g2x(g)
    #if
    else :
        g = int(Loc.Main)
        main, offset = C.g2x(g)
    #else
    return (main, offset, g)
#__getcoords

def mrnaSplit(mrnaAcc) :
    """
        Split the full NM accession number into the number without the version
        and the version number.
        
        Arguments:
            mrnaAcc ; The NM accession number with version
            
        Returns:
            tuple  ; The NM accession number without version and version number.
    """

    accno = mrnaAcc
    version = 0
    if '.' in mrnaAcc :
        accno = mrnaAcc.split('.')[0]
        version = int(mrnaAcc.split('.')[1])
    return (accno, version)
#mrnaSplit


def makeCrossmap(build, acc, Conf) :
    '''
        Build the crossmapper

        Arguments:
            build   ; The human genome build
            acc     ; The NM accession number, including version.

        Returns:
            Cross ; A Crossmap object.
        
        
    '''

    # Make a connection to the MySQL database with the username / db
    #   information from the configuration file.
    D = Db.Mapping(build, Conf.Db)

    # Get the rest of the input variables.
    accno = mrnaSplit(acc)[0]
    version = mrnaSplit(acc)[1]
    
    # Check whether the NM version number is in the database.
    db_version = D.get_NM_version(accno)
    if not db_version :
        return None
    if db_version != version :
        return "Reference sequence version not found. " \
              "Available: %s.%i" % (accno, db_version)
    # Retrieve info on the NM accession number.
    result = D.get_NM_info(accno)
    del D
    
    exonStarts = __sl2il(result[0].split(',')[:-1]) # Get all the exon start 
                                                    # sites.
    exonEnds = __sl2il(result[1].split(',')[:-1])   # Get all the end sites.
    cdsStart = int(result[2])                     # The CDS start.
    cdsEnd = int(result[3])                       # The CDS stop.
    strand = result[4]                            # The orientation.
    
    # Convert the exonStarts and exonEnds lists to an RNA splice sites list and
    #   a CDS splice sites list.
    mRNA = []
    
    CDS = [] 
    if cdsStart != cdsEnd :
        CDS = [cdsStart + 1]              # The counting from 0 conversion.
        CDS.append(cdsEnd)
    #if
    
    for i in range(len(exonStarts)) :
        mRNA.append(exonStarts[i] + 1)    # This is an interbase conversion.
        mRNA.append(exonEnds[i])
    #for
    
    # Convert the strand information to orientation.
    orientation = 1
    if strand == '-' :
        orientation = -1
    
    # Build the crossmapper.
    return Crossmap.Crossmap(mRNA, CDS, orientation)
#makeCrossmap

def mainMapping(build, acc, var, C, O) :
    """
        One of the entry points (called by the HTML publisher).
            
        Arguments:
            build   ; The human genome build.
            acc     ; The NM accession number and version.
            var     ; The variant.
            C       ; Conf object
            O       ; Output object
    
        Returns:
            V       ; ClassSerializer object with types startmain, startoffset,
                      endmain, endoffset, start_g, end_g and mutation type
    """

    Cross = makeCrossmap(build, acc, C)
    if not Cross :
        return None
    if isinstance(Cross, str) :
        return Cross

    # Make a parsetree for the given variation.
    P = Parser.Nomenclatureparser(O)
    parsetree = P.parse("NM_0000:" + var) # This NM number is bogus.
    del P
    if not parsetree :
        return None

    # initiate ClassSerializer object
    V = Mapping()

    if parsetree :
        # Get the coordinates in both c. and g. notation.
        startmain, startoffset, start_g = \
            __getcoords(Cross, parsetree.RawVar.StartLoc.PtLoc, 
                        parsetree.RefType)
        
        # Assume there is no end position given.
        end_g = start_g
        endmain = startmain
        endoffset = startoffset
        
        # If there is an end position, calculate the coordinates.
        if parsetree.RawVar.EndLoc :
            endmain, endoffset, end_g = \
                __getcoords(Cross, parsetree.RawVar.EndLoc.PtLoc, 
                            parsetree.RefType)
        #if
        
        # Assign these values to the Mapping V types.
        V.startmain = startmain
        V.startoffset = startoffset
        V.endmain = endmain
        V.endoffset = endoffset
        V.start_g = start_g
        V.end_g = end_g
        V.mutationType = parsetree.RawVar.MutationType
    #if
    del Cross
    return V
#main_Mapping       

def mainTranscript(build, acc, C, Cross) :
    """
        One of the entry points (called by the HTML publisher).
            
        Arguments:
            build   ; The human genome build.
            acc     ; The full NM accession number (including version).
            C       ; Conf object
    
        Returns:
            T       ; ClassSerializer object with the types trans_start,
                      trans_stop and CDS_stop.

    """

    # Initiate ClassSerializer object
    T = Transcript()
    # Return transcription start, transcription end and
    # CDS stop in c. notation.
    info = Cross.info()
    # Assign these values to the Transcript T types.
    T.trans_start = info[0]
    T.trans_stop  = info[1]
    T.CDS_stop    = info[2]

    return T
#mainTranscript       

def __extractChange(L, var) :
    '''
        Extracts the position and the mutation description (part behind
        the positions) from a (complete) HGVS notation.

        Arguments:
            L   ; output object
            var ; The variant in complete HGVS notation
         
        Returns:
            tuple: position and the part behind the positions
                in HGVS notation
    '''
        
    P = Parser.Nomenclatureparser(L)
    parsetree = P.parse(var)
    del P
    if not parsetree :
        return None
        
    position = parsetree.RawVar.StartLoc.PtLoc.Main # Start position.
    if parsetree.RawVar.Arg1 :
        change = parsetree.RawVar.Arg1
        if parsetree.RawVar.Arg2 :
            change += '>' + parsetree.RawVar.Arg2
    #if
    else :
        change = ''

    if parsetree.RawVar.MutationType not in ('subst', 'del') :
        changeSuffix = parsetree.RawVar.MutationType + change
    else :
        if parsetree.RawVar.MutationType == 'del' :
            changeSuffix = parsetree.RawVar.MutationType
        else :
            changeSuffix = change
    #else
    del parsetree, change

    return position, changeSuffix
#__extractChange


def cTog(build, var, C, D, L) :
    """
       Converts a complete HGVS c. notation into a g. notation

       Arguments:
           build    ; The human genome build.
           var      ; The variant in HGVS c.notation.
            C       ; Conf object
            D       ; A handle to the database.
            L       ; Output object for logging
        
       Returns:
           var_in_g ; The variant in HGVS g. notation (string).
        
    """

    mrnaAcc = var.split(':')[0] # the transcript (NM_..) accession number 
    mut = var.split(':')[1]  # the variant

    # extract the change (e.g. A>T)
    change = __extractChange(L, var)
    if not change :
        return None
    mapped_mrnaAcc = mainMapping(build, mrnaAcc, mut, C, L)
    if not mapped_mrnaAcc :
        return None
    if isinstance(mapped_mrnaAcc, str) :
        return mapped_mrnaAcc
    # you will have to deal with the orientation
    accno = mrnaSplit(mrnaAcc)[0]
    accInfo = D.get_NM_info(accno)
    strand = accInfo[4] # The orientation.
    del accno, accInfo
    
    # construct the variant description
    if str(mapped_mrnaAcc.start_g) != str(mapped_mrnaAcc.end_g) :
        if strand == '+' :
            var_in_g = "g." + str(mapped_mrnaAcc.start_g) \
                    + "_" + str(mapped_mrnaAcc.end_g) + change[1]
        #if
        else :
            var_in_g = "g." + str(mapped_mrnaAcc.end_g) \
                    + "_" + str(mapped_mrnaAcc.start_g) + change[1]
        #else
    #if
    else :
        var_in_g = "g." + str(mapped_mrnaAcc.start_g) + change[1]
    
    del change, mapped_mrnaAcc, strand
    return var_in_g
#cTog


def gToc(build, var, C, D, L) :
    """
        Conversion of (chromosomal) g. to c. numbering
        
        Arguments:
            build   ; The human genome build
            var     ; The (complete) variant in HGVS g. notation
            C       ; Conf object
            D       ; A handle to the database.
            L       ; Output object for logging
            
        Returns:
            list of variants in HGVS c. notation for all transcripts at that
            position
    """

    mut = var.split(':')[1]  # the variant
    accNo = var.split(':')[0] # the chromosomal accession number (NC_)
    # get the chromosome name
    chrom = D.chromName(accNo)
    if not chrom :
        return None

    # extract the start position and change
    change =  __extractChange(L, var)
    if not change :
        return None
    upPos = int(change[0]) - 5000
    downPos = int(change[0]) + 5000

    transcripts = D.get_Transcripts(chrom, upPos, downPos,1)
    HGVS_notations = []

    for mrnaAcc in transcripts:
    
        # Retrieve info on the NM accession number.
        result = D.get_NM_info(mrnaSplit(mrnaAcc)[0])
        strand = result[4] # The orientation.
        del result
    
        mapped_mrnaAcc = mainMapping(build, mrnaAcc, mut, C, L)
        if isinstance(mapped_mrnaAcc, str) :
            return mapped_mrnaAcc
        Cross = makeCrossmap(build, mrnaAcc, C)
        vStart = Cross.tuple2string((mapped_mrnaAcc.startmain,
                                     mapped_mrnaAcc.startoffset))
        vEnd = Cross.tuple2string((mapped_mrnaAcc.endmain,
                                   mapped_mrnaAcc.endoffset))
        info_mrnaAcc = mainTranscript(build, mrnaAcc, C, Cross)
    
        # Use n. or c. numbering
        if info_mrnaAcc.trans_start == "1" and \
            info_mrnaAcc.trans_stop == info_mrnaAcc.CDS_stop :
            numbType = ":n."
        #if
        else :
            numbType = ":c."
    
        if strand == '+' :
            var_in_c = mrnaAcc + numbType + vStart
            if (vStart == vEnd or mapped_mrnaAcc.mutationType == 'subst') :
                var_in_c += change[1]
            else :
                var_in_c += "_" + vEnd + change[1]
        #if
        else :
            # reversed orientation
            var_in_c = mrnaAcc + numbType + vEnd
            if (vStart == vEnd or mapped_mrnaAcc.mutationType == 'subst') :
                var_in_c += change[1]
            else :
                var_in_c += "_" + vStart + change[1]
        #else
        del vStart, vEnd, info_mrnaAcc, mapped_mrnaAcc, numbType
        HGVS_notations.append(var_in_c)
    #for
    return HGVS_notations
#gToC
