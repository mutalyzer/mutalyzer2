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

# 19-04-2010 by gerard
from soaplib.serializers.primitive import String, Integer
from soaplib.serializers.clazz import ClassSerializer


class Complex2(ClassSerializer) :
    '''
        Extended ClassSerializer object with mixed types of attributes
        
        Attributes:
            trans_start ; Define the type of trans_start
            trans_stop  ; Define the type of trans_stop
            CDS_stop    ; Define the type of CDS_stop
    '''
    class types :
        trans_start = Integer
        trans_stop  = Integer
        CDS_stop    = Integer
    #types
#Complex
Complex2.typecode = TC.Struct(Complex2, [ TC.Integer('trans_start'),
                                        TC.Integer('trans_stop'),
                                        TC.Integer('CDS_stop') ], 'Complex2')

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

def __process(LOVD_ver, build, acc, Conf, O) :
    # Make a connection to the MySQL database with the username / db
    #   information from the configuration file.
    Database = Db.Db(Conf, "local") # Open the database.
    
    # Get the rest of the input variables.
    accno = acc
    version = 0
    if '.' in acc :
        accno = acc.split('.')[0]        # The NM accession number.
        version = int(acc.split('.')[1]) # The version of the accession 
                                         #   number.
    
    # Check whether the NM version number is in the database.
    db_version = Database.get_NM_version(accno)
    if not db_version :
        O.LogMsg(__file__, "EARG %s" % accno)
        raise Fault(Fault.Client, "EARG",
            detail = "The accno argument %s was not a valid "\
                     "NM accession number." % accno)
        return
    #if
    if db_version != version :
        O.ErrorMsg(__file__, 
            "Reference sequence version not found. Available: %s.%i" % (
            accno, db_version))
        raise Fault(Fault.Client, "EARG",
            detail = "The Reference sequence version (%s.%i) was not found. "\
                     "Available: %s.%i" % (accno, version, accno, db_version))
        return
    #if

    # Retrieve info on the NM accession number.
    result = Database.get_NM_info(accno)
    del Database
    
    exonStarts = __sl2il(result[0].split(',')[:-1]) # Get all the exon start 
                                                    # sites.
    exonEnds = __sl2il(result[1].split(',')[:-1])   # Get all the end sites.
    cdsStart = int(result[2])                     # The CDS start.
    cdsEnd = int(result[3])                       # The CDS stop.
    strand = result[4]                            # The orientation.
    
    # Convert the exonStarts and exonEnds lists to an RNA splice sites list and
    #   a CDS splice sites list.
    mRNA = []
    
    CDS = [cdsStart + 1]                  # The counting from 0 conversion.
    CDS.append(cdsEnd)
    
    for i in range(len(exonStarts)) :
        mRNA.append(exonStarts[i] + 1)    # This is an interbase conversion.
        mRNA.append(exonEnds[i])
    #for
    
    # Convert the strand information to orientation.
    orientation = 1
    if strand == '-' :
        orientation = -1
    
    # Build the crossmapper.
    Cross = Crossmap.Crossmap(mRNA, CDS, orientation)
    
    # Return transcription start, transcription end and
    # CDS stop in c. notation.
    T = Complex2()
    info = Cross.info()
    T.trans_start = info[0]
    T.trans_stop  = info[1]
    T.CDS_stop    = info[2]
    print "%i\n%i\n%i" % info
    return T

#__map
def main(LOVD_ver, build, acc, var) :
    """
        The entry point (called by the HTML publisher).
            
        Arguments:
            LOVD_ver ; The LOVD version (ignored for now).
            build    ; The human genome build (ignored for now, hg19 assumed).
            acc      ; The NM accession number and version.
    
        Returns:
            trans_start  ; Transcription start in c. notation.
            trans_stop   ; Transcription stop in c. notation.
            CDS_stop     ; CDS stop in c. notation.

        On error an exception is raised:
            detail       ; Human readable description of the error.
            faultstring: ; A code to indicate the type of error.
                EARG   ; The argument was not valid.
                ERANGE ; An invalid range was given.
    """

    Conf = Config.Config() # Read the configuration file.
    O = Output.Output(Conf, __file__)

    O.LogMsg(__file__, "Received %s:%s (LOVD_ver %s, build %s)" % (
        acc, var, LOVD_ver, build))

    result = __process(LOVD_ver, build, acc, Conf, O)

    O.LogMsg(__file__, "Finished processing %s:%s (LOVD_ver %s, build %s)" % (
        acc, var, LOVD_ver, build))
    del O
    del Conf
    return result
#main        

if __name__ == "__main__" :
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

        

