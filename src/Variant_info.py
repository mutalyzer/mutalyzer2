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

    Arguments:
        $1 ; The LOVD version (ignored for now).
        $2 ; The human genome build (ignored for now, hg19 assumed).
        $3 ; The NM accession number and version.
        $4 ; The variant, or empty.

    Returns:
        start_main   ; The main coordinate of the start position in c. 
                       (non-star) notation.
        start_offset ; The offset coordinate of the start position in c. 
                       notation (intronic position).
        end_main     ; The main coordinate of the end position in c. (non-star)
                       notation.
        end_offset   ; The offset coordinate of the end position in c. 
                       notation (intronic position).
        start_g      ; The g. notation of the start position.
        end_g        ; The g. notation of the end position.
        type         ; The mutation type.

    Returns (alternative):
        trans_start  ; Transcription start in c. notation.
        trans_stop   ; Transcription stop in c. notation.
        CDS_stop     ; CDS stop in c. notation.
"""

import sys      # argv
import Config   # Config()
import Db       # Db(), get_NM_version(), get_NM_info()
import Crossmap # Crossmap(), g2x(), x2g(), main2int(), offset2int(), info()
import Parser   # Nomenclatureparser(), parse()
import Output   # Output(), LogMsg()

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


def __process(LOVD_ver, build, acc, var, Conf, O) :
    # Make a connection to the MySQL database with the username / db
    #   information from the configuration file.
    Database = Db.Db(Conf) # Open the database.
    
    # Get the rest of the input variables.
    accno = acc.split('.')[0]        # The NM accession number.
    version = int(acc.split('.')[1]) # The version of the accession 
                                     #   number.
    
    # Check whether the NM version number is in the database.
    db_version = Database.get_NM_version(accno)
    if not db_version :
        O.ErrorMsg(__file__, "Reference sequence not found.")
        return
    #if
    if db_version != version :
        O.ErrorMsg(__file__, 
            "Reference sequence version not found. Available: %s.%i" % (
            accno, db_version))
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
    
    # If no variant is given, return transcription start, transcription end and
    #   CDS stop in c. notation.
    if not var :
        info = Cross.info()
        print "%i\n%i\n%i" % info
        return
    #if
    
    # Make a parsetree for the given variation.
    P = Parser.Nomenclatureparser(O)
    parsetree = P.parse("NM_0000:" + var) # This NM number is bogus.
    del P

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
        
        # And return the output.
        print "%i\n%i\n%i\n%i\n%i\n%i\n%s" % (startmain, startoffset, endmain, 
                                              endoffset, start_g, end_g, 
                                              parsetree.RawVar.MutationType)
    #if
    del Cross
#__process

def main(LOVD_ver, build, acc, var) :
    Conf = Config.Config() # Read the configuration file.
    O = Output.Output(Conf, __file__)

    O.LogMsg(__file__, "Received %s:%s (LOVD_ver %s, build %s)" % (
        acc, var, LOVD_ver, build))

    __process(LOVD_ver, build, acc, var, Conf, O)

    O.LogMsg(__file__, "Finished processing %s:%s (LOVD_ver %s, build %s)" % (
        acc, var, LOVD_ver, build))
    del O
    del Conf
#main        

if __name__ == "__main__" :
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
