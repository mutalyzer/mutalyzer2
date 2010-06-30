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

import sys                   # argv
from Modules import Db       # Db(), get_NM_version(), get_NM_info()
from Modules import Crossmap # Crossmap(), g2x(), x2g(), main2int(), 
                             # offset2int(), info()
from Modules import Parser   # Nomenclatureparser(), parse()
from Modules import Output   # Output(), LogMsg()
from Modules import Config
from Modules import Mapper

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

def main(LOVD_ver, build, acc, var) :
    """
        The entry point (called by the HTML publisher).
            
        Arguments:
            LOVD_ver ; The LOVD version (ignored for now).
            build    ; The human genome build.
            acc      ; The NM accession number and version.
            var      ; The variant, or empty.
    
        Returns:
            start_main   ; The main coordinate of the start position in c. 
                           (non-star) notation.
            start_offset ; The offset coordinate of the start position in c. 
                           notation (intronic position).
            end_main     ; The main coordinate of the end position in c. 
                           (non-star) notation.
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

    C = Config.Config()
    O = Output.Output(__file__, C.Output)

    O.addMessage(__file__, -1, "INFO", 
                 "Received %s:%s (LOVD_ver %s, build %s)" % (acc, var,
                 LOVD_ver, build))

    Cross = Mapper.makeCrossmap(build, acc, C)

    # If no variant is given, return transcription start, transcription end and
    #   CDS stop in c. notation.
    if not var :
        info = Cross.info()
        print "%i\n%i\n%i" % info
        return
    #if
    
    #V = Mapper.makeParsetree(O, Cross, var)
    V = Mapper.mainMapping(build, acc, var, C, O)

    O.addMessage(__file__, -1, "INFO", 
                 "Finished processing %s:%s (LOVD_ver %s, build %s)" % (acc,
                 var, LOVD_ver, build))
    del O, C, Cross
    # And return the output.
    print "%i\n%i\n%i\n%i\n%i\n%i\n%s" % (V.startmain, V.startoffset, V.endmain,
                                          V.endoffset, V.start_g, V.end_g, 
                                          V.mutationType)
        
#main        

if __name__ == "__main__" :
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
