#!/usr/bin/python

"""
Search for an NM number in the MySQL database, if the version number
matches, get the start and end positions in a variant and translate these
positions to I{g.} notation if the variant is in I{c.} notation and vice versa.
    - If no end position is present, the start position is assumed to be the
      end position.
    - If the version number is not found in the database, an error message is
      generated and a suggestion for an other version is given.
    - If the reference sequence is not found at all, an error is returned.
    - If no variant is present, the transcription start and end and CDS end
      in I{c.} notation is returned.
    - If the variant is not accepted by the nomenclature parser, a parse error
      will be printed.

@requires: sys
@requires: Modules.Db
@requires: Modules.Crossmap
@requires: Modules.Parser
@requires: Modules.Output
@requires: Modules.Config
@requires: Modules.Mapper

@todo: documentation
"""

import sys                   # argv
from mutalyzer import Db       # Db(), get_NM_version(), get_NM_info()
from mutalyzer import Crossmap # Crossmap(), g2x(), x2g(), main2int(),
                             # offset2int(), info()
from mutalyzer import Parser   # Nomenclatureparser(), parse()
from mutalyzer import Output   # Output(), LogMsg()
from mutalyzer import Config
from mutalyzer import Mapper

def __sl2il(l) :
    """
    Convert a list of strings to a list of integers.
    
    @arg l: A list of strings
    @type l: list

    @return: A list of integers
    @rtype: list
    """

    for i in range(len(l)) :
        l[i] = int(l[i])
    return l
#__sl2il

def __getcoords(C, Loc, Type) :
    """
    Return main, offset and g positions given either a position in
    I{c.} or in I{g.} notation.

    @arg C: A crossmapper
    @type C: object
    @arg Loc: Either a location in I{g.} or I{c.} notation
    @type Loc: object
    @arg Type: The reference type
    @type Type: character
    
    @return: triple:
        - 0 ; Main coordinate in I{c.} notation
        - 1 ; Offset coordinate in I{c.} notation
        - 2 ; Position in I{g.} notation
    @rtype: triple (integer, integer, integer)
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

    Returns:
      - start_main   ; The main coordinate of the start position in I{c.}
                       (non-star) notation.
      - start_offset ; The offset coordinate of the start position in I{c.}
                       notation (intronic position).
      - end_main     ; The main coordinate of the end position in I{c.}
                       (non-star) notation.
      - end_offset   ; The offset coordinate of the end position in I{c.}
                       notation (intronic position).
      - start_g      ; The I{g.} notation of the start position.
      - end_g        ; The I{g.} notation of the end position.
      - type         ; The mutation type.

    Returns (alternative):
      - trans_start  ; Transcription start in I{c.} notation.
      - trans_stop   ; Transcription stop in I{c.} notation.
      - CDS_stop     ; CDS stop in I{c.} notation.

    @arg LOVD_ver: The LOVD version (ignored for now)
    @type LOVD_ver: string
    @arg build: The human genome build
    @type build: string
    @arg acc: The NM accession number and version
    @type acc: string
    @arg var: The variant, or empty
    @type var: string
    
    @return: 
    @rtype: 
    """

    C = Config.Config()
    O = Output.Output(__file__, C.Output)

    O.addMessage(__file__, -1, "INFO",
                 "Received %s:%s (LOVD_ver %s, build %s)" % (acc, var,
                 LOVD_ver, build))

    Converter = Mapper.Converter(build, C, O)

    #V = Mapper.makeParsetree(O, Cross, var)

    # If no variant is given, return transcription start, transcription end and
    #   CDS stop in c. notation.
    if var :
        ret = Converter.mainMapping(acc, var)
        #for i in Converter.crossmap.RNA :
        #    print i, Converter.crossmap.g2c(i)
    else :
        ret = Converter.giveInfo(acc)
        if ret:
            return "%i\n%i\n%i" % ret

    if not getattr(ret, "startmain", None) :
        output = O.getOutput("LOVDERR")
        if output:
            return output[0]
        else:
            #print "\n".join(O.getMessages())
            return "Unknown error occured"

    O.addMessage(__file__, -1, "INFO",
                 "Finished processing %s:%s (LOVD_ver %s, build %s)" % (acc,
                 var, LOVD_ver, build))
    del O, C
    # And return the output.
    return "%i\n%i\n%i\n%i\n%i\n%i\n%s" % (ret.startmain, ret.startoffset,
        ret.endmain, ret.endoffset, ret.start_g, ret.end_g, ret.mutationType)

#main

if __name__ == "__main__" :
    print main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
