#!/usr/bin/python

"""
    Search for an NM number in the MySQL database, if the version number 
    matches, get the start and end positions in a variant and translate these
    positions to g. notation.

    - If no end position is present, the start position is assumed to be the
      end position. 
    - If the version number is not found in the database, an error message is
      generated and a suggestion for an other version is given.
    - If the reference sequence is not found at all, an error is returned.
    - If the variant is not accepted by the nomenclature parser, a parse error
      will be printed.

    Arguments:
        $1 ; The LOVD version (ignored for now).
        $2 ; The human genome build (ignored for now, hg19 assumed).
        $3 ; The NM accession number and version.
        $4 ; The variant.

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
"""

import sys      # argv
import Config   # Config()
import Db       # Db(), get_NM_version(), get_NM_info()
import Crossmap # Crossmap(), star2abs() off2int(), c2g()
import Parser   # Nomenclatureparser(), parse()

def sl2il(l) :
    """
        Convert a list of strings to a list of integers.

        Arguments: l ; A list of strings.

        Returns: list ; A list of integers.
    """

    for i in range(len(l)) :
        l[i] = int(l[i])
    return l
#sl2il

# Make a connection to the MySQL database with the username / db information
#   from the configuration file.
Conf = Config.Config() # Read the configuration file.
Database = Db.Db(Conf) # Open the database.

# Get all the input variables.
LOVD_ver = sys.argv[1]                 # The LOVD version (not used).
build = sys.argv[2]                    # The human genome build (hg 19 assumed).
accno = sys.argv[3].split('.')[0]        # The NM accession number.
version = int(sys.argv[3].split('.')[1]) # The version of the accession number.
variant = sys.argv[4]                    # The variant in HGVS notation.

# Check whether the NM version number is in the database.
db_version = Database.get_NM_version(accno)
if not db_version :
    print "Error: Reference sequence not found."
    exit()
#if
if db_version != version :
    print "Error: Reference sequence version not found. Available: %s.%i" % (
          accno, db_version)
    exit()
#if

# Retrieve info on the NM accession number.
result = Database.get_NM_info(accno)

exonStarts = sl2il(result[0].split(',')[:-1]) # Get all the exon start sites.
exonEnds = sl2il(result[1].split(',')[:-1])   # Get all the exon end sites.
cdsStart = int(result[2])                     # The CDS start.
cdsEnd = int(result[3])                       # The CDS stop.
strand = result[4]                            # The orientation.

# Convert the exonStarts and exonEnds lists to an RNA splice sites list and
#   a CDS splice sites list.
mRNA = []
CDS = [cdsStart + 1]                  # The counting from 0 conversion.
for i in range(len(exonStarts)) :
    mRNA.append(exonStarts[i] + 1)    # This is an interbase conversion.
    mRNA.append(exonEnds[i])
    if exonStarts[i] >= cdsStart and exonStarts[i] <= cdsEnd :
        CDS.append(exonStarts[i] + 1) # Also an interbase conversion.
    if exonEnds[i] >= cdsStart and exonEnds[i] <= cdsEnd :
        CDS.append(exonEnds[i])
#for
CDS.append(cdsEnd)

# Convert the strand information to orientation.
orientation = 1
if strand == '-' :
    orientation = -1

# Now we can build the crossmapper.
Cross = Crossmap.Crossmap(mRNA, CDS, orientation)

# Make a parsetree for the given variation.
P = Parser.Nomenclatureparser()
parsetree = P.parse("NM_0000:" + variant) # This NM number is bogus.

# Get the main and offset coordinates of the start position in non-star 
#   notation.
startmain = Cross.star2abs(parsetree.RawVar.StartLoc.PtLoc.MainSgn,
                           parsetree.RawVar.StartLoc.PtLoc.Main)
startoffset =  Cross.off2int(parsetree.RawVar.StartLoc.PtLoc.OffSgn,
                             parsetree.RawVar.StartLoc.PtLoc.Offset)
print "%i\n%i" % (startmain, startoffset)

# Do the same for the end position if it is present, otherwise repeat the
#   start position.
if parsetree.RawVar.EndLoc :
    print Cross.star2abs(parsetree.RawVar.EndLoc.PtLoc.MainSgn,
                         parsetree.RawVar.EndLoc.PtLoc.Main)
    print Cross.off2int(parsetree.RawVar.EndLoc.PtLoc.OffSgn,
                        parsetree.RawVar.EndLoc.PtLoc.Offset)
#if                        
else :
    print "%i\n%i" % (startmain, startoffset)

# Calculate the g. notation of the start position.
start_g = Cross.c2g(parsetree.RawVar.StartLoc.PtLoc.MainSgn,
                    parsetree.RawVar.StartLoc.PtLoc.Main,
                    parsetree.RawVar.StartLoc.PtLoc.OffSgn,
                    parsetree.RawVar.StartLoc.PtLoc.Offset)
print start_g

# Do the same for the end position if it is present, otherwise repeat the
#   start position.
if parsetree.RawVar.EndLoc :
    print Cross.c2g(parsetree.RawVar.EndLoc.PtLoc.MainSgn,
                    parsetree.RawVar.EndLoc.PtLoc.Main,
                    parsetree.RawVar.EndLoc.PtLoc.OffSgn,
                    parsetree.RawVar.EndLoc.PtLoc.Offset)
else :
    print start_g

# Print the variant type.
print parsetree.RawVar.MutationType
