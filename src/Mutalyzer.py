#!/usr/bin/python

"""
    The nomenclature checker.
"""

import sys
import math
import types
import Bio

import Bio.Seq
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import seq3
from Bio import Restriction

from Modules import Retriever
from Modules import GenRecord
from Modules import Crossmap
from Modules import Parser
from Modules import Db
from Modules import Mutator
from Modules import Output
from Modules import Config

from operator import itemgetter, attrgetter

#TODO: SET TO FALSE DEBUG FLAG
DEBUG = False

def __formatRange(pos1, pos2) :
    """
        Simplify a range to one position when applicable.

        Arguments:
            pos1 ; First coordinate of a range.
            pos2 ; Second coordinate of a range.

        Returns:
            string ; pos1_pos2 in case of a real range, pos1 otherwise.
    """

    if pos1 == pos2 :
        return str(pos1)
    return "%i_%i" % (pos1, pos2)
#__formatRange

def __intronicPosition(Loc) :
    """
        Check whether a location is intronic.

        Arguments:
            Loc ; A location from the Parser module.

        Returns:
            boolean ; True if the location is intronic, False otherwise.
    """

    if not Loc :
        return False
    if not Loc.PtLoc :
        return False
    if not Loc.PtLoc.Offset :
        return False
    return True
#__intronicPosition

def __checkIntronPosition(main, offset, transcript) :
    """
        Check whether a c. position is really in an intron: The main coordinate
        must be a splice site and the offset coordinate must have the correct
        sign.

        Arguments:
            main       ; Main coordinate of the position.
            offset     ; Offset coordinate of the position.
            transcript ; Transcript under scrutiny.

        Returns:
            boolean ; True if the combination (main, offset) is valid for this
                      transcript. False otherwise.
    """

    main_g = transcript.CM.x2g(main, 0)
    rnaList = transcript.CM.RNA

    if offset :
        #print main_g, offset, rnaList
        orientedOffset = offset * transcript.CM.orientation
        if main_g in rnaList :          # The main coordinate is a splice site.
            if rnaList.index(main_g) % 2 == 0 : # Splice donor.
                if orientedOffset > 0 :         # So the sign must be '+'.
                    return False
            else :                              # Splice acceptor.
                if orientedOffset < 0 :         # So the sign must be '-'.
                    return False
        #if
        else :
            return False
    #if

    return True
#__checkIntronPosition            

def __roll(ref, start, stop) :
    """
        Determine the variability of a variant by looking at cyclic
        permutations. Not all cyclic permutations are tested at each time, it
        is sufficient to check ``aW'' if ``Wa'' matches (with ``a'' a letter,
        ``W'' a word) when rolling to the left for example.

        Arguments:
            ref   ; A reference sequence.
            start ; Start position of the pattern in the reference sequence.
            stop  ; End position of the pattern in the reference sequence.

        Returns:
            tuple: 
                left  ; Amount of positions that the pattern can be shifted to 
                        the left.
                right ; Amount of positions that the pattern can be shifted to
                        the right.
    """

    pattern = ref[start - 1:stop] # Extract the pattern.
    patternLength = len(pattern)  

    # Keep rolling to the left as long as a cyclic permutation matches.
    minimum = start - 2
    j = patternLength - 1
    while minimum > -1 and ref[minimum] == pattern[j % patternLength] :
        j -= 1
        minimum -= 1
    #while

    # Keep rolling to the right as long as a cyclic permutation matches.
    maximum = stop
    j = 0
    while maximum < len(ref) and ref[maximum] == pattern[j % patternLength] :
        j += 1
        maximum += 1
    #while

    return start - minimum - 2, maximum - stop
#__roll

def __palinsnoop(string) :
    """
        Check a sequence for a reverse-complement-palindromic prefix (and
        suffix). If one is detected, return the length of this prefix. If the
        string equals its reverse complement, return -1.
        
        Arguments:
            string ; A nucleotide sequence.

        Returns:
            integer ; The number of elements that are palindromic or -1 if the
                      string is a ``palindrome''.
    """

    revcomp = Bio.Seq.reverse_complement(string)

    for i in range(int(math.ceil(len(string) / 2.0))) :
        if string[i] != revcomp[i] :
            return i # The first i elements are ``palindromic''.
    return -1        # Perfect ``palindrome''.
#__palinsnoop

def __bprint(s, O, where) :
    # FIXME obsoleted function (replaced by __bprint2()), but still used.
    """
    """

    if not s :
        return

    block = 10
    line = 6 * block

    m = int(math.floor(math.log(len(s), 10)) + 1)
    o = 1
    output = "%s " % str(o).rjust(m)
    for i in range(0, len(s), block) :
        output += ' ' + s[i:i + block]
        if not (i + block) % line and i + block < len(s) :
            o += line
            O.addOutput(where, output)
            output = "%s " % str(o).rjust(m)
        #if
    #for
    O.addOutput(where, output)
#__bprint

def __insertTag(s, pos1, pos2, tag1, tag2) :
    """
        Insert two tags (tag1 and tag2) in string s at positions pos1 and pos2
        respectively if the positions are within the length of s. If not,
        either insert one tag or do nothing. If pos1 equals pos2, don't do
        anything either.

        Arguments:
            s    ; A sequence.
            pos1 ; Position of tag1.
            pos2 ; Position of tag2.
            tag1 ; Content of tag1.
            tag2 ; Content of tag2.

        Returns:
            string ; The original sequence, or a sequence with eiter tag1,
                     tag2 or both tags inserted.
    """

    output = s
    block = len(s)

    if pos1 != pos2 :               # Only do something if pos1 != pos2.
        if 0 <= pos1 < block :
            output = output[:pos1] + tag1 + output[pos1:] # Insert tag1.
        if 0 <= pos2 < block :
            output = output[:-(block - pos2)] + tag2 + \
                     output[-(block - pos2):]             # Insert tag2.
    #if

    return output
#__insertTag

def __bprint2(s, pos1, pos2, O, where) :
    """
        Make a fancy representation of a protein and put it in the Output
        object under the name ``where''.

        Arguments:
            s     ; A protein sequence.
            pos1  ; First position to highlight.
            pos2  ; Last position to highlight.
            O     ; The Output object.
            where ; Location in the Output object to store the representation.
    """

    if not s :
        return

    block = 10       # Each block consists of 10 amino acids.
    line = 6 * block # Each line consists of 6 blocks.

    tag1 = "<b style=\"color:#FF0000\">" # Use this tag for highlighting.
    tag2 = "</b>"                        # And this one to end highlighting.

    # The maximum length for positions is the 10_log of the length of the
    #   protein.
    m = int(math.floor(math.log(len(s), 10)) + 1) 
    o = 1
    output = "%s " % str(o).rjust(m)     # Add the first position.
    for i in range(0, len(s), block) :   # Add the blocks.
        output += ' ' + __insertTag(s[i:i + block], pos1 - i,
                                    pos2 - i, tag1, tag2)
        if not (i + block) % line and i + block < len(s) :
            o += line                    # One line done.
            O.addOutput(where, output)   # Add it to the output.
            # And add the next line (while escaping any potential highlighting).
            output = \
                "<tt style = \"color:000000;font-weight:normal\">%s</tt> " % \
                str(o).rjust(m)
        #if
    #for
    O.addOutput(where, output)
#__bprint2

def __PtLoc2main(Loc) :
    """
        Convert the main coordinate in a location (from the Parser) to an
        integer.

        Arguments:
            Loc ; A location.

        Returns:
            integer ; Integer representation of the main coordinate.
    """

    main = int(Loc.Main)
    if Loc.MainSgn == '-' :
        return -main

    return main
#__PtLoc2main

def __PtLoc2offset(Loc) :
    """
        Convert the offset coordinate in a location (from the Parser) to an 
        integer.

        Arguments:
            Loc ; A location.

        Returns;
            integer ; Integer representation of the offset coordinate.
    """

    if Loc.Offset :
        if Loc.Offset == '?' : # This is highly debatable.
            return 0
        offset = int(Loc.Offset)
        if Loc.OffSgn == '-' :
            return -offset
        return offset
    #if

    return 0
#__PtLoc2offset

def __splice(string, splice_sites) :
    """
        Construct the transcript or the coding sequence from a record and
        a list of splice sites.

        Arguments:
            record       ; A GenBank record (see the BioPython documentation).
            splice_sites ; A list of even length of integers.

        Returns:
            String ; The concatenation of slices from the sequence that is
                     present in the GenBank record.
    """

    transcript = ""

    for i in range(0, len(splice_sites), 2) :
        transcript += string[splice_sites[i] - 1:splice_sites[i + 1]]

    return transcript
#__splice

def __nsplice(string, splice_sites, CDS, orientation) :
    #FIXME document this.
    """
    """

    transcript = ""
    if orientation == 1 :
        for i in range(0, len(splice_sites), 2) :
            if CDS[0] >= splice_sites[i] and CDS[0] <= splice_sites[i + 1] :
                transcript += string[CDS[0] - 1:splice_sites[i + 1]]
            else :
                if splice_sites[i] > CDS[0] :
                    transcript += \
                        string[splice_sites[i] - 1:splice_sites[i + 1]]
        #for
    #if
    else :
        for i in range(0, len(splice_sites), 2) :
            if CDS[1] >= splice_sites[i] and CDS[1] <= splice_sites[i + 1] :
                transcript += string[splice_sites[i] - 1: CDS[1]]
            else :
                if splice_sites[i] < CDS[1] :
                    transcript += \
                        string[splice_sites[i] - 1:splice_sites[i + 1]]
        #for
    #else

    return transcript
#__nsplice

def __cdsLen(splice_sites) :
    """
        Calculate the length of a CDS.

        Arguments:
            splice_sites ; The coordinates of the CDS including internal splice
                           sites.
                           
        Returns:                           
            integer ; Length of the CDS.
    """

    l = 0

    for i in range(0, len(splice_sites), 2) :
        l += splice_sites[i + 1] - splice_sites[i] + 1
    return l
#__cdsLen

def __checkDNA(arg) :
    """
        Check whether a string is a DNA string.

        Arguments:
            arg ; Any string.

        Returns:
            boolean ; True if the string is a DNA string, False otherwise.
    """

    for i in str(arg) :
        if not i in IUPAC.unambiguous_dna.letters :
            return False
    return True
#__checkDNA

def __checkOptArg(ref, p1, p2, arg, O) :
    """
        Do several checks for the optional argument of a variant.


        Arguments:
            ref ; The reference sequence.
            p1  ; Start position of the variant.
            p2  ; End position of the variant.
            arg ; The optional argument.
            O   ; The Output object.

        Returns:
            boolean ; True if the optional argument is correct, False otherwise.
    """

    if arg : # The argument is optional, if it is not present, it is correct.
        if arg.isdigit() :         # If it is a digit (3_9del7 for example),
            length = int(arg)      #   the digit must be equal to the length
            interval = p2 - p1 + 1 #   of the given range.
            if length != interval :
                O.addMessage(__file__, 3, "EARGLEN",
                    "The length (%i) differed from that of the range (%i)." % (
                    length, interval))
                return False
            #if
        #if
        else :
            if not __checkDNA(arg) : # If it is not a digit, it muse be DNA.
                O.addMessage(__file__, 4, "ENODNA",
                    "Invalid letters in argument.")
                return False
            #if
            # And the DNA must match the reference sequence.
            ref_slice = str(ref[p1 - 1:p2])
            if ref_slice != str(arg) : # FIXME more informative.
                O.addMessage(__file__, 3, "EREF",
                    "%s not found at position %s, found %s instead." % (
                    arg, __formatRange(p1, p2), ref_slice))
                return False
            #if
        #else
    #if
    return True
#__checkOptArg

def __lcp(str1, str2) :
    """
        Calculate the length of the longest common prefix of two strings.

        Arguments:
            str1 ; The first string.
            str2 ; The second string.

        Returns:
            integer ; The length of the longest common prefix of str1 and str2.
    """

    pos = 0
    s1l = len(str1) # Use the lengths to make sure we don't exceed the length
    s2l = len(str2) # of the strings.

    while pos < s1l and pos < s2l and str1[pos] == str2[pos] :
        pos += 1

    return pos
#__lcp

def __lcs(str1, str2) :
    """
        Calculate the length of the longest common suffix of two strings.

        Arguments:
            str1 ; The first string.
            str2 ; The second string.

        Returns:
            integer ; The length of the longest common suffix of str1 and str2.
    """

    t1 = str1[::-1] # Invert str1.
    t2 = str2[::-1] # Invert str2.

    # The lcp of the two inverted strings is the lcs of the original strings.
    return __lcp(t1, t2)
#__lcs

def findInFrameDescription(str1, str2) :
    """
        Give a description of an inframe difference of two proteins. Also give
        the position at which the proteins start to differ and the positions at
        which they are the same again.

        Arguments:
            str1 ; The original protein.
            str2 ; The mutated protein.

        Retuns:
            vector:
                string  ; Protein description of the change.
                integer ; Start position of the change.
                integer ; End position of the change in the first protein.
                integer ; End position of the change in the second protein.
    """

    # Nothing happened.
    if str1 == str2 :
        return ("p.(=)", 0, 0, 0)

    lcp = __lcp(str1, str2)
    lcs = __lcs(str1[lcp:], str2[lcp:])
    str1_end = len(str1) - lcs
    str2_end = len(str2) - lcs

    # Insertion / Duplication / Extention.
    if not str1_end - lcp :
        if len(str1) == lcp :
            return ("p.(*%i%sext*%i)" % (len(str1) + 1, seq3(str2[len(str1)]),
                abs(len(str1) - len(str2))), len(str1), len(str1), len(str2))
        inLen = str2_end - lcp

        if lcp - inLen >= 0 and str1[lcp - inLen:lcp] == str2[lcp:str2_end] :
            if inLen == 1 :
                return ("p.(%s%idup)" % (seq3(str1[lcp - inLen]),
                    lcp - inLen + 1), 
                        lcp, lcp, lcp + 1)
            return ("p.(%s%i_%s%idup)" % (seq3(str1[lcp - inLen]),
                lcp - inLen + 1, seq3(str1[lcp - 1]), lcp), lcp, lcp, 
                lcp + inLen)
        #if
        return ("p.(%s%i_%s%iins%s)" % (seq3(str1[lcp - 1]), lcp,
            seq3(str1[lcp]), lcp + 1, seq3(str2[lcp:str2_end])), lcp, lcp, 
            str2_end)
    #if

    # Deletion / Inframe stop.
    if not str2_end - lcp :
        if len(str2) == lcp :
            return ("p.(%s%i*)" % (seq3(str1[len(str2)]), len(str2) + 1), 
                0, 0, 0)

        if lcp + 1 == str1_end :
            return ("p.(%s%idel)" % (seq3(str1[lcp]), lcp + 1), 
                lcp, lcp + 1, lcp)
        return ("p.(%s%i_%s%idel)" % (seq3(str1[lcp - 1]), lcp + 1,
            seq3(str1[str1_end - 1]), str1_end), lcp, str1_end, lcp)
    #if

    # Substitution.
    if str1_end == str2_end and str1_end == lcp + 1 :
        return ("p.(%s%i%s)" % (seq3(str1[lcp]), lcp + 1, seq3(str2[lcp])), 
            lcp, lcp + 1, lcp + 1)

    # InDel.
    if lcp + 1 == str1_end :
        return ("p.(%s%idelins%s)" % (seq3(str1[lcp]), lcp + 1,
            seq3(str2[lcp:str2_end])), lcp, lcp + 1, str2_end)
    return ("p.(%s%i_%s%idelins%s)" % (seq3(str1[lcp]), lcp + 1,
        seq3(str1[str1_end - 1]), str1_end, seq3(str2[lcp:str2_end])), lcp, 
        str1_end, str2_end)
#findInFrameDescription

def findFrameShift(str1, str2) :
    """
        Give the description of an out of frame difference between two
        proteins. Give a description of an inframe difference of two proteins.
        Also give the position at which the proteins start to differ and the
        end positions (to be compatible with the findInFrameDescription()
        function).

        Arguments:
            str1 ; The original protein.
            str2 ; The mutated protein.

        Retuns:
            vector:
                string  ; Protein description of the change.
                integer ; Start position of the change.
                integer ; End position of the first protein.
                integer ; End position of the second protein.
    """

    lcp = __lcp(str1, str2)

    if lcp == len(str2) : # NonSense mutation.
        return ("p.(%s%i*)" % (seq3(str1[lcp]), lcp + 1), lcp, len(str1), lcp)
    if lcp == len(str1) :
        return ("p.(*%i%sext*%i)" % (len(str1) + 1, seq3(str2[len(str1)]),
            abs(len(str1) - len(str2))), len(str1), len(str1), len(str2))
    return ("p.(%s%i%sfs*%i)" % (seq3(str1[lcp]), lcp + 1, seq3(str2[lcp]),
        len(str2) - lcp + 1), lcp, len(str1), len(str2))
#findFrameShift

def __toProtDescr(CDSStop, orig, trans) :
    """
        Wrapper function for the findInFrameDescription() and findFrameShift()
        functions. It uses the value CDSStop to decide which one to call.

        Arguments:
            CDSStop ; Position of the stop codon in c. notation (CDS length).
            orig    ; The original protein.
            trans   ; The mutated protein.

        Retuns:
            vector:
                string  ; Protein description of the change.
                integer ; Start position of the change.
                integer ; End position of the change in the first protein.
                integer ; End position of the change in the second protein.
    """

    if CDSStop % 3 :
        ret = findFrameShift(str(orig), str(trans))
    else :
        ret = findInFrameDescription(str(orig), str(trans))
    if str(orig[0]) != str(trans[0]) :         # Mutation in start codon.
        return ("p.?", ret[1], ret[2], ret[3])
    return ret
#__toProtDescr

def __trim2(str1, str2) :
    """
        Given two strings, trim the lcp and the lcs.
        
        Arguments:
            str1 ; A string.
            str2 ; An other string.

        Returns:
            tuple:
                string ; Trimmed version of str1.
                string ; Trimmed version of str2.
    """

    lcp = __lcp(str1, str2)
    lcs = __lcs(str1[lcp:], str2[lcp:])
    return str1[lcp:len(str1) - lcs], str2[lcp:len(str2) - lcs], lcp, lcs
#__trim2

def __rangeToC(M, g1, g2) :
    # FIXME apparently obsolete.
    """
        Convert a genomic range to a CDS oriented range.

        Arguments:
            M  ;
            g1 ;
            g2 ;

        Returns:
            tuple:
                string ;
                string ;
    """

    if M.orientation == -1 :
        return M.g2c(g2), M.g2c(g1)
    return M.g2c(g1), M.g2c(g2)
#__rangeToC

def _createBatchOutput(O):
    #TODO More documentation.
    """
        Format the results to a batch output.

        Filter the mutalyzer output

    """
    goi, toi = O.getOutput("geneSymbol")[-1] # Two strings [can be empty]
    tList   = []                             # Temporary List
    tDescr  = []                             # Temporary Descr

    reference = O.getOutput("reference")[-1]
    recordType = O.getOutput("recordType")[0]
    descriptions = O.getOutput("NewDescriptions")
        #iName, jName, mType, cDescr, pDescr, gAcc, cAcc, pAcc, 
        #fullDescr, fullpDescr

    if len(descriptions) == 0:
        #No descriptions generated [unlikely]
        return
    if O.Summary()[0]:
        #There were errors during the run, return.
        return
    for descr in descriptions:
        if goi in descr[0] and toi in descr[1]: # Gene and Transcript
            if tDescr:
                # Already inserted a value in the tDescr
                tDescr, tList = [], descriptions
                break
            tDescr = descr

    tList = descriptions

    var = O.getOutput("variant")[-1]

    # Generate output
    outputline = ""
    if tDescr: #Filtering worked, only one Description left
        (gName, trName, mType, cDescr,
            pDescr, gAcc, cAcc, pAcc, fullD, fullpD) = tDescr

        gene = "%s_v%.3i" % (gName, int(trName))

        outputline += "%s\t%s\t%s\t" % (reference, gene, var)

        #Add genomic Description
        outputline += "%s\t" % O.getOutput("gDescription")[0]

        #Add coding Description & protein Description
        outputline += "%s\t%s\t" % (cDescr, pDescr)

        gc = cDescr and "%s:%s" % (gene, cDescr)
        gp = pDescr and "%s:%s" % (gene, pDescr)

        #Add mutation with GeneSymbols
        outputline += "%s\t%s\t" % (gc, gp)

        #Add References, should get genomic ref from parsed data
        if recordType == "LRG":
            gAcc = reference
        if recordType == "GB":
            geno = ["NC", "NG", "AC", "NT", "NW", "NZ", "NS"]
            for g in geno:
                if reference.startswith(g):
                    gAcc = reference
                    break
        outputline += "%s\t%s\t%s\t" % (gAcc or "", cAcc or "", pAcc or "")

    else:
        outputline += "\t"*11

    #Add list of affected transcripts "|" seperator
    if tList:
        outputline += "%s\t" % "|".join(e[-2] for e in tList)
        outputline += "%s\t" % "|".join(e[-1] for e in tList)
    else:
        outputline += "\t"*2

    #Link naar additional info:
    #outputline+="http://localhost/mutalyzer2/redirect?mutationName=%s" %\
    #        "todovariant"


    O.addOutput("batchDone", outputline)
#_createBatchOutput

def checkSubstitution(start_g, Arg1, Arg2, MUU, GenRecordInstance, O) :
    """
        Do a semantic check for substitutions, do the actual substitution
        and give it a name.

        Arguments:
            start_g           ; Genomic location of the substitution.
            Arg1              ; Nucleotide in the reference sequence.
            Arg2              ; Nucleotide in the mutated sequence.
            MUU               ; A Mutator object.
            GenRecordInstance ; A GenRecord object.
            O                 ; The Output object.
    """

    if not __checkDNA(Arg2) : # It must be DNA.
        #O.addMessage(__file__, 4, "ENODNA", "Invalid letter in input")
        return
    if Arg1 == Arg2 :         # And there must be a real change.
        O.addMessage(__file__, 3, "ENOVAR",
            "No mutation given (%c>%c) at position %i." % (
            Arg1, Arg1, start_g))

    MUU.subM(start_g, Arg2)
    GenRecordInstance.name(start_g, start_g, "subst", MUU.orig[start_g - 1],
                           Arg2, None)
#checkSubstitution

def checkDeletionDuplication(start_g, end_g, mutationType, MUU,
                             GenRecordInstance, O) :
    """
        Do a semantic check for a deletion or duplication, do the actual 
        deletion/duplication and give it a name.

        Arguments:
            start_g           ; Genomic start position of the del/dup.
            end_g             ; Genomic end position of the del/dup.
            mutationType      ; The type (del or dup).
            MUU               ; A Mutator object.
            GenRecordInstance ; A GenRecord object.
            O                 ; The Output object.
    """

    roll = __roll(MUU.orig, start_g, end_g)

    shift = roll[1]
    if GenRecordInstance.record.molType == 'n' :
        mRNA = GenRecordInstance.record.geneList[0].transcriptList[0
            ].mRNA.positionList
        for i in mRNA :
            if end_g <= i and end_g + roll[1] > i :
                print "ALARM"
                shift = i - end_g
                print shift
                break
            #if
        #for
    #if


    if shift : # FIXME, The warning may not be apropriate.
        newStart = start_g + shift
        newStop = end_g + shift
        O.addMessage(__file__, 2, "WROLL",
            "Sequence \"%s\" at position %s was given, however, " \
            "the HGVS notation prescribes that it should be \"%s\" at " \
            "position %s." % (
            MUU.visualiseLargeString(str(MUU.orig[start_g - 1:end_g])),
            __formatRange(start_g, end_g), 
            MUU.visualiseLargeString(str(MUU.orig[newStart - 1:newStop])),
            __formatRange(newStart, newStop)))
    #if
    if mutationType == "del" :
        MUU.delM(start_g, end_g)
    else :
        MUU.dupM(start_g, end_g)
    GenRecordInstance.name(start_g, end_g, mutationType, "", "",
                           (roll[0], shift))
#checkDeletionDuplication

def checkInversion(start_g, end_g, MUU, GenRecordInstance, O) :
    """
    """

    snoop = __palinsnoop(MUU.orig[start_g - 1:end_g])
    if snoop :
        if snoop == -1 :
            O.addMessage(__file__, 2, "WNOCHANGE",
                "Sequence \"%s\" at position %i_%i is a palindrome " \
                "(its own reverse complement)." % (
                MUU.visualiseLargeString(str(MUU.orig[start_g - 1:end_g])), 
                start_g, end_g))
            return
        #if
        else :
            O.addMessage(__file__, 2, "WNOTMINIMAL",
                "Sequence \"%s\" at position %i_%i is a partial " \
                "palindrome (the first %i nucleotide(s) are the reverse " \
                "complement of the last one(s)), the HGVS notation " \
                "prescribes that it should be \"%s\" at position %i_%i." % (
                MUU.visualiseLargeString(str(MUU.orig[start_g - 1:end_g])),
                start_g, end_g, snoop,
                MUU.visualiseLargeString(
                    str(MUU.orig[start_g + snoop - 1: end_g - snoop])),
                start_g + snoop, end_g - snoop))
            start_g += snoop
            end_g -= snoop
        #else
    #if
    MUU.invM(start_g, end_g)
    if start_g == end_g :
        O.addMessage(__file__, 2, "WWRONGTYPE", "Inversion at position "\
            "%i is actually a substitution." % start_g)
        GenRecordInstance.name(start_g, start_g, "subst", MUU.orig[start_g - 1],
            Bio.Seq.reverse_complement(MUU.orig[start_g - 1]), None)
    #if
    else :
        GenRecordInstance.name(start_g, end_g, "inv", "", "", None)
#checkInversion

def checkInsertion(start_g, end_g, Arg1, MUU, GenRecordInstance, O) :
    """
    """

    if start_g + 1 != end_g :
        O.addMessage(__file__, 3, "EINSRANGE",
            "%i and %i are not consecutive positions." % (start_g, end_g))
        return
    #if
    if not Arg1 or not __checkDNA(Arg1) :
        O.addMessage(__file__, 3, "EUNKVAR", "Although the syntax of this " \
            "variant is correct, the effect can not be analysed.")
        return
    #if

    MUU.insM(start_g, Arg1)
    insertionLength = len(Arg1)
    newStart = MUU.shiftpos(start_g)
    newStop = MUU.shiftpos(start_g) + insertionLength
    roll = __roll(MUU.mutated, newStart + 1, newStop)

    shift = roll[1]
    if GenRecordInstance.record.molType == 'n' :
        mRNA = GenRecordInstance.record.geneList[0].transcriptList[0
            ].mRNA.positionList
        for i in mRNA :
            if newStop <= i and newStop + roll[1] > i :
                print "ALARM"
                shift = i - newStop
                print shift
                break
            #if
        #for
    #if

    if roll[0] + shift >= insertionLength :
        O.addMessage(__file__, 2, "WINSDUP",
            "Insertion of %s at position %i_%i was given, " \
            "however, the HGVS notation prescribes that it should be a " \
            "duplication of %s at position %i_%i." % (
            MUU.mutated[newStart:newStop], start_g, start_g + 1,
            MUU.mutated[newStart:newStop], start_g + shift, 
            start_g + shift + insertionLength - 1))
        end_g += shift - 1
        start_g = end_g - insertionLength + 1
        GenRecordInstance.name(start_g, end_g, "dup", "", "",
                               (roll[0] + shift - insertionLength, 0))
    #if
    else :
        if shift :
            O.addMessage(__file__, 2, "WROLL", "Insertion of %s at position " \
                "%i_%i was given, however, the HGVS notation prescribes " \
                "that it should be an insertion of %s at position %i_%i." % (
                MUU.mutated[newStart:newStop], start_g, start_g + 1,
                MUU.mutated[newStart + shift:newStop + shift], 
                newStart + shift, newStart + shift + 1))
        GenRecordInstance.name(start_g, start_g + 1, "ins", 
            MUU.mutated[newStart + shift:newStop + shift] , "", 
            (roll[0], shift))
#checkInsertion

def __ivs2g(location, transcript) :
    """
    """

    ivsNumber = int(location.IVSNumber)

    if ivsNumber < 1 or ivsNumber > transcript.CM.numberOfIntrons() :
        return None

    if location.OffSgn == '+' :
        return transcript.CM.getSpliceSite(ivsNumber * 2 - 1) + \
            transcript.CM.orientation * int(location.Offset)
    return transcript.CM.getSpliceSite(ivsNumber * 2) - \
        transcript.CM.orientation * int(location.Offset)
#__ivs2g

def __ex2g(location, transcript) :
    """
    """

    numberOfExons = transcript.CM.numberOfExons()

    exNumberStart = int(location.EXNumberStart)
    if exNumberStart < 1 or exNumberStart > transcript.CM.numberOfExons() :
        return None
    start_g = transcript.CM.getSpliceSite(exNumberStart * 2 - 2)

    if location.EXNumberStop :
        exNumberStop = int(location.EXNumberStop)
        if exNumberStop < 1 or exNumberStop > transcript.CM.numberOfExons() :
            return None
        stop_g = transcript.CM.getSpliceSite(exNumberStop * 2 - 1)
    else :
        stop_g = transcript.CM.getSpliceSite(exNumberStart * 2 - 1)

    return start_g, stop_g
#__ex2g

def __normal2g(RawVar, transcript) :
    """
    """

    if not RawVar.StartLoc.PtLoc.Main.isdigit() :
        return None, None # For ? in a position.

    start_g = int(RawVar.StartLoc.PtLoc.Main)
    end_g = start_g
    if RawVar.EndLoc :
        if not RawVar.EndLoc.PtLoc.Main.isdigit() : # For ? in a position.
            return None, None
        end_g = transcript.CM.main2int(
            RawVar.EndLoc.PtLoc.MainSgn + RawVar.EndLoc.PtLoc.Main)
    #if


    # If it is not, convert it to g. notation.
    if transcript :
        start_main = transcript.CM.main2int(RawVar.StartLoc.PtLoc.MainSgn + \
                                            RawVar.StartLoc.PtLoc.Main)
        #if not RawVar.StartLoc.PtLoc.Offset.isdigit() :
        #    return

        start_offset = __PtLoc2offset(RawVar.StartLoc.PtLoc)

        if not __checkIntronPosition(start_main, start_offset, transcript) :
            return None, None

        start_g = transcript.CM.x2g(start_main, start_offset)
        end_g = start_g
        if RawVar.EndLoc :
            end_main = transcript.CM.main2int(RawVar.EndLoc.PtLoc.MainSgn + \
                                           RawVar.EndLoc.PtLoc.Main)
            #if not RawVar.EndLoc.PtLoc.Offset.isdigit() :
            #    return
            end_offset = __PtLoc2offset(RawVar.EndLoc.PtLoc)
            if not __checkIntronPosition(end_main, end_offset, transcript) :
                return None, None
            end_g = transcript.CM.x2g(end_main, end_offset)
        #if
        if transcript.CM.orientation == -1 :
            start_g, end_g = end_g, start_g
    #if

    return start_g, end_g
#__normal2g    

def __rv(MUU, RawVar, GenRecordInstance, parts, O, transcript) :
    """
    """

    # FIXME check this
    # First assume that the variant is given in g. notation.
    #print RawVar.StartLoc.PtLoc.MainSgn + RawVar.StartLoc.PtLoc.Main
    #print __PtLoc2offset(RawVar.StartLoc.PtLoc)

    Arg1 = RawVar.Arg1
    Arg2 = RawVar.Arg2
    
    if RawVar.EXLoc :
        start_g, end_g = __ex2g(RawVar.EXLoc, transcript)
        if not start_g :
            O.addMessage(__file__, 3, "EPOS", "Invalid EX position given.")
            return
        #if
        if end_g < start_g : # FIXME
            start_g, end_g = end_g, start_g
    #if
    else :
        if RawVar.StartLoc.IVSLoc :
            if GenRecordInstance.record.molType != 'g' :
                O.addMessage(__file__, 3, "ENOINTRON", "Intronic position " \
                    "given for a non-genomic reference sequence.")
                return
            start_g = __ivs2g(RawVar.StartLoc.IVSLoc, transcript)
            if not start_g :
                O.addMessage(__file__, 3, "EPOS", "Invalid IVS position given.")
                return
            #if
            end_g = start_g
            if RawVar.EndLoc and RawVar.EndLoc.IVSLoc : # FIXME
                end_g = __ivs2g(RawVar.EndLoc.IVSLoc, transcript)
                if end_g < start_g :
                    start_g, end_g = end_g, start_g
            #if
        #if
        else :        
            if GenRecordInstance.record.molType != 'g' and \
               (__intronicPosition(RawVar.StartLoc) or 
                __intronicPosition(RawVar.EndLoc)) :
                O.addMessage(__file__, 3, "ENOINTRON", "Intronic position " \
                    "given for a non-genomic reference sequence.")
                return
            start_g, end_g = __normal2g(RawVar, transcript)
            if not start_g :
                O.addMessage(__file__, 3, "ESPLICE", "Invalid intronic " \
                    "position given.")
                return
        #else
    #else            

    if start_g < 1 :
        O.addMessage(__file__, 4, "ERANGE", "Position %i is out of range." %
                     start_g)
        return
    #if
    if end_g > len(MUU.orig) :
        O.addMessage(__file__, 4, "ERANGE", "Position %s is out of range." %
                     end_g)
        return
    #if

    if transcript and transcript.CM.orientation == -1 :
        Arg1 = Bio.Seq.reverse_complement(RawVar.Arg1)
        Arg2 = Bio.Seq.reverse_complement(RawVar.Arg2)

    if RawVar.MutationType in ["del", "dup", "subst", "delins"] :
        __checkOptArg(MUU.orig, start_g, end_g, Arg1, O)

    if RawVar.MutationType == "subst" :
        checkSubstitution(start_g, Arg1, Arg2, MUU, GenRecordInstance, O)
    if RawVar.MutationType in ["del", "dup"] :
        checkDeletionDuplication(start_g, end_g, RawVar.MutationType, MUU,
                                 GenRecordInstance, O)
    if RawVar.MutationType == "inv" :
        checkInversion(start_g, end_g, MUU, GenRecordInstance, O)
    if RawVar.MutationType == "ins" :
        checkInsertion(start_g, end_g, Arg1, MUU, GenRecordInstance, O)


    # DelIns.
    if RawVar.MutationType == "delins" :
        if not Arg1 :
            Arg1 = MUU.orig[start_g - 1:end_g]

        if str(Arg1) == str(Arg2) :
            O.addMessage(__file__, 2, "WNOCHANGE",
                "Sequence \"%s\" at position %i_%i is identical to the " \
                "variant." % (
                MUU.visualiseLargeString(str(MUU.orig[start_g - 1:end_g])),
                start_g, end_g))
            return
        #if

        del_part, ins_part, lcp, lcs = __trim2(Arg1, Arg2)
        if not len(del_part) :
            O.addMessage(__file__, 2, "WWRONGTYPE", "The given DelIns " \
                         "is actually an insertion.")
            checkInsertion(start_g + lcp - 1, start_g + lcp, ins_part, MUU,
                           GenRecordInstance, O)
            return
        #if
        if len(del_part) == 1 and len(ins_part) == 1 :
            O.addMessage(__file__, 2, "WWRONGTYPE", "The given DelIns " \
                         "is actually a substitution.")
            checkSubstitution(start_g + lcp, del_part, ins_part, MUU,
                              GenRecordInstance, O)
            return
        #if
        if not len(ins_part) :
            O.addMessage(__file__, 2, "WWRONGTYPE", "The given DelIns " \
                         "is actually a deletion.")
            checkDeletionDuplication(start_g + lcp, end_g - lcs, "del",
                                     MUU, GenRecordInstance, O)
            return
        #if
        if str(Bio.Seq.reverse_complement(del_part)) == ins_part :
            O.addMessage(__file__, 2, "WWRONGTYPE", "The given DelIns " \
                         "is actually an inversion.")
            checkInversion(start_g + lcp, end_g - lcs, MUU,
                           GenRecordInstance, O)
            return
        #if

        MUU.delinsM(start_g + lcp, end_g - lcs, ins_part)

        GenRecordInstance.name(start_g, end_g, "delins", ins_part, "", None)
    #if
#__rv

def __ppp(MUU, parts, GenRecordInstance, O) :
    if parts.RawVar or parts.SingleAlleleVarSet :
        if parts.RefType == 'r' :
            O.addMessage(__file__, 4, "ERNA", "Descriptions on RNA level " \
                "are not supported.")
        if parts.RefType in ['c', 'n'] :
            GS, W = None, None
            goi, toi = O.getOutput("geneSymbol")[-1]
            if parts.LrgAcc:                   # LRG
                GS = GenRecordInstance.record.geneList[0] #LRG pick top gene
                if toi:
                    W = GS.findLocus(toi)
                    if not W:
                        O.addMessage(__file__, 4, "ENOTRANSCRIPT",
                            "Transcripts found for gene %s. Please " \
                            "choose from: %s" %(GS.name,
                                ", ".join(GS.listLoci())))
                else:                       # No transcript id given
                    if len(GS.transcriptList) == 1:
                        #No transcript given, only 1 found
                        W = GS.transcriptList[0]
                    else:
                        O.addMessage(__file__, 4, "ENOTRANSCRIPT",
                            "No transcript given for gene %s. Please " \
                            "choose from: %s" %(GS.name,
                                ", ".join(GS.listLoci())))

            #if
            else:
                # gene of interest
                genes = GenRecordInstance.record.listGenes()
                toi = toi and "%.3i" % int(toi)

                if goi in genes: #we found our gene
                    GS = GenRecordInstance.record.findGene(goi)
                elif (len(genes) == 1) and not(goi):
                    #There is only one gene in the Record, message?
                    GS = GenRecordInstance.record.geneList[0]
                else:
                    O.addMessage(__file__, 4, "EINVALIDGENE",
                        "Gene %s not found. Please choose from: %s" % (
                        goi, ", ".join(genes)))

                if GS:
                    #Find Transcript
                    transcripts = GS.listLoci()
                    if toi in transcripts:
                        W = GS.findLocus(toi)
                    elif (len(transcripts) == 1) and not(toi):
                        W = GS.transcriptList[0]
                    else:
                        O.addMessage(__file__, 4, "ENOTRANSCRIPT",
                            "Transcripts found for gene %s. Please " \
                            "choose from: %s" %(GS.name,
                            ", ".join(GS.listLoci())))
            #else

            # Add seletcted geneSymbol to output
            O.addOutput("geneSymbol", (GS and GS.name or "", W and W.name or ""))

            # Return if no transcript is selected
            if not W:
                #Skip all BatchJobs with the same preColon data
                O.addOutput("BatchFlags", ("S2",
                    O.getOutput("preColon")[-1]))
                return None #Explicit return in case of an error
        #if
        else :
            W = None
        #if W and not W.location :
        #    W = None
        if W and not W.transcribe :
            return

        if parts.SingleAlleleVarSet:
            for i in parts.SingleAlleleVarSet :
                __rv(MUU, i.RawVar, GenRecordInstance, parts, O, W)
        else :
            __rv(MUU, parts.RawVar, GenRecordInstance, parts, O, W)


        if not W : # Genomic given or error with transcript
            return
        if not GenRecordInstance.record.geneList : # EST
            return

        for i in range(0, W.CM.numberOfExons() * 2, 2) :
            exonStart = W.CM.getSpliceSite(i)
            exonStop = W.CM.getSpliceSite(i + 1)
            O.addOutput("exonInfo", [exonStart, exonStop,
                W.CM.g2c(exonStart), W.CM.g2c(exonStop)])

        O.addOutput("cdsStart_g", W.CM.x2g(1, 0))
        O.addOutput("cdsStart_c", 1)
        cdsStop = W.CM.info()[2]
        O.addOutput("cdsStop_g", W.CM.x2g(cdsStop, 0))
        O.addOutput("cdsStop_c", cdsStop)

        if W.transcribe :
            O.addOutput("myTranscriptDescription", W.description)

        if W.translate :
            cds = Seq(str(__splice(MUU.orig, W.CDS.positionList)),
                      IUPAC.unambiguous_dna)
            cdsm = Seq(str(__nsplice(MUU.mutated,
                                     MUU.newSplice(W.mRNA.positionList),
                                     MUU.newSplice(W.CDS.location),
                                     W.CM.orientation)),
                       IUPAC.unambiguous_dna)
            if W.CM.orientation == -1 :
                cds = Bio.Seq.reverse_complement(cds)
                cdsm = Bio.Seq.reverse_complement(cdsm)
            #if

            if '*' in cds.translate()[:-1] :
                O.addMessage(__file__, 3, "ESTOP", "In frame stop codon found.")
                return
            #if
            orig = cds.translate(table = W.txTable, to_stop = True)
            O.addOutput("oldprotein", orig + '*')
            trans = cdsm.translate(table = W.txTable, to_stop = True)

            if not trans or trans[0] != 'M' :
                __bprint(orig + '*', O, "oldProteinFancy")
                if str(cdsm[0:3]) in \
                    Bio.Data.CodonTable.unambiguous_dna_by_id[
                        W.txTable].start_codons :
                    O.addOutput("newprotein", '?')
                    __bprint('?', O, "newProteinFancy")
                    O.addOutput("altStart", str(cdsm[0:3]))
                    if str(orig[1:]) != str(trans[1:]) :
                        O.addOutput("altProtein", 'M' + trans[1:] + '*')
                        __bprint('M' + trans[1:] + '*', O, "altProteinFancy")
                #if
                else :
                    O.addOutput("newprotein", '?')
                    __bprint('?', O, "newProteinFancy")
                #else
            else :
                cdsLen = __cdsLen(MUU.newSplice(W.CDS.positionList))
                descr = __toProtDescr(cdsLen, orig, trans)
                print descr
                O.addOutput("myProteinDescription", descr[0])

                __bprint2(orig + '*', descr[1], descr[2], O, 
                    "oldProteinFancy")
                if str(orig) != str(trans) :
                    O.addOutput("newprotein", trans + '*')
                    __bprint2(trans + '*', descr[1], descr[3], O, 
                        "newProteinFancy")
            #else
        #if
    #if
#__ppp

def process(cmd, C, O) :
    parser = Parser.Nomenclatureparser(O)
    O.addOutput("inputvariant", cmd)
    ParseObj = parser.parse(cmd)
    del parser
    if not ParseObj :
        #Parsing went wrong
        return None     #Excplicit return of None in case of an error

    if ParseObj.Version :
        RetrieveRecord = ParseObj.RefSeqAcc + '.' + ParseObj.Version
    else :
        RetrieveRecord = ParseObj.RefSeqAcc

    D = Db.Cache(C.Db)
    if ParseObj.LrgAcc :
        filetype = "LRG"
        RetrieveRecord = ParseObj.LrgAcc
        geneSymbol = ("", ParseObj.LRGTranscriptID)
        retriever = Retriever.LargeRetriever(C.Retriever, O, D)
    else :
        if ParseObj.Gene:
            geneSymbol = (ParseObj.Gene.GeneSymbol or "",
                    ParseObj.Gene.TransVar or "")
            if ParseObj.Gene.ProtIso :
                O.addMessage(__file__, 4, "EPROT", "Indexing by protein " \
                    "isoform is not supported.")
        else:
            geneSymbol = ("", "")
        retriever = Retriever.GenBankRetriever(C.Retriever, O, D)
        filetype = "GB"

    # Store the recordType for output formatting
    O.addOutput("recordType", filetype)

    # Note concerning objects in outputObject, example:
    # O.getOutput('reference')[-1] countains the last added value
    # O.getOutput('reference')[0] countains the first added value
    # These can refer to the same element
    O.addOutput("reference", RetrieveRecord)

    # The geneSymbol[0] is used as a filter for batch runs
    O.addOutput("geneSymbol", geneSymbol) #tuple(Gene, TransV)

    # preColon is used to filter out Batch entries
    # that will result in identical errors
    O.addOutput("preColon", cmd.split(":")[0])
    O.addOutput("variant", cmd.split(":")[-1])

    record = retriever.loadrecord(RetrieveRecord)
    #if record and record.version and not '.' in RetrieveRecord : #FIXME
    #    O.addOutput("reference", RetrieveRecord + '.' + record.version)
    #else :

    if not record :
        return
    del retriever
    del D

    GenRecordInstance = GenRecord.GenRecord(O, C.GenRecord)
    GenRecordInstance.record = record
    GenRecordInstance.checkRecord()
    #NOTE:  GenRecordInstance is carrying the sequence in   .record.seq
    #       so is the Mutator.Mutator instance MUU          .orig

    MUU = Mutator.Mutator(GenRecordInstance.record.seq, C.Mutator, O)
    __ppp(MUU, ParseObj, GenRecordInstance, O)

    # PROTEIN
    for i in GenRecordInstance.record.geneList :
        #if i.location :
        for j in i.transcriptList :
            if not ';' in j.description and j.CDS and j.translate :
                cds = Seq(str(__splice(MUU.orig, j.CDS.positionList)), 
                          IUPAC.unambiguous_dna)
                cdsm = Seq(str(__nsplice(MUU.mutated, 
                                         MUU.newSplice(j.mRNA.positionList), 
                                         MUU.newSplice(j.CDS.location), 
                                         j.CM.orientation)),
                           IUPAC.unambiguous_dna)
                if j.CM.orientation == -1 :
                    cds = Bio.Seq.reverse_complement(cds)
                    cdsm = Bio.Seq.reverse_complement(cdsm)
                #if

                #if '*' in cds.translate()[:-1] :
                #    O.addMessage(__file__, 3, "ESTOP", 
                #                 "In frame stop codon found.")
                #    return
                ##if

                if not len(cds) % 3 :
                    orig = cds.translate(table = j.txTable, cds = True, 
                        to_stop = True)
                    trans = cdsm.translate(table = j.txTable, 
                        to_stop = True)

                    cdsLen = __cdsLen(MUU.newSplice(j.CDS.positionList))
                    j.proteinDescription = __toProtDescr(cdsLen, orig, 
                        trans)[0]
                #if
                else :
                    O.addMessage(__file__, 2, "ECDS", "CDS length is " \
                        "not a multiple of three in gene %s, transcript " \
                        "variant %s." % (i.name, j.name))
                    j.proteinDescription = '?'
    # /PROTEIN

    reference = O.getOutput("reference")[-1]
    if ';' in GenRecordInstance.record.description :
        descr = '['+GenRecordInstance.record.description+']'
    else:
        descr = GenRecordInstance.record.description

    O.addOutput("genomicDescription", "%s:%c.%s" % (reference,
        GenRecordInstance.record.molType, descr))
    O.addOutput("gDescription", "%c.%s" % (
        GenRecordInstance.record.molType, descr))
    O.addOutput("molType", GenRecordInstance.record.molType)

    if GenRecordInstance.record.chromOffset :
        if ';' in GenRecordInstance.record.chromDescription :
            chromDescr = '['+GenRecordInstance.record.chromDescription+']'
        else:
            chromDescr = GenRecordInstance.record.chromDescription

        O.addOutput("genomicChromDescription", "%s:%c.%s" % (
            GenRecordInstance.record.recordId, 
            GenRecordInstance.record.molType, chromDescr))
    #if

    if GenRecordInstance.record._sourcetype == "LRG": #LRG record
        for i in GenRecordInstance.record.geneList:
            for j in sorted(i.transcriptList, key = attrgetter("name")) :
                (iName, jName, mType, cDescr, pDescr,
                        gAcc, cAcc, pAcc, fullDescr, fullpDescr) =\
                    (i.name, j.name, j.molType, "", "", "", "", "", "", "")

                if ';' in j.description:
                    descr = '['+j.description+']'
                else:
                    descr = j.description

                if j.name:
                    fullDescr =\
                        "%st%s:%c.%s" % (reference, j.name, j.molType, descr)
                    O.addOutput("descriptions", fullDescr)
                #if
                else:
                    O.addOutput("descriptions", (i.name))

                if j.molType == 'c':
                    cDescr = "c.%s" % descr
                    pDescr = j.proteinDescription
                    fullpDescr = "%sp%s:%s" % (reference, j.name, pDescr)
                    O.addOutput("protDescriptions", fullpDescr)
                    cAcc, pAcc = j.transcriptID, j.proteinID
                #if

                O.addOutput("NewDescriptions", (
                    iName, jName, mType, cDescr, pDescr, gAcc,
                    cAcc, pAcc, fullDescr, fullpDescr))
            #for
        #for
    #if
    else :
        for i in GenRecordInstance.record.geneList :
            for j in sorted(i.transcriptList, key = attrgetter("name")) :
                (iName, jName, mType, cDescr, pDescr,
                        gAcc, cAcc, pAcc, fullDescr, fullpDescr) =\
                    (i.name, j.name, j.molType, "", "", "", "", "", "", "")

                if ';' in j.description :
                    descr = '['+j.description+']'
                else:
                    descr = j.description

                fullDescr = "%s(%s_v%s):%c.%s" % (reference,\
                        iName, jName, mType, descr)
                O.addOutput("descriptions", fullDescr)

                if (j.molType == 'c') :
                    cDescr = "c.%s" % descr
                    pDescr = j.proteinDescription
                    fullpDescr = "%s(%s_i%s):%s" % (
                        reference, iName, jName, pDescr)
                    O.addOutput("protDescriptions", fullpDescr)
                    cAcc, pAcc = j.transcriptID, j.proteinID
                #if

                O.addOutput("NewDescriptions", (
                    iName, jName, mType, cDescr, pDescr, gAcc,
                    cAcc, pAcc, fullDescr, fullpDescr))
            #for
        #for
    #else


    # LEGEND
    for i in GenRecordInstance.record.geneList :
        for j in sorted(i.transcriptList, key = attrgetter("name")) :

            if not j.name: continue #Exclude nameless transcripts

            O.addOutput("legends", ["%s_v%s" % (i.name, j.name),
                        j.transcriptID, j.locusTag, 
                        j.transcriptProduct, j.linkMethod])
            if j.translate :
                O.addOutput("legends", ["%s_i%s" % (i.name, j.name),
                    j.proteinID, j.locusTag, 
                    j.proteinProduct, j.linkMethod])
        #for

    #Add GeneSymbol and Transcript Var to the Output object for batch
    if ParseObj.Gene:
        O.addOutput("geneOfInterest", dict(ParseObj.Gene.items()))
    else:
        O.addOutput("geneOfInterest", dict())

    _createBatchOutput(O)

    del MUU

    return GenRecordInstance
    #if
#process

def main(cmd) :
    C = Config.Config()
    O = Output.Output(__file__, C.Output)

    O.addMessage(__file__, -1, "INFO", "Received variant " + cmd)

    RD = process(cmd, C, O)

    O.addMessage(__file__, -1, "INFO", "Finished processing variant " + cmd)

    ### OUTPUT BLOCK ###
    gn = O.getOutput("genename")
    if gn :
        print "Gene Name: " + gn[0]
    tv = O.getOutput("transcriptvariant")
    if tv :
        print "Transcript variant: " + tv[0]
        print
    #if

    for i in O.getMessages() :
        print i
    errors, warnings, summary = O.Summary()
    print summary
    print

    #if not errors :
    if not errors or DEBUG:
        visualisation = O.getOutput("visualisation")
        if visualisation :
            for i in range(len(visualisation)) :
                if i and not i % 3 :
                    print
                print visualisation[i]
            #for
            print
        #if

        reference = O.getOutput("reference")[-1]
        for i in O.getOutput("descriptions") :
            print i
        print
        for i in O.getOutput("protDescriptions") :
            print i
        print

        if RD.record and RD.record._sourcetype == "LRG": #LRG record
            from collections import defaultdict
            toutput = defaultdict(list)
            poutput = defaultdict(list)
            for i in RD.record.geneList:
                for j in i.transcriptList:
                    d = j.description
                    d = ';' in d and '['+d+']' or d
                    if j.name:
                        toutput[i.name].append(
                            "%st%s:%c.%s" % (reference, j.name, j.molType, d))
                    else:
                        pass
                    if j.molType == 'c':
                        poutput[i.name].append(
                                "%sp%s:%s" % (reference, j.name,
                                    j.proteinDescription))
                        poutput[i.name].sort()
                toutput[i.name].sort()

            #Transcript Notation
            print "Following transcripts were affected:"
            for key, values in toutput.items():
                print key
                for value in values:
                    print "\t"+value

            #Protein Notation
            print "\nFollowing proteins were affected:"
            for key, values in poutput.items():
                print key
                for value in values:
                    print "\t"+value
            #for
        #if
        else :
            for i in RD.record.geneList :
                for j in i.transcriptList :
                    if ';' in j.description :
                        print "%s(%s_v%s):%c.[%s]" % (reference, i.name, j.name,
                                                      j.molType, j.description)
                    else :
                        print "%s(%s_v%s):%c.%s" % (reference, i.name, j.name,
                                                    j.molType, j.description)
                        if (j.molType == 'c') :
                            print "%s(%s_i%s):%s" % (reference, i.name, j.name,
                                                     j.proteinDescription)
                    #else
                #for
            #for
        #else

        #Genomic Notation
        rdrd = RD.record.description
        gdescr = ';' in rdrd and '['+rdrd+']' or rdrd
        print "\nGenomic notation:\n\t%s:g.%s" % (reference, gdescr)
        print O.getOutput("genomicChromDescription")

        op = O.getOutput("oldprotein")
        if op :
            print "\nOld protein:"
            #__bprint(op[0], O)
            for i in O.getOutput("oldProteinFancy") :
                print i
            print
        #if
        np = O.getOutput("newprotein")
        if np :
            print "\nNew protein:"
            #__bprint(np[0], O)
            for i in O.getOutput("newProteinFancy") :
                print i
            print
        #if
        ap = O.getOutput("altProtein")
        if ap :
            print "\nAlternative protein using start codon %s:" % \
                O.getOutput("altstart")[0]
            #__bprint(ap[0], O)
            for i in O.getOutput("altProteinFancy") :
                print i
            print
        #if

        for i in O.getOutput("exonInfo") :
            print i
        print
        print O.getOutput("cdsStart")
        print O.getOutput("cdsStop")
        print

        for i in O.getOutput("legends") :
            print i

        print
        print "Restriction sites:"
        for i in O.getOutput("restrictionSites") :
            print i

        print "+++ %s" % O.getOutput("myTranscriptDescription")

    #if
    ### OUTPUT BLOCK ###
    del O
#main

if __name__ == "__main__" :
    if len(sys.argv) > 1:
        main(sys.argv[1])
#if
