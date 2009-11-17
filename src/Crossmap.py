#!/usr/bin/python

class Crossmap() :
    """
        Convert from g. to c. or n. notation or vice versa.
        
        Private variables:
            __STOP         ; CDS stop in c. notation.
            __crossmapping ; A list that contains either c. or n. positions
                             corresponding to the g. positions in the RNA list.

        Public variables:
            RNA         ; The list of RNA splice sites.
            CDS         ; The CDS list (if present).
            orientation ; The orientation of the transcript:  1 = forward
                                                             -1 = reverse.

        Special methods:
            __init__(RNA, CDS, orientation) ; Initialise the class and do the 
                                              cross mapping of the splice 
                                              sites.

        Private methods:
            __plus(a, b)   ; A protected '+' that skips 0 if 
                             a <= 0 and a + b >= 0.
            __minus(a, b)  ; A protected '-' that skips 0 if
                             a >= 0 and a - b <= 0.
            __minusr(a, b) ; A protected '-' that skips 0 if
                             a > 0 and b < 0.
            __star(a)      ; Translate from __STOP to '*' notation.
            __rstar(s)     ; Translate from '*' to __STOP notation.
            __crossmap_splice_sites() ; Calculate the __crossmapping list.
            g2x(a)    ; Translate from g. notation to c. or n. notation.
            x2g(a, b) ; Translate c. or n. notation to g. notation.
    """

    def __init__(self, RNA, CDS, orientation) :
        """
            Initialise the class and do the cross mapping of the splice sites.

            Arguments:
                RNA         ; The list of RNA splice sites.
                CDS         ; The CDS list (may be empty).
                orientation ; The orientation of the transcript.

            Private variables (altered):
                __STOP         ; CDS stop in c. notation.
                __crossmapping ; A list that contains either c. or n. positions
                                 corresponding to the g. positions in the RNA 
                                 list.

            Public variables (altered):
                RNA         ; The list of RNA splice sites.
                CDS         ; The CDS list (if present).
                orientation ; The orientation of the transcript:  1 = forward
                                                                 -1 = reverse.
        """
        self.__STOP = 10000000000
        self.__crossmapping = len(RNA) * [None]
        self.RNA = list(RNA)
        self.CDS = list(CDS)
        self.orientation = orientation

        self.__crossmap_splice_sites()
    #__init__

    def __plus(self, a, b) :
        """
            This method returns a + b unless a is smaller than zero and the
            result is larger than zero, in that case it returns a + b + 1. 
            In effect the number 0 is skipped while adding.
            
            Arguments: 
                a ; First argument of the addition.
                b ; Second argument of the addition.
            
            Returns:
                integer ; a + b or a + b + 1.

        """
        r = a + b

        if a <= 0 and r >= 0 :
            return r + 1

        return r
    #__plus

    def __minus(self, a, b) :
        """
            This method returns a - b unless a is larger than zero and the
            result is smaller than zero, in that case it returns (a - b) - 1. 
            In effect the number 0 is skipped while subtracting.
            
            Arguments: 
                a ; First argument of the subtraction.
                b ; Second argument of the subtraction.
            
            Returns:
                integer ; a - b or (a - b) - 1.
        """
        r = a - b

        if a >= 0 and r <= 0 :
            return r - 1

        return r
    #__minus

    def __minusr(self, a, b) :
        """
            This method returns a - b unless a is larger than zero and b is
            smaller than zero, in that case it returns (a - b) - 1. 
            In effect the number 0 is skipped while subtracting.
            
            Arguments: 
                a ; First argument of the subtraction.
                b ; Second argument of the subtraction.
            
            Returns:
                integer ; a - b or (a - b) - 1.
        """
        r = a - b

        if a > 0 and b < 0 :
            return r - 1

        return r
    #__minusr

    def __star(self, a) :
        """
            This method converts the __STOP notation to the '*' notation.

            Arguments:
                a ; An integer in __STOP notation.

            Returns:
                string ; The converted notation (may be unaltered).
        """
        if a > self.__STOP :
            return '*' + str(a - self.__STOP)

        return str(a)
    #__star

    def __rstar(self, s) :
        """
            This method converts the '*' notation to the __STOP notation.

            Arguments:
                s ; A string in '*' notation.

            Returns:
                integer ; The converted notation (may be unaltered).
        """
        if s[0] == '*' :
            return self.__STOP + int(s[1:])

        return int(s)
    #__rstar

    def __crossmap_splice_sites(self) :
        """
            This method calculates either:
              1: The c. notation of the CDS start and stop, including splice
                 sites.
              2: The c. notation of the RNA splice sites.
              3: The n. notation of the RNA splice sites.
            
            For option 1 only provide an list with CDS splice sites.
            For option 2 provide an list with RNA and one with CDS splice
              sites.
            For option 3 only provide an list with RNA splice sites.
            
            Examples:
            - Get the n. notation of the RNA splice sites. The input is in
              forward notation.
            Crossmap(RNA, [], 1)
            - Get the c. notation of the CDS start and stop, and the internal
              splice sites. The input is in forward notation.
            Crossmap(CDS, [], 1)
            - Get the c. notation of the RNA splice sites. The input is in
              reverse complement.
            Crossmap(RNA, CDS, -1)
            
            The output is straightforward, except for the c. notation of the
            downstream RNA splice sites. This is denoted by __STOP + the
            distance to the stop codon, as an alternative to the *-notation.

            Private variables (altered):
                __crossmapping ; A list that contains either c. or n. positions
                                 corresponding to the g. positions in the RNA 
                                 list.

            Private variables:
                __STOP         ; A large number to indicate positions after CDS
                                 stop.

            Public variables:
                RNA         ; The list of RNA splice sites.
                CDS         ; The CDS list (if present).
                orientation ; The orientation of the transcript:  1 = forward
                                                                 -1 = reverse.
        """
        CDSlen = len(self.CDS)
        RNAlen = len(self.RNA)
        cPos = 1 # This value stays one unless we both have mRNA and CDS.
        d = self.orientation
        c = (d - 1) / -2         # c, x and z are used to unify forward and
        x = c * (CDSlen - 3) + 1 # reverse complement
        y = c * (RNAlen - 1)

        if CDSlen : # If we both have mRNA and CDS, we have to search for start.
            i = y - c
            while d * self.RNA[i] < d * self.CDS[x] :  # Find CDS start.
                i += d
            cPos = d * (self.CDS[x] - self.CDS[x - d] + 
                        (d * 2)) # Get the right boundary.

            while d * i > d * y :                      # Go back to exon 1.
                cPos = self.__minus(cPos, d * 
                       (self.RNA[i] - self.RNA[i - d] + d))
                i -= d * 2
            #while
        #if

        i = y - c
        while d * i < d * (y - c) + RNAlen :
            self.__crossmapping[i + c] = cPos;
            if i % 2 :                      # We are skipping an intron.
                cPos = self.__plus(cPos, 1) # Only add 1 (mind the 0).
            else : # We are skipping an exon, so add the length.
                cPos = self.__plus(cPos, self.RNA[i + 1] - self.RNA[i])

            # Set __STOP when we find CDS stop.
            if CDSlen and cPos < self.__STOP and \
               d * self.RNA[i - c + 1] >= d * self.CDS[d * (CDSlen - 2) + x] :
                self.__STOP = cPos - (d * (self.RNA[i - c + 1] - 
                                      self.CDS[d * (CDSlen - 2) + x]))
            i += d
        #while
    #__crossmap_splice_sites

    def g2x(self, a) :
        """
            This function calculates either:
              1: The n. notation from a g. notation.
              2: The c. notation from a g. notation.
            
            For option 1 only provide an array with mRNA splice sites and one
            with the c. notation of the splice sites.
            For option 2 provide an array with mRNA splice sites, one with the
            c. notation of the splice sites and an array with the CDS splice
            sites.
            
            Examples:
            - Get the n. notation of a g. position i. The input is in forward
              notation. 
            Crossmap(RNA, [], 1)
            g2x(i)
            - Get the c. notation of a g. position i. The input is in reverse
              notation. 
            Crossmap(mRNA, CDS, -1);
            g2x(i);
            
            The output is fully compatible with the HVGS nomenclature as
            defined on 01-07-2009.

            Arguments:
                a ; The genomic position that must be translated.

            Private variables:
                __crossmapping ; A list that contains either c. or n. positions
                                 corresponding to the g. positions in the RNA 
                                 list.
                __STOP         ; A large number to indicate positions after CDS
                                 stop.

            Public variables:
                RNA         ; The list of RNA splice sites.
                CDS         ; The CDS list (if present).
                orientation ; The orientation of the transcript:  1 = forward
                                                                 -1 = reverse.

            Returns:
                string ; The c. or n. notation of position a.
        """
        CDSlen = len(self.CDS)
        RNAlen = len(self.RNA)
        d = self.orientation
        c = (d - 1) / -2     # c, x and z are used to unify forward and
        x = c * (CDSlen - 1) # reverse complement.
        y = c * (RNAlen - 1)

        if d * a < d * self.RNA[y] : # A position before the first exon.
            return str(self.__crossmapping[y]) + "-u" + \
                   str(d * (self.RNA[y] - a))
        if d * a > d * self.RNA[RNAlen - y - 1] : # After the last exon.
            return self.__star(self.__crossmapping[RNAlen - y - 1]) + "+d" + \
                               str(d * (a - self.RNA[RNAlen - y - 1]))

        for i in xrange(RNAlen) : # A "normal" position.
            if i % 2 :            # We're checking the intron positions.
                if self.RNA[i] < a and a < self.RNA[i + 1] : # Intron.
                    if d * (a - self.RNA[i]) > d * (self.RNA[i + 1] - a) :
                        # The position was closer to the next exon.
                        return self.__star(self.__crossmapping[i + 1 - c]) + \
                                           '-' + \
                                           str(d * (self.RNA[i + 1 - c] - a))
                    # The position was closer to the previous exon.
                    return self.__star(self.__crossmapping[i + c]) + '+' + \
                                       str(d * (a - self.RNA[i + c]))
                #if
            else :                # We're checking the exon positions.
                if self.RNA[i] <= a and a <= self.RNA[i + 1] :
                    if CDSlen and d * a > d * self.CDS[CDSlen - x - 1] :
                        # Use the next splice site after CDS STOP.
                        return self.__star(self.__minus( \
                                           self.__crossmapping[i + 1 - c], \
                                           d * (self.RNA[i + 1 - c] - a)))
                    # Use the previous splice site in other cases.
                    return self.__star(self.__plus(self.__crossmapping[i + c], \
                                                   d * (a - self.RNA[i + c])))
                #if
        #for
    #g2x

    def x2g(self, a, b) :
        """
            This function calculates either:
              1: The g. notation from a n. notation.
              2: The g. notation from a c. notation.
            
            Whether option 1 or 2 applies depends on the content of mRNAm.
            
            
            Examples:
            - Get the g. notation of a n. position i. The input is in forward
              notation.
            Crossmap(RNA, [], 1)
            x2g(i)
            - Get the g. notation of a c. position i with offset j. The input
              is in reverse notation.
            Crossmap(mRNA, CDS, -1);
            x2g(i, j);
            
            Arguments: 
                a ; The n. or c. position to be translated.
                b ; The offset of position a.

            Private variables:
                __crossmapping ; A list that contains either c. or n. positions
                                 corresponding to the g. positions in the RNA 
                                 list.

            Public variables:
                RNA         ; The list of RNA splice sites.
                orientation ; The orientation of the transcript:  1 = forward
                                                                 -1 = reverse.
            
            Returns:
                integer ; A g. position.
        """
        d = self.orientation
        c = (-d - 1) / -2 # Used to unify forward and reverse complement.
        RNAlen = len(self.RNA)

        # Assume a position before exon 1.
        ret = self.RNA[0] - d * (self.__crossmapping[0] - a)
        if d * a > d * self.__crossmapping[RNAlen - 1] :
            # It is after the last exon.
            ret = self.RNA[RNAlen - 1] + \
                  d * (a - self.__crossmapping[RNAlen - 1])
        for i in range(0, RNAlen, 2) : # Is it in an exon?
            if d * self.__crossmapping[i] <= d * a and \
               d * a <= d * self.__crossmapping[i + 1] :
                ret = self.RNA[i + c] - d * \
                      self.__minusr(self.__crossmapping[i + c], a)
        ret += d * b # Add the intron count.

        return ret
    #x2g

    def c2str(self, mainsgn, main, offsgn, offset) :
        return str(mainsgn) + str(main) + str(offsgn) + str(offset)
    #printc

    def star2abs(self, mainsgn, main) :
        return self.__rstar(mainsgn + str(main))
    #star2abs

    def off2int(self, offsgn, offset) :
        if not offset :
            return 0
        ret = int(offset)
        if offsgn == '-' :
            return -ret
        return ret
    #off2int

    def c2g(self, mainsgn, main, offsgn, offset) :
        m = mainsgn + main
        o = 0
        if offset :
            o = int(offset)
            if offsgn == '-' :
                o = -o
    
        return self.x2g(self.__rstar(m), o)
    #c2g
#Crossmap

#
# Unit test.
#
if __name__ == "__main__" :
    # Build a crossmapper for a hypothetical gene.
    mRNAf = [5002, 5125, 27745, 27939, 58661, 58762, 74680, 74767, 103409, 
             103528, 119465, 119537, 144687, 144810, 148418, 149215]
    CDSf = [27925, 27939, 58661, 58762, 74680, 74736]
    Cf = Crossmap(mRNAf, CDSf, 1)

    # Build a crossmapper for a hypothetical gene on the reverse strand.
    mRNAr = [2000, 2797, 6405, 6528, 31678, 31750, 47687, 47806, 76448, 76535, 
             92453, 92554, 123276, 123470, 146090, 146213]
    CDSr = [76479, 76535, 92453, 92554, 123276, 123290]
    Cr = Crossmap(mRNAr, CDSr, -1)

    # Check whether the gene on the forward strand has the right splice sites
    #   in c. notation.
    if Cf._Crossmap__crossmapping != [-304, -181, -180, 15, 16, 117, 118, 205, 
                                      206, 325, 326, 398, 399, 522, 523, 1320] :
        print "Crossmapping list not built correctly (c. notation)."

    # Check whether the splice sites in c. notation are the same as the ones on
    #   the reverse strand.
    if Cf._Crossmap__crossmapping != Cr._Crossmap__crossmapping[::-1] :
        print "Forward/reverse discrepancy (c. notation)."

    # Do some g. to c. conversion checking for the gene on the forward strand.
    if Cf.g2x(5001) != "-304-u1" or\
       Cf.g2x(5124) != "-182" or \
       Cf.g2x(5126) != "-181+1" or \
       Cf.g2x(27924) != "-1" or \
       Cf.g2x(27925) != "1" or \
       Cf.g2x(58660) != "16-1" or \
       Cf.g2x(74736) != "174" or \
       Cf.g2x(74737) != "*1" or \
       Cf.g2x(103408) != "*32-1" or \
       Cf.g2x(103410) != "*33" or \
       Cf.g2x(149216) != "*1146+d1" :
        print "Forward g. to c. conversion incorrect."

    # Do some g. to c. conversion checking for the gene on the reverse strand.
    if Cr.g2x(146214) != "-304-u1" or \
       Cr.g2x(146091) != "-182" or \
       Cr.g2x(146089) != "-181+1" or \
       Cr.g2x(123291) != "-1" or \
       Cr.g2x(123290) != "1" or \
       Cr.g2x(92555) != "16-1" or \
       Cr.g2x(76479) != "174" or \
       Cr.g2x(76478) != "*1" or \
       Cr.g2x(47807) != "*32-1" or \
       Cr.g2x(47805) != "*33" or \
       Cr.g2x(1999) != "*1146+d1" :
        print "Reverse g. to c. conversion incorrect."

    # Do some c. to g. conversion checking for the gene on the forward strand.
    if Cf.x2g(-304, -1) != 5001 or \
       Cf.x2g(-182, 0) != 5124 or \
       Cf.x2g(-181, 1) !=  5126 or \
       Cf.x2g(-1, 0) != 27924 or \
       Cf.x2g(1, 0) != 27925 or \
       Cf.x2g(16, -1) != 58660 or \
       Cf.x2g(174, 0) != 74736 or \
       Cf.x2g(Cf._Crossmap__rstar("*1"), 0) != 74737 or \
       Cf.x2g(Cf._Crossmap__rstar("*32"), -1) != 103408 or \
       Cf.x2g(Cf._Crossmap__rstar("*33"), 0) != 103410 or \
       Cf.x2g(Cf._Crossmap__rstar("*1146"), 1) != 149216 :
        print "Forward c. to g. conversion incorrect."

    # Do some c. to g. conversion checking for the gene on the reverse strand.
    if Cr.x2g(-304, -1) != 146214 or \
       Cr.x2g(-182, 0) != 146091 or \
       Cr.x2g(-181, 1) !=  146089 or \
       Cr.x2g(-1, 0) != 123291 or \
       Cr.x2g(1, 0) != 123290 or \
       Cr.x2g(16, -1) != 92555 or \
       Cr.x2g(174, 0) != 76479 or \
       Cr.x2g(Cf._Crossmap__rstar("*1"), 0) != 76478 or \
       Cr.x2g(Cf._Crossmap__rstar("*32"), -1) != 47807 or \
       Cr.x2g(Cf._Crossmap__rstar("*33"), 0) != 47805 or \
       Cr.x2g(Cf._Crossmap__rstar("*1146"), 1) != 1999 :
        print "Reverse c. to g. conversion incorrect."

    # Build a crossmapper for the hypothetical gene, missing the first exon and
    #   the last two exons.
    mRNAf2 = [27745, 27939, 58661, 58762, 74680, 74767, 103409, 
             103528, 119465, 119537]
    CDSf2 = [27925, 27939, 58661, 58762, 74680, 74736]
    Cf2 = Crossmap(mRNAf2, CDSf2, 1)

    # Build a crossmapper for the hypothetical gene on the reverse strand, 
    #   missing the first exon and the last two exons.
    mRNAr2 = [31678, 31750, 47687, 47806, 76448, 76535, 
              92453, 92554, 123276, 123470]
    CDSr2 = [76479, 76535, 92453, 92554, 123276, 123290]
    Cr2 = Crossmap(mRNAr2, CDSr2, -1)

    if Cf.g2x(27925) != Cf2.g2x(27925) or \
       Cr.g2x(123290) != Cr2.g2x(123290) :
        print "Failed missing UTR exon test."

    # Build a crossmapper for the hypothetical gene, but now non-coding.
    Cf3 = Crossmap(mRNAf, [], 1)

    # Build a crossmapper for the hypothetical gene on the reverse strand, but 
    #   now non-coding.
    Cr3 = Crossmap(mRNAr, [], -1)

    # Check whether the gene on the forward strand has the right splice sites
    #   in n. notation.
    if Cf3._Crossmap__crossmapping != [1, 124, 125, 319, 320, 421, 422, 509, 
                                       510, 629, 630, 702, 703, 826, 827, 
                                       1624] :
        print "Crossmapping list not built correctly (n. notation)."

    # Check whether the splice sites in n. notation are the same as the ones on
    #   the reverse strand.
    if Cf3._Crossmap__crossmapping != Cr3._Crossmap__crossmapping[::-1] :
        print "Forward/reverse discrepancy (n. notation)."

    # Do some g. to n. conversion checking for the gene on the forward strand.
    if Cf3.g2x(5001) != "1-u1" or \
       Cf3.g2x(5002) != "1" or \
       Cf3.g2x(5126) != "124+1" or \
       Cf3.g2x(149216) != "1624+d1" :
        print "Forward g. to n. conversion incorrect."

    # Do some g. to n. conversion checking for the gene on the reverse strand.
    if Cr3.g2x(146214) != "1-u1" or \
       Cr3.g2x(146213) != "1" or \
       Cr3.g2x(146089) != "124+1" or \
       Cr3.g2x(1999) != "1624+d1" :
        print "Forward g. to n. conversion incorrect."

    # Do some n. to g. conversion checking for the gene on the forward strand.
    if Cf3.x2g(1, -1) != 5001 or \
       Cf3.x2g(1, 0) != 5002 or \
       Cf3.x2g(124, 1) !=  5126 or \
       Cf3.x2g(1624, 1) != 149216 :
        print "Forward n. to g. conversion incorrect."

    # Do some n. to g. conversion checking for the gene on the reverse strand.
    if Cr3.x2g(1, -1) != 146214 or \
       Cr3.x2g(1, 0) != 146213 or \
       Cr3.x2g(124, 1) !=  146089 or \
       Cr3.x2g(1624, 1) != 1999 :
        print "Reverse n. to g. conversion incorrect."
#if
