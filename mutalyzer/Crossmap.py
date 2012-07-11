#!/usr/bin/python

"""
Module for conversion from genomic coordinates to coding sequence
orientated coordinates and vice versa.
The conversions are done based upon a list of splice sites, the CDS start
and stop and the orientation of a transcript.

"""
#Public classes:
#    - Crossmap ; Convert from g. to c. or n. notation or vice versa.

class Crossmap() :
    """
    Convert from I{g.} to I{c.} or I{n.} notation or vice versa.

    Private variables:
        - __STOP         ; CDS stop in I{c.} notation.
        - __crossmapping ; A list that contains either I{c.} or I{n.} positions
                           corresponding to the I{g.} positions in the RNA list.

    Public variables:
        - RNA         ; The list of RNA splice sites.
        - CDS         ; CDS start and stop (if present).
        - orientation ; The orientation of the transcript:  1 = forward
                                                           -1 = reverse.

    Special methods:
        - __init__(RNA, CDS, orientation) ; Initialise the class and do the
                                            cross mapping of the splice sites.

    Private methods:
        - __plus(a, b)              ; A protected '+' that skips 0 if
                                      a <= 0 and a + b >= 0.
        - __minus(a, b)             ; A protected '-' that skips 0 if
                                      a >= 0 and a - b <= 0.
        - __minusr(a, b)            ; A protected '-' that skips 0 if
                                      a > 0 and b < 0.
        - __crossmap_splice_sites() ; Calculate the __crossmapping list.

    Public methods:
        - int2main(a) ; Translate from __STOP to '*' notation.
        - main2int(s) ; Translate from '*' to __STOP notation.
        - g2x(a) ; Translate from I{g.} notation to I{c.} or I{n.} notation.
        - x2g(a, b) ; Translate I{c.} or I{n.} notation to I{g.} notation.
        - int2offset(t) ; Convert a tuple of integers to offset-notation.
        - offset2int(s) ; Convert an offset in HGVS notation to an integer.
        - tuple2string(t) ; Convert a tuple (main, offset) in __STOP notation
            to I{c.} notation.
        - g2c(a) ; Uses both g2x() and tuple2string() to translate a genomic
            position to __STOP notation to I{c.} notation.
        - info() ; Return transcription start, transcription end and CDS stop.
        - getSpliceSite(number) ; Return the coordinate of a splice site.
        - numberOfIntrons() ; Returns the number of introns.
        - numberOfExons() ; Returns the number of exons.
    """

    def __init__(self, RNA, CDS, orientation) :
        """
        Initialise the class and do the cross mapping of the splice sites.

        Private variables (altered):
            - __STOP         ; CDS stop in I{c.} notation.
            - __crossmapping ; A list that contains either I{c.} or I{n.}
                               positions corresponding to the I{g.} positions in
                               the RNA list.
            - __trans_start  ; Transcription start site in I{c.} notation.
            - __trans_end    ; Transcription end side in I{c.} notation.

        @arg RNA: The list of RNA splice sites
        @type RNA: list
        @arg CDS: CDS start and stop (if present, may be empty)
        @type CDS: list
        @arg orientation: The orientation of the transcript
            - 1 = forward
            - E{-}1 = reverse
        @type orientation: integer
        """
#        Public variables (altered):
#            - RNA         ; The list of RNA splice sites.
#            - CDS         ; CDS start and stop (if present).
#            - orientation ; The orientation of the transcript:  1 = forward
#                                                             -1 = reverse.

        self.__STOP = None
        self.__crossmapping = len(RNA) * [None]
        self.RNA = list(RNA)
        self.CDS = list(CDS)
        self.orientation = orientation

        self.__crossmap_splice_sites()

        start = (orientation - 1) / 2
        self.__trans_start = self.__crossmapping[start]
        self.__trans_end = self.__crossmapping[start - self.orientation]

        if not self.__STOP :
            self.__STOP = self.__trans_end
    #__init__

    def __plus(self, a, b) :
        """
        This method returns a + b unless a is smaller than zero and the
        result is larger than zero, in that case it returns a + b + 1.
        In effect the number 0 is skipped while adding.

        @arg a: First argument of the addition
        @type a: integer
        @arg b: Second argument of the addition
        @type b: integer

        @return: a + b or a + b + 1
        @rtype: integer
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

        @arg a: First argument of the subtraction
        @type a: integer
        @arg b: Second argument of the subtraction
        @type b: integer

        @return: a - b or (a - b) - 1
        @rtype: integer
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

        @arg a: First argument of the subtraction
        @type a: integer
        @arg b: Second argument of the subtraction
        @type b: integer

        @return: a - b or (a - b) - 1
        @rtype: integer
        """

        r = a - b

        if a > 0 and b < 0 :
            return r - 1

        return r
    #__minusr

    def __crossmap_splice_sites(self) :
        """
        This method calculates either:
            1. The I{c.} notation of the CDS start and stop, including splice
               sites.
            2. The I{c.} notation of the RNA splice sites.
            3. The I{n.} notation of the RNA splice sites.

        For option 1 only provide an list with CDS splice sites.
        For option 2 provide an list with RNA splice sites and one with
        the CDS start and stop.
        For option 3 only provide an list with RNA splice sites.

        Examples:

        Crossmap(RNA, [], 1)
            - Get the I{n.} notation of the RNA splice sites. The input is in
              forward notation.

        Crossmap(CDS, [], 1)
            - Get the I{c.} notation of the CDS start and stop, and the internal
              splice sites. The input is in forward notation.

        Crossmap(RNA, CDS, -1)
            - Get the I{c.} notation of the RNA splice sites. The input is in
              reverse complement.


        The output is straightforward, except for the I{c.} notation of the
        downstream RNA splice sites. This is denoted by __STOP + the
        distance to the stop codon, as an alternative to the *-notation.

        Private variables (altered):
            - __crossmapping ; A list that contains either I{c.} or I{n.} positions
                               corresponding to the I{g.} positions in the RNA
                               list.

        Private variables:
            - __STOP         ; A large number to indicate positions after CDS
                               stop.

        Public variables:
            - RNA         ; The list of RNA splice sites.
            - CDS         ; CDS start and stop (if present).
            - orientation ; The orientation of the transcript:  1 = forward,
                                                               -1 = reverse.
        """

        RNAlen = len(self.RNA)
        cPos = 1 # This value stays one unless we both have mRNA and CDS.
        d = self.orientation
        c = (d - 1) / -2   # c, x and y are used to unify forward and
        x = (-d - 1) / -2  # reverse complement.
        y = c * (RNAlen - 1)

        if self.CDS : # If we both have mRNA and CDS, we have to search for
                      # CDS start.
            i = y - c
            # Find CDS start.
            while d * (self.RNA[i] - ((i + 1) % 2)) < d * self.CDS[c] + c :
                i += d
            cPos = d * (self.RNA[i] - self.CDS[c] + (d * 2)) # Get the right
                                                             # boundary.

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
            if self.CDS and not self.__STOP and \
               d * self.RNA[i - c + 1] >= d * self.CDS[x] :
                self.__STOP = cPos - (d * (self.RNA[i - c + 1] - self.CDS[x]))
            i += d
        #while
    #__crossmap_splice_sites

    def g2x(self, a) :
        """
        This function calculates either:
            1. The I{n.} notation from a I{g.} notation.
            2. The I{c.} notation from a I{g.} notation.

        For option 1 only provide an array with mRNA splice sites and one
        with the I{c.} notation of the splice sites.
        For option 2 provide an array with mRNA splice sites, one with the
        I{c.} notation of the splice sites and an array with the CDS start and
        stop.

        Examples:

        Crossmap(RNA, [], 1)
        g2x(i)
            - Get the I{n.} notation of a I{g.} position i. The input is in forward
            notation.

        Crossmap(mRNA, CDS, -1);
        g2x(i);
            - Get the I{c.} notation of a I{g.} position i. The input is in reverse
              notation.

        The output is fully compatible with the HVGS nomenclature as defined
        on 01-07-2009.

        Private variables:
            - __crossmapping ; A list that contains either I{c.} or I{n.} positions
                               corresponding to the I{g.} positions in the RNA
                               list.
            - __STOP         ; A large number to indicate positions after CDS
                               stop.

        Public variables:
            - RNA         ; The list of RNA splice sites.
            - CDS         ; CDS start and stop (if present).
            - orientation ; The orientation of the transcript:  1 = forward,
                                                               -1 = reverse.

        @arg a: The genomic position that must be translated
        @type a: integer

        @return: The I{c.} or I{n.} notation of position a
        @rtype: string
        """
        # TODO update documentation.

        RNAlen = len(self.RNA)
        d = self.orientation
        c = (d - 1) / -2     # c and y are used to unify forward and reverse
        y = c * (RNAlen - 1) # complement.

        if d * a < d * self.RNA[y] : # A position before the first exon.
            return ((self.__crossmapping[y]), -d * (self.RNA[y] - a))
        if d * a > d * self.RNA[RNAlen - y - 1] : # After the last exon.
            return (self.__crossmapping[RNAlen - y - 1],
                    d * (a - self.RNA[RNAlen - y - 1]))

        for i in xrange(RNAlen) : # A "normal" position.
            if i % 2 :            # We're checking the intron positions.
                if self.RNA[i] < a and a < self.RNA[i + 1] : # Intron.
                    if d * (a - self.RNA[i]) > d * (self.RNA[i + 1] - a) :
                        # The position was closer to the next exon.
                        return (self.__crossmapping[i + 1 - c],
                                -d * (self.RNA[i + 1 - c] - a))
                    # The position was closer to the previous exon.
                    return (self.__crossmapping[i + c],
                            d * (a - self.RNA[i + c]))
                #if
            else :                # We're checking the exon positions.
                if self.RNA[i] <= a and a <= self.RNA[i + 1] :
                    return (self.__plus(self.__crossmapping[i + c],
                                        d * (a - self.RNA[i + c])), 0)
                #if
        #for
    #g2x

    def x2g(self, a, b) :
        """
        This function calculates either:
            1. The I{g.} notation from a I{n.} notation.
            2. The I{g.} notation from a I{c.} notation.

        Whether option 1 or 2 applies depends on the content of mRNAm.


        Examples:

        Crossmap(RNA, [], 1)
        x2g(i)
            - Get the I{g.} notation of a I{n.} position i. The input is in forward
            notation.

        Crossmap(mRNA, CDS, -1);
        x2g(i, j);
            - Get the I{g.} notation of a I{c.} position i with offset j. The input
            is in reverse notation.

        Private variables:
            - __crossmapping ; A list that contains either I{c.} or I{n.} positions
                             corresponding to the I{g.} positions in the RNA
                             list.

        Public variables:
            - RNA         ; The list of RNA splice sites.
            - orientation ; The orientation of the transcript:  1 = forward
                                                             -1 = reverse.

        @arg a: The I{n.} or I{c.} position to be translated
        @type a: integer
        @arg b: The offset of position a
        @type b: integer

        @return: A I{g.} position
        @rtype: integer
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

        if a < 0 and self.__crossmapping[d - c] == 1 : # Patch for CDS start on
            ret += d                               # first nucleotide of exon 1.

        return ret
    #x2g

    def int2main(self, a) :
        """
        This method converts the __STOP notation to the '*' notation.

        Private variables:
            - __STOP ; CDS stop in I{c.} notation.

        @arg a: An integer in __STOP notation
        @type a: integer

        @return: The converted notation (may be unaltered)
        @rtype: string
        """

        if a > self.__STOP :
            return '*' + str(a - self.__STOP)

        return str(a)
    #int2main

    def main2int(self, s) :
        """
        This method converts the '*' notation to the __STOP notation.

        Private variables:
            - __STOP ; CDS stop in I{c.} notation.

        @arg s: A string in '*' notation
        @type s: string

        @return: The converted notation (may be unaltered)
        @rtype: integer
        """
        if s[0] == '*':
            return self.__STOP + int(s[1:])

        return int(s)
    #main2int

    def int2offset(self, t, fuzzy=False):
        """
        Convert a tuple of integers to offset-notation. This adds a `+',
        and `u' or `d' to the offset when appropriate. The main value is
        not returned.

        @arg t: A tuple of integers: (main, offset) in __STOP notation
        @type t: tuple
        @kwarg fuzzy: Denotes that the coordinate is fuzzy (i.e. offset is
            unknown).
        @type fuzzy: bool

        @return: The offset in HGVS notation
        @rtype: string
        """

        if t[1] > 0 :                      # The exon boundary is downstream.
            if fuzzy: return '+?'
            if t[0] >= self.__trans_end :  # It is downstream of the last exon.
                return "+d" + str(t[1])
            return '+' + str(t[1])
        #if
        if t[1] < 0 :                       # The exon boundary is uptream.
            if fuzzy: return '-?'
            if t[0] <= self.__trans_start : # It is upstream of the first exon.
                return "-u" + str(-t[1])
            return str(t[1])
        #if
        return ''                           # No offset was given.
    #int2offset

    def offset2int(self, s) :
        """
        Convert an offset in HGVS notation to an integer. This removes
        `+', `u' and `d' when present. It also converts a `?' to something
        sensible.

        @arg s: An offset in HGVS notation
        @type s: string

        @return: The offset as an integer
        @rtype: integer
        """

        if not s :    # No offset given.
            return 0
        if s == '?' : # Here we ignore an uncertainty.
            return 0  # FIXME, this may have to be different.
        if s[1] == 'u' or s[1] == 'd' : # Remove `u' or `d'.
            if s[0] == '-' :            # But save the `-'.
                return -int(s[2:])
            return int(s[2:])
        #if
        if s[1:] == '?' : # Here we ignore an unvertainty in the intron.
            return 0      # FIXME, this may have to be different.
        if s[0] == '-' :
            return -int(s[1:])          # Save the `-' here too.
        return int(s[1:])
    #offset2int

    def tuple2string(self, t, fuzzy=False) :
        """
        Convert a tuple (main, offset) in __STOP notation to I{c.} notation.

        @arg t: A tuple (main, offset) in __STOP notation
        @type t: tuple
        @kwarg fuzzy: Denotes that the coordinate is fuzzy (i.e. offset is
            unknown).
        @type fuzzy: bool

        @return: The position in HGVS notation
        @rtype: string
        """

        if t[0] >= self.__trans_end or t[0] <= self.__trans_start:
            return str(self.int2main(self.__minus(t[0], -t[1])))
        return str(self.int2main(t[0])) + str(self.int2offset(t, fuzzy))
    #tuple2string

    def g2c(self, a, fuzzy=False) :
        """
        Uses both g2x() and tuple2string() to translate a genomic position
        to __STOP notation to I{c.} notation.

        @arg a: The genomic position that must be translated
        @type a: integer
        @kwarg fuzzy: Denotes that the coordinate is fuzzy (i.e. offset is
            unknown).
        @type fuzzy: bool

        @return: The position in HGVS notation
        @rtype: string
        """
        return self.tuple2string(self.g2x(a), fuzzy)
    #g2c

    def info(self) :
        """
        Return transcription start, transcription end and CDS stop.

        @return: (trans_start, trans_stop, CDS_stop)
        @rtype: triple
        """

        return (self.__trans_start, self.__trans_end, self.__STOP)
    #info

    def getSpliceSite(self, number) :
        """
        Return the coordinate of a splice site.

        @arg number: the number of the RNA splice site counting from
            transcription start.
        @type number: integer
        @return: coordinate of the RNA splice site.
        @rtype: integer
        """

        if self.orientation == 1 :
            return int(self.RNA[number])
        return int(self.RNA[len(self.RNA) - number - 1])
    #getSpliceSite

    def numberOfIntrons(self) :
        """
        Returns the number of introns.

        @return: number of introns
        @rtype: integer
        """

        return len(self.RNA) / 2 - 1
    #numberOfIntrons

    def numberOfExons(self) :
        """
        Returns the number of exons.

        @return: number of exons
        @rtype: integer
        """

        return len(self.RNA) / 2
    #numberOfExons
#Crossmap
