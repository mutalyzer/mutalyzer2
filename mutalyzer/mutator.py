"""
Module for mutating a string.

Mutations are described in the original coordinates. These coordinates are
transfered to the mutated coordinates with the aid of an internal shift
list, which keeps track of the sizes of changes. Using the original
coordinates greatly simplifies combined mutations in a variant. A
visualisation of each raw variant within a combined variant is made and
effects on restriction sites are also analysed.

The original as well as the mutated string are stored here.
"""


from Bio import Restriction
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import reverse_complement

from mutalyzer import util
from mutalyzer import config


class Mutator() :
    """
    Mutate a string and register all shift points. For each mutation a
    visualisation is made (on genomic level) and the addition or deletion
    of restriction sites is detected. Output for each raw variant is stored
    in the output object as 'visualisation', 'deletedRestrictionSites' and
    'addedRestrictionSites' respectively.

    Private variables:
        - __output           ; The output object.
        - __shift            ; A sorted list of tuples (position, shiftsize)
                               where the modifications in length are stored.
                               Each first element of the tuples in this list
                               is unique, each second element is non-zero.
        - __removed_sites    ; Set of splice sites to ignore in mutated
                               string.
        - __restrictionBatch ;

    Public variables:
        - orig    ; The original string.
        - mutated ; The mutated string.

    Special methods:
        - __init__(orig) ; Initialise the class with the original string.

    Private methods:
        - __sortins(tuple)      ; Insert a tuple in a sorted list, after
                                  insertion the list stays sorted.
        - __makeRestrictionSet()
        - __mutate(pos1, pos2, ins) ; A general mutation function that does a
                                      delins on interbase coordinates of the
                                      original string.

    Public methods:
        - shiftpos(position)       ; Calculate the position in the mutated
                                     string given the position in the
                                     original string.
        - newSplice(sites)         ; Generate a list of new splice sites.
        - delM(pos1, pos2)         ; Delete a range from non-interbase
                                     position pos1 to pos2.
        - insM(pos, ins)           ; Insert a string at interbase position
                                     pos.
        - delimsM(pos1, pos2, ins) ; Delete a range from non-interbase
                                     position pos1 to pos2 and insert ins.
        - subM(pos, nuc)           ; Substitute a nucleotite at non-interbase
                                     position pos for nuc.
        - invM(pos1, pos2)         ; Invert a range from non-interbase
                                     position pos1 to pos2.
        - dupM(pos1, pos2)         ; Duplicate a range from non-interbase
                                     position pos1 to pos2.
    """

    def __init__(self, orig, output) :
        """
        Initialise the class with the original string.

        Private variables (altered):
            - __output           ; Initialised with the output object.
            - __shift            ; Initialised to the empty list.
            - __restrictionBatch ; Initialised to a default set of
                                   restriction enzymes.

        Public variables (altered):
            - orig    ; Initialised to the parameter orig.
            - mutated ; Initialised to the parameter orig.

        @arg orig:   The original string before mutation
        @type orig: string
        @arg output: The output object
        @type output: object
        """
        self.__output = output
        self.__shift = []
        self.__removed_sites = set()
        self.__restrictionBatch = Restriction.RestrictionBatch([], ['N'])

        self.orig = orig
        self.mutated = orig
    #__init__

    def __sortins(self, tuple) :
        """
        Insert a tuple in a sorted list, the list is sorted on the first
        element of the tuples. After insertion the list stays sorted.
        If a tuple is inserted where tuple[0] already exists, this entry
        is altered.
        If an altered entry has zero as its second element, the entry is
        removed.

        Private variables (altered):
            - __shift ; A tuple can be added, removed or altered.

        @arg tuple: An ordered pair where tuple[0] denotes a position and
                    tuple[1] denotes the change in shift at this position
        @type tuple: tuple (integer)
        """

        if not tuple[1] : # Only non-zero shift sizes are relevant.
            return

        for i in range(len(self.__shift)) : # Look where to insert this tuple
            if self.__shift[i][0] == tuple[0] : # If it already exists,
                self.__shift[i][1] += tuple[1]  # alter it.
                if not self.__shift[i][1] :     # If it results in a zero,
                    self.__shift.pop(i)         # remove it.
                return
            #if

            if self.__shift[i][0] > tuple[0] : # We found a successor, so
                self.__shift.insert(i, tuple)  # insert it before the successor.
                return
            #if
        #for

        self.__shift.append(tuple) # If we couldn't find a successor, so this
                                   # entry will be the last one in the list.
    #__sortins

    def __makeRestrictionList(self, seq) :
        """
        Return a set of restriction enzymes that can bind in a certain
        sequence.

        Private variables:
            - __restrictionBatch ; A RestrictionBatch object.

        @arg seq: The sequence to be analysed
        @type seq: string

        @return: A list of restriction enzymes
        @rtype: list
        """

        restrictionAnalysis = Restriction.Analysis(self.__restrictionBatch, seq)

        d = restrictionAnalysis.with_sites()
        ret = []

        for i in d.keys() :
            for _ in d[i] :
                ret.append(str(i))

        return ret
    #__makeRestrictionSet

    def __restrictionDiff(self, list1, list2) :
        """
        Compare two lists, and count those elements which are only present
        in list1.

        @arg list1: some list
        @type list1: list
        @arg list2: some (other) list
        @type list2: list

        @return: the elements only present in list 1, together with the number
        of occurrences, if more than once present
        @rtype: list
        """

        tempList = list(list1)
        for i in list2 :
            if i in tempList :
                tempList.remove(i)

        ret = []
        tempList.sort()
        for i in set(tempList) :
            c = tempList.count(i)
            if c > 1 :
                ret.append("%s (%i)" % (i, c))
            else :
                ret.append(i)
        #for

        return ret
    #__restrictionDiff

    def __mutate(self, pos1, pos2, ins) :
        """
        A general mutation function that does a delins on interbase
        coordinates of the original string. The change in length (if any)
        is stored by calling the __sortins() function.
        The coordinates are those of the original string, so we use the
        __shifsize() function to map them to the mutated string, on which
        we perform the alteration.

        Private variables:
            - __output ; Visualisation information is added.

        Public variables (altered):
            - mutated ; This string will reflect the result of the given
                        delins.

        @arg pos1:  The first interbase position of the deletion
        @type pos1: integer
        @arg pos2:  The second interbase position of the deletion
        @type pos2: integer
        @arg ins:   The insertion
        @type ins:  string

        @return: visualisation
        @rtype: string
        """

        #
        # This part is for visualisation.
        #

        loflank = self.orig[max(pos1 - config.get('flanksize'), 0):pos1]
        roflank = self.orig[pos2:pos2 + config.get('flanksize')]
        delPart = self.orig[pos1:pos2]
        #odel = delPart
        #if len(odel) > self.__config.maxvissize :
        #    odel = "%s [%ibp] %s" % (odel[:self.__config.flankclipsize],
        #        len(odel) - self.__config.flankclipsize * 2,
        #        odel[-self.__config.flankclipsize:])
        odel = self.visualiseLargeString(delPart)

        bp1 = self.shiftpos(pos1)
        bp2 = self.shiftpos(pos2)
        lmflank = self.mutated[max(bp1 - config.get('flanksize'), 0):bp1]
        rmflank = self.mutated[bp2:bp2 + config.get('flanksize')]

        #insvis = ins
        #if len(ins) > self.__config.maxvissize :
        #    insvis = "%s [%ibp] %s" % (ins[:self.__config.flankclipsize],
        #        len(ins) - self.__config.flankclipsize * 2,
        #        ins[-self.__config.flankclipsize:])
        insvis = self.visualiseLargeString(ins)
        fill = abs(len(odel) - len(insvis))
        if len(odel) > len(ins) :
            visualisation = ["%s %s %s" % (loflank, odel, roflank),
                "%s %s%s %s" % (lmflank, insvis, '-' * fill, rmflank)]
        else :
            visualisation = ["%s %s%s %s" % (loflank, odel, '-' * fill,
                roflank), "%s %s %s" % (lmflank, insvis, rmflank)]

        #
        # End visualisation part.
        #

        #
        # Restriction site analysis:
        #

        list1 = self.__makeRestrictionList(loflank + delPart + roflank)
        list2 = self.__makeRestrictionList(lmflank + ins + rmflank)
        self.__output.addOutput("restrictionSites",
            [self.__restrictionDiff(list2, list1),
             self.__restrictionDiff(list1, list2)])
            #[str(list(set1 - set2))[1:-1], str(list(set2 - set1))[1:-1]])

        #
        # End restriction site analysis:
        #

        self.mutated = self.mutated[:self.shiftpos(pos1)] + ins + \
                       self.mutated[self.shiftpos(pos2):]
        self.__sortins([pos2 + 1, len(ins) + pos1 - pos2])

        return visualisation
    #__mutate

    def visualiseLargeString(self, string) :
        """
        If the length of a sequence is larger than a certain maxvissize, the
        string is clipped; otherwise the string is just returned.

        @arg string: DNA sequence
        @type string: string

        @return: either the original sequence, or an abbreviation of it
        @rtype:  string
        """

        if len(string) > config.get('maxvissize'):
            return "%s [%ibp] %s" % (string[:config.get('flankclipsize')],
                len(string) - config.get('flankclipsize') * 2,
                string[-config.get('flankclipsize'):])
        return string
    #visualiseIns

    def shift_minus_at(self, position):
        """
        Indicates if the position-shift gets smaller at exactly the given
        position.

        @arg  position: Position in the original string.
        @type position: int

        @return: True if the position-shift gets smaller at exactly the
                 given position, False otherwise.
        @rtype:  bool

        @todo: Since the __shift list is sorted we could optimize this a
               bit.
        """
        return reduce(lambda b,s: b or (s[0] == position and s[1] < 0),
                      self.__shift, False)
    #shift_minus_at

    def shiftpos(self, position) :
        """
        Calculate the position in the mutated string, given a position in
        the original string.

        Private variables:
            - __shift ; Used to calculate the shift.

        @arg position:  The position in the original string for which we want the
                        shift size
        @type position: integer

        @return: The position in the mutated string
        @rtype:  integer
        """
        ret = position

        for i in range(len(self.__shift)) :
            if self.__shift[i][0] > position :
                return ret

            ret += self.__shift[i][1]
        #for

        return ret
    #shiftpos

    def add_removed_sites(self, sites):
        """
        Add sites to the set of splice sites to ignore in the mutated string.

        @arg sites:  A list of splice sites to ignore.
        @type sites: list of int

        @todo: Resulting list of ignored sites should always be even.
        @todo: Don't remove CDS start/stop, as happens e.g. with
               AL449423.14(CDKN2A_v002):c.5_400del.
        """
        for site in sites:
            self.__removed_sites.add(site)
    #add_ignore_sites

    def newSplice(self, sites) :
        """
        Generate a list of new splice sites.

        @arg  sites: A list of old splice sites.
        @type sites: list of int

        @return: A list of new splice sites.
        @rtype:  list of int


        Example 1 (DNA): NG_012772.1(BRCA2_v001)

                  ...---------[=========]----------...
                              ^         ^
                            18964     19013

          Variant           Expected new location for splice site 18964
          g.18963del        18963
          g.18964del        18964
          g.18963_18964ins  18964
          g.18964_18965ins  18964

          Variant           Expected new location for splice site 19013
          g.19013del        19012
          g.19014del        19013
          g.19013_19014ins  19014


        Example 2 (RNA): NM_000088.3

                      ...============][==============...
                                     /\
                                  229  230

          Variant           Expected new location for splice sites 229,230
          n.228del          228,229
          n.229del          228,229
          n.230del          229,230
          n.231del          229,230
          n.228_229ins      230,231
          n.229_230ins      229,230 or 230,231
          n.230_231ins      229,230
        """

        # We use shiftpos(i+1)-1 instead of shiftpos(i) (and its mirror)
        # to make sure insertions directly before or after an exon are
        # placed inside the exon.
        #
        # Example:
        #
        #   -----SPLICE[======]SPLICE----------SPLICE[=======]SPLICE-----
        #                      ^                    ^
        #                      ins                  ins
        #
        #   These two insertions should be mapped inside the exons because
        #   they are before and after (respectively) their exons and don't
        #   hit the (biological) splice sites.
        #
        # This also makes sure deletions of the last exon base are really
        # removed from the exon. The problem is that positions following
        # (but not including) the deletion get a shift, but the splice site
        # is stored by the position of the last exon base. So the splice
        # site position would not be decremented without the +1-1 dance.

        new_sites = []

        prev_donor = None
        filtered_sites = filter(lambda s: s not in self.__removed_sites, sites)
        for acceptor, donor in util.grouper(filtered_sites):

            # We don't want to do the -1+1 dance if
            # 1) there is a deletion directly before the exon, or
            # 2) there is another exon directly before this exon, or
            # 3) this is the first site in the list.
            #
            # A consequence of check 2) is that insertions between two
            # directly adjacent exons are seen as insertions in the first
            # exon.
            #
            # Condition 3) makes sure we don't include insertions directly
            # in front of CDS start in the CDS. It also affects translation
            # start, but this should be no problem.
            if not prev_donor or prev_donor == acceptor - 1 or \
                   self.shift_minus_at(acceptor):
                new_sites.append(self.shiftpos(acceptor))
            else:
                new_sites.append(self.shiftpos(acceptor - 1) + 1)

            # Should never happen since splice sites come in pairs.
            if not donor: continue

            # We don't want to do the +1-1 dance if this is the last site
            # in the list. This makes sure we don't include insertions
            # directly at CDS end in the CDS. It also affects translation
            # end, but this should be no problem.
            if donor == sites[-1]:
                new_sites.append(self.shiftpos(donor))
            else:
                new_sites.append(self.shiftpos(donor + 1) - 1)

            prev_donor = donor

        return new_sites
    #newSplice

    def delM(self, pos1, pos2) :
        """
        Delete a range from non-interbase position pos1 to pos2.

        Private variables:
            - __output ; Visualisation information is added.

        @arg pos1:  The first nucleotide of the range to be deleted
        @type pos1: integer
        @arg pos2:  The last nucleotide of the range to be deleted
        @type pos2: integer
        """

        if pos1 == pos2 :
            visualisation = ["deletion of %i" % pos1]
        else :
            visualisation = ["deletion of %i to %i" % (pos1, pos2)]

        visualisation.extend(self.__mutate(pos1 - 1, pos2, ''))
        self.__output.addOutput("visualisation", visualisation)
    #delM

    def insM(self, pos, ins) :
        """
        Insert a string at interbase position pos.

        Private variables:
           -  __output ; Visualisation information is added.

        @arg pos:  The interbase position where the insertion should take place
        @type pos: integer
        @arg ins:  The insertion
        @type ins: string
        """
        visualisation = ["insertion between %i and %i" % (pos, pos + 1)]
        visualisation.extend(self.__mutate(pos, pos, ins))
        self.__output.addOutput("visualisation", visualisation)
    #insM

    def delinsM(self, pos1, pos2, ins) :
        """
        Delete a range from non-interbase position pos1 to pos2 and insert
        ins.

        @arg pos1:  The first nucleotide of the range to be deleted
        @type pos1: integer
        @arg pos2:  The last nucleotide of the range to be deleted.
        @type pos2: integer
        @arg ins:   The insertion
        @type ins:  string
        """

        visualisation = ["delins from %i to %i" % (pos1, pos2)]
        visualisation.extend(self.__mutate(pos1 - 1, pos2, ins))
        self.__output.addOutput("visualisation", visualisation)
    #delinsM

    def subM(self, pos, nuc) :
        """
        Substitute a nucleotide at non-interbase position pos for nuc.

        Private variables:
            - __output ; Visualisation information is added.

       @arg pos:  The position where the substitution should take place
       @type pos: integer
       @arg nuc:  The new nucleotide
       @type nuc: string
        """

        visualisation = ["substitution at %i" % pos]
        visualisation.extend(self.__mutate(pos - 1, pos, nuc))
        self.__output.addOutput("visualisation", visualisation)
    #subM

    def invM(self, pos1, pos2) :
        """
        Invert a range from non-interbase position pos1 to pos2.

        Public variables:
            - orig ; The original string.

        @arg pos1:  The first nucleotide of the range to be inverted
        @type pos1: integer
        @arg pos2:  The last nucleotide of the range to be inverted
        @type pos2: integer
        """

        visualisation = ["inversion between %i and %i" % (pos1, pos2)]
        visualisation.extend(self.__mutate(pos1 - 1, pos2, \
            reverse_complement(self.orig[pos1 - 1:pos2])))
        self.__output.addOutput("visualisation", visualisation)
    #invM

    def dupM(self, pos1, pos2) :
        """
        Duplicate a range from non-interbase position pos1 to pos2.

        Public variables:
            - orig ; The original string.

        @arg pos1:  The first nucleotide of the range to be duplicated
        @type pos1: integer
        @arg pos2:  The last nucleotide of the range to be duplicated
        @type pos2: integer
        """

        visualisation = ["duplication from %i to %i" % (pos1, pos2)]
        visualisation.extend(self.__mutate(pos2, pos2,
                                           self.orig[pos1 - 1:pos2]))
        self.__output.addOutput("visualisation", visualisation)
    #dupM
#Mutator
