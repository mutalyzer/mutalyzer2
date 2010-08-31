#!/usr/bin/python

from Bio import Restriction
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import reverse_complement # reverse_complement()

"""
    Module for mutating a string.

    Mutations are described in the original coordinates. These coordinates are
    transfered to the mutated coordinates with the aid of an internal shift
    list, which keeps track of the sizes of changes. Using the original
    coordinates greatly simplifies combined mutations in a variant. A
    visualisation of each raw variant within a combined variant is made and
    effects on restriction sites are also analysed.

    The original as well as the mutated string are stored here.

    Public classes:
        Mutator ; Mutate a string and register all shift points.
"""

class Mutator() :
    """
        Mutate a string and register all shift points. For each mutation a
        visualisation is made (on genomic level) and the addition or deletion
        of restriction sites is detected. Output for each raw variant is stored
        in the output object as 'visualisation', 'deletedRestrictionSites' and
        'addedRestrictionSites' respectively.

        Private variables:
            __config           ; Configuration variables of this class:

            __output           ; The output object.
            __shift            ; A sorted list of tuples (position, shiftsize)
                                 where the modifications in length are stored.
                                 Each first element of the tuples in this list
                                 is unique, each second element is non-zero.
            __restrictionBatch ;

        Public variables:
            orig    ; The original string.
            mutated ; The mutated string.

        Special methods:
            __init__(orig) ; Initialise the class with the original string.

        Private methods:
            __sortins(tuple)      ; Insert a tuple in a sorted list, after
                                    insertion the list stays sorted.
            __makeRestrictionSet()
            __mutate(pos1, pos2, ins) ; A general mutation function that does a
                                        delins on interbase coordinates of the
                                        original string.

        Public methods:
            shiftpos(position)       ; Calculate the position in the mutated
                                       string given the position in the
                                       original string.
            newSplice(sites)         ; Generate a list of new splice sites.
            delM(pos1, pos2)         ; Delete a range from non-interbase
                                       position pos1 to pos2.
            insM(pos, ins)           ; Insert a string at interbase position
                                       pos.
            delimsM(pos1, pos2, ins) ; Delete a range from non-interbase
                                       position pos1 to pos2 and insert ins.
            subM(pos, nuc)           ; Substitute a nucleotite at non-interbase
                                       position pos for nuc.
            invM(pos1, pos2)         ; Invert a range from non-interbase
                                       position pos1 to pos2.
            dupM(pos1, pos2)         ; Duplicate a range from non-interbase
                                       position pos1 to pos2.
    """

    def __init__(self, orig, config, output) :
        """
            Initialise the class with the original string.

            Arguments:
                orig   ; The original string before mutation.
                config ; Configuration variables.
                output ; The output object.

            Private variables (altered):
                __config           ; Initialised with the configuration
                                     variables.
                __output           ; Initialised with the output object.
                __shift            ; Initialised to the empty list.
                __restrictionBatch ; Initialised to a default set of
                                     restriction enzymes.

            Public variables (altered):
                orig    ; Initialised to the parameter orig.
                mutated ; Initialised to the parameter orig.
        """

        self.__config = config
        self.__output = output
        self.__shift = []
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

            Arguments:
                tuple ; An ordered pair where tuple[0] denotes a position and
                        tuple[1] denotes the change in shift at this position.

            Private variables (altered):
                __shift ; A tuple can be added, removed or altered.
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

            Arguments:
                seq ; The sequence to be analysed.

            Returns:
                list ; A list of restriction enzymes.

            Private variables:
                __restrictionBatch ; A RestrictionBatch object.
        """

        restrictionAnalysis = Restriction.Analysis(self.__restrictionBatch, seq)

        d = restrictionAnalysis.with_sites()
        ret = []

        for i in d.keys() :
            for j in d[i] :
                ret.append(str(i))

        return ret
    #__makeRestrictionSet

    def __restrictionDiff(self, list1, list2) :
        #TODO documentation
        """
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

            Arguments:
                pos1 ; The first interbase position of the deletion.
                pos2 ; The second interbase position of the deletion.
                ins  ; The insertion.

            Private variables:
                __config ; The variables maxvissize, flanksize and flankclipsize
                           are used in the visualisation.
                __output ; Visualisation information is added.

            Public variables (altered):
                mutated ; This string will reflect the result of the given
                          delins.
        """

        #
        # This part is for visualisation.
        #

        loflank = self.orig[max(pos1 - self.__config.flanksize, 0):pos1]
        roflank = self.orig[pos2:pos2 + self.__config.flanksize]
        delPart = self.orig[pos1:pos2]
        #odel = delPart
        #if len(odel) > self.__config.maxvissize :
        #    odel = "%s [%ibp] %s" % (odel[:self.__config.flankclipsize],
        #        len(odel) - self.__config.flankclipsize * 2,
        #        odel[-self.__config.flankclipsize:])
        odel = self.visualiseLargeString(delPart)

        bp1 = self.shiftpos(pos1)
        bp2 = self.shiftpos(pos2)
        lmflank = self.mutated[max(bp1 - self.__config.flanksize, 0):bp1]
        rmflank = self.mutated[bp2:bp2 + self.__config.flanksize]

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
        self.__sortins([pos1 + 1, len(ins) + pos1 - pos2])

        return visualisation
    #__mutate

    def visualiseLargeString(self, string) :
        """
        """

        if len(string) > self.__config.maxvissize :
            return "%s [%ibp] %s" % (string[:self.__config.flankclipsize],
                len(string) - self.__config.flankclipsize * 2,
                string[-self.__config.flankclipsize:])
        return string
    #visualiseIns

    def shiftpos(self, position) :
        """
            Calculate the position in the mutated string, given a position in
            the original string.

            Arguments:
                position ; The position in the original string for which we
                           want the shift size.

            Private variables:
                __shift ; Used to calculate the shift.

            Returns:
                integer ; The position in the mutated string.
        """
        ret = position

        for i in range(len(self.__shift)) :
            if self.__shift[i][0] > position :
                return ret

            ret += self.__shift[i][1]
        #for

        return ret
    #shiftpos

    def newSplice(self, sites) :
        """
            Generate a list of new splice sites.

            Arguments:
                sites ; A list of old splice sites.

            Returns:
                list ; A list of new splice sites.
        """

        ret = []

        j = 0
        for i in sites :
            if (j % 2) :
                #ret.append(self.shiftpos(i + 1) - 1)
                ret.append(self.shiftpos(i))
            else :
                #ret.append(self.shiftpos(i - 1) + 1)
                ret.append(self.shiftpos(i))
            j += 1
        #for

        return ret
    #newSplice

    def delM(self, pos1, pos2) :
        """
            Delete a range from non-interbase position pos1 to pos2.

            Arguments:
                pos1 ; The first nucleotide of the range to be deleted.
                pos2 ; The last nucleotide of the range to be deleted.

            Private variables:
                __output ; Visualisation information is added.
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

            Arguments:
                pos ; The interbase position where the insertion should take
                      place.
                ins ; The insertion, a string.

            Private variables:
                __output ; Visualisation information is added.
        """

        visualisation = ["insertion between %i and %i" % (pos, pos + 1)]
        visualisation.extend(self.__mutate(pos, pos, ins))
        self.__output.addOutput("visualisation", visualisation)
    #insM

    def delinsM(self, pos1, pos2, ins) :
        """
            Delete a range from non-interbase position pos1 to pos2 and insert
            ins.

            Arguments:
                pos1 ; The first nucleotide of the range to be deleted.
                pos2 ; The last nucleotide of the range to be deleted.
                ins  ; The insertion, a string.
        """

        visualisation = ["delins from %i to %i" % (pos1, pos2)]
        visualisation.extend(self.__mutate(pos1 - 1, pos2, ins))
        self.__output.addOutput("visualisation", visualisation)
    #delinsM

    def subM(self, pos, nuc) :
        """
            Substitute a nucleotite at non-interbase position pos for nuc.

            Arguments:
                pos ; The position where the substitution should take place.
                nuc ; The new nucleotide.

            Private variables:
                __output ; Visualisation information is added.
        """

        visualisation = ["substitution at %i" % pos]
        visualisation.extend(self.__mutate(pos - 1, pos, nuc))
        self.__output.addOutput("visualisation", visualisation)
    #subM

    def invM(self, pos1, pos2) :
        """
            Invert a range from non-interbase position pos1 to pos2.

            Arguments:
                pos1 ; The first nucleotide of the range to be inverted.
                pos2 ; The last nucleotide of the range to be inverted.

            Public variables:
                orig ; The original string.
        """

        visualisation = ["inversion between %i and %i" % (pos1, pos2)]
        visualisation.extend(self.__mutate(pos1 - 1, pos2, \
            reverse_complement(self.orig[pos1 - 1:pos2])))
        self.__output.addOutput("visualisation", visualisation)
    #invM

    def dupM(self, pos1, pos2) :
        """
            Duplicate a range from non-interbase position pos1 to pos2.

            Arguments:
                pos1 ; The first nucleotide of the range to be duplicated.
                pos2 ; The last nucleotide of the range to be duplicated.

            Public variables:
                orig ; The original string.
        """

        visualisation = ["duplication from %i to %i" % (pos1, pos2)]
        visualisation.extend(self.__mutate(pos2, pos2,
                                           self.orig[pos1 - 1:pos2]))
        self.__output.addOutput("visualisation", visualisation)
    #dupM
#Mutator

#
# Unit test.
#
if __name__ == "__main__" :
    pass
#if
