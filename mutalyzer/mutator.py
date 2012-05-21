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


from collections import defaultdict

from Bio import Restriction
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import reverse_complement

from mutalyzer import util
from mutalyzer import config


class Mutator():
    """
    Mutate a string and register all shift points. For each mutation a
    visualisation is made (on genomic level) and the addition or deletion
    of restriction sites is detected. Output for each raw variant is stored
    in the output object as 'visualisation', 'deletedRestrictionSites' and
    'addedRestrictionSites' respectively.
    """
    def __init__(self, orig, output):
        """
        Initialise the instance with the original sequence.

        @arg orig: The original sequence before mutation.
        @type orig: str
        @arg output: The output object.
        @type output: mutalyzer.Output.Output
        """
        self._shifts = defaultdict(int)
        self._removed_sites = set()
        self._restriction_batch = Restriction.RestrictionBatch([], ['N'])

        self._output = output
        self.orig = orig

        self.mutated = orig
    #__init__

    def _restriction_count(self, sequence):
        """
        Return the count per restriction enzyme that can bind in a certain
        sequence.

        @arg sequence: The sequence to be analysed
        @type sequence: str

        @return: A mapping of restriction enzymes to counts.
        @rtype: dict
        """
        analysis = Restriction.Analysis(self._restriction_batch, sequence)
        return dict((str(k), len(v)) for k, v in analysis.with_sites().items())
    #_restriction_count

    def _counts_diff(self, counts1, counts2):
        """
        Compare two lists, and count those elements which are only present
        in list1.

        @arg list1: some list
        @type list1: list
        @arg list2: some (other) list
        @type list2: list

        @return: The elements only present in list1, together with the number
            of occurrences if more than once present.
        @rtype: list
        """
        enzymes = set(counts1.keys()) - set(counts2.keys())

        diff = []
        for enzyme in sorted(enzymes):
            count = counts1[enzyme]
            diff.append('%s (%i)' % (enzyme, count) if count > 1 else enzyme)

        return diff
    #_counts_diff

    def _visualise(self, pos1, pos2, ins):
        """
        Create visualisation and do a restriction site analysis on the given
        indel.

        @arg pos1: First interbase position of the deleted sequence.
        @type pos1: int
        @arg pos2: Second interbase position of the deleted sequence.
        @type pos2: int
        @arg ins: Inserted sequence.
        @type ins: str

        @return: Visualisation.
        @rtype: str
        """
        loflank = self.orig[max(pos1 - config.get('flanksize'), 0):pos1]
        roflank = self.orig[pos2:pos2 + config.get('flanksize')]
        delPart = self.orig[pos1:pos2]
        odel = util.visualise_sequence(delPart, config.get('maxvissize'),
                                       config.get('flankclipsize'))

        bp1 = self.shift(pos1)
        bp2 = self.shift(pos2)
        lmflank = self.mutated[max(bp1 - config.get('flanksize'), 0):bp1]
        rmflank = self.mutated[bp2:bp2 + config.get('flanksize')]

        insvis = util.visualise_sequence(ins, config.get('maxvissize'),
                                         config.get('flankclipsize'))
        fill = abs(len(odel) - len(insvis))
        if len(odel) > len(ins):
            visualisation = ['%s %s %s' % (loflank, odel, roflank),
                             '%s %s%s %s' % (lmflank, insvis, '-' * fill, rmflank)]
        else:
            visualisation = ['%s %s%s %s' % (loflank, odel, '-' * fill, roflank),
                             '%s %s %s' % (lmflank, insvis, rmflank)]

        # Todo: This part is for restriction site analysis. It doesn't really
        #     belong in this method, but since it uses many variables computed
        #     for the visualisation, we leave it here for the moment.
        counts1 = self._restriction_count(loflank + delPart + roflank)
        counts2 = self._restriction_count(lmflank + ins + rmflank)
        self._output.addOutput('restrictionSites',
                               [self._counts_diff(counts2, counts1),
                                self._counts_diff(counts1, counts2)])

        return visualisation
    #_visualise

    def _add_shift(self, position, shift):
        """
        Add a shift to the shift list.

        @arg position: Position in the original string.
        @type position: int
        @arg shift: Shift size.
        @type shift: int
        """
        self._shifts[position] += shift
    #_add_shift

    def _shift_minus_at(self, position):
        """
        Indicates if the position-shift gets smaller at exactly the given
        position.

        @arg position: Position in the original string.
        @type position: int

        @return: True if the position-shift gets smaller at exactly the
            given position, False otherwise.
        @rtype: bool
        """
        return self._shifts[position] < 0
    #_shift_minus_at

    def shift_at(self, position):
        """
        Calculate the shift given a position in the original string.

        @arg position: Position in the original string.
        @type position: int

        @return: Shift for the given position.
        @rtype: int
        """
        return sum(s for p, s in self._shifts.items() if p <= position)
    #shift_at

    def shift(self, position):
        """
        Calculate the position in the mutated string, given a position in the
        original string.

        @arg position: Position in the original string.
        @type position: int

        @return: Position in the mutated string.
        @rtype: int
        """
        return position + self.shift_at(position)
    #shift

    def add_removed_sites(self, sites):
        """
        Add sites to the set of splice sites to ignore in the mutated string.

        @arg sites: List of splice sites to ignore.
        @type sites: list(int)

        @todo: Resulting list of ignored sites should always be even.
        @todo: Don't remove CDS start/stop, as happens e.g. with
            AL449423.14(CDKN2A_v002):c.5_400del.
        """
        for site in sites:
            self._removed_sites.add(site)
    #add_removed_sites

    def shift_sites(self, sites):
        """
        Calculate the list of splice sites on the mutated string, given a list
        of splice sites on the original string.

        @arg sites: List of splice sites on the original string.
        @type sites: list(int)

        @return: List of splice sites on the mutated string.
        @rtype: list(int)


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
        filtered_sites = [s for s in sites if s not in self._removed_sites]
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
                    self._shift_minus_at(acceptor):
                new_sites.append(self.shift(acceptor))
            else:
                new_sites.append(self.shift(acceptor - 1) + 1)

            # Should never happen since splice sites come in pairs.
            if not donor: continue

            # We don't want to do the +1-1 dance if this is the last site
            # in the list. This makes sure we don't include insertions
            # directly at CDS end in the CDS. It also affects translation
            # end, but this should be no problem.
            if donor == sites[-1]:
                new_sites.append(self.shift(donor))
            else:
                new_sites.append(self.shift(donor + 1) - 1)

            prev_donor = donor

        return new_sites
    #shift_sites

    def _mutate(self, pos1, pos2, ins):
        """
        A general mutation function that does a delins on interbase
        coordinates of the original string. The change in length (if any)
        is stored in the shift list.

        The coordinates are those of the original string, so we use the shift
        list to map them to the mutated string, on which we perform the
        alteration.

        @arg pos1: First interbase position of the deleted sequence.
        @type pos1: int
        @arg pos2: Second interbase position of the deleted sequence.
        @type pos2: int
        @arg ins: Inserted sequence.
        @type ins: str
        """
        correct = 1 if pos1 == pos2 else 0
        self.mutated = (self.mutated[:self.shift(pos1 + 1) - 1] +
                        ins +
                        self.mutated[self.shift(pos2 + correct) - correct:])

        self._add_shift(pos2 + 1, pos1 - pos2 + len(ins))
    #_mutate

    def deletion(self, pos1, pos2):
        """
        Delete a range from non-interbase position pos1 to pos2.

        @arg pos1: First nucleotide of the deleted sequence.
        @type pos1: int
        @arg pos2: Last nucleotide of the deleted sequence.
        @type pos2: int
        """
        if pos1 == pos2:
            visualisation = ['deletion of %i' % pos1]
        else:
            visualisation = ['deletion of %i to %i' % (pos1, pos2)]

        visualisation.extend(self._visualise(pos1 - 1, pos2, ''))
        self._output.addOutput('visualisation', visualisation)

        self._mutate(pos1 - 1, pos2, '')
    #deletion

    def insertion(self, pos, ins):
        """
        Insert a string at interbase position pos.

        @arg pos: Interbase position where the insertion should take place.
        @type pos: int
        @arg ins: Inserted sequence.
        @type ins: str
        """
        visualisation = ['insertion between %i and %i' % (pos, pos + 1)]
        visualisation.extend(self._visualise(pos, pos, ins))
        self._output.addOutput('visualisation', visualisation)

        self._mutate(pos, pos, ins)
    #insertion

    def delins(self, pos1, pos2, ins):
        """
        Delete a range from non-interbase position pos1 to pos2 and insert
        sequence ins.

        @arg pos1: First nucleotide of the deleted sequence.
        @type pos1: int
        @arg pos2: Last nucleotide of the deleted sequence.
        @type pos2: int
        @arg ins: Inserted sequence.
        @type ins: str
        """
        visualisation = ['delins from %i to %i' % (pos1, pos2)]
        visualisation.extend(self._visualise(pos1 - 1, pos2, ins))
        self._output.addOutput('visualisation', visualisation)

        self._mutate(pos1 - 1, pos2, ins)
    #delins

    def substitution(self, pos, nuc):
        """
        Substitute a nucleotide at non-interbase position pos for nuc.

        @arg pos: Position of the substitution.
        @type pos: int
        @arg nuc: Substituted nucleotide.
        @type nuc: str
        """
        visualisation = ['substitution at %i' % pos]
        visualisation.extend(self._visualise(pos - 1, pos, nuc))
        self._output.addOutput('visualisation', visualisation)

        self._mutate(pos - 1, pos, nuc)
    #substitution

    def inversion(self, pos1, pos2):
        """
        Invert a range from non-interbase position pos1 to pos2.

        @arg pos1: First nucleotide of the inverted sequence.
        @type pos1: int
        @arg pos2: Last nucleotide of the inverted sequence.
        @type pos2: int
        """
        visualisation = ['inversion between %i and %i' % (pos1, pos2)]
        visualisation.extend(
            self._visualise(pos1 - 1, pos2,
                            reverse_complement(self.orig[pos1 - 1:pos2])))
        self._output.addOutput('visualisation', visualisation)

        self._mutate(pos1 - 1, pos2,
                     reverse_complement(self.orig[pos1 - 1:pos2]))
    #inversion

    def duplication(self, pos1, pos2):
        """
        Duplicate a range from non-interbase position pos1 to pos2.

        @arg pos1: First nucleotide of the duplicated sequence.
        @type pos1: int
        @arg pos2: Last nucleotide of the duplicated sequence.
        @type pos2: int
        """
        visualisation = ['duplication from %i to %i' % (pos1, pos2)]
        visualisation.extend(
            self._visualise(pos2, pos2, self.orig[pos1 - 1:pos2]))
        self._output.addOutput('visualisation', visualisation)

        self._mutate(pos1 - 1, pos1 - 1, self.orig[pos1 - 1:pos2])
    #duplication
#Mutator
