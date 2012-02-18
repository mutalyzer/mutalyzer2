"""
The HGVS variant nomenclature checker.

Entrypoint is the check_variant() function.

Notes about naming positions:
* CDS -> use start/stop
* splice sites or exons -> acceptor/donor
* translation -> begin/end
* any range of bases -> first/last
* interbase position (if two numbers are used) -> before/after
"""


from operator import itemgetter, attrgetter

import Bio
import Bio.Seq
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from mutalyzer import util
from mutalyzer.grammar import Grammar
from mutalyzer.mutator import Mutator
from mutalyzer.mapping import Converter
from mutalyzer import Retriever
from mutalyzer import GenRecord
from mutalyzer import Db


# Exceptions used (privately) in this module.
class _VariantError(Exception): pass
class _RawVariantError(_VariantError): pass
class _UnknownPositionError(_RawVariantError): pass
class _NotDNAError(_RawVariantError): pass
class _PositionsNotConsecutiveError(_RawVariantError): pass
class _LengthMismatchError(_RawVariantError): pass
class _ReferenceMismatchError(_RawVariantError): pass
class _RangeInsertionError(_RawVariantError): pass
class _OffsetSignError(_RawVariantError):
    def __init__(self, main, offset, acceptor):
        self.main = main
        self.offset = offset
        self.acceptor = acceptor
class _OffsetNotFromBoundaryError(_RawVariantError):
    def __init__(self, main):
        self.main = main
class _InvalidExonError(_RawVariantError):
    def __init__(self, exon):
        self.exon = exon
class _InvalidIntronError(_RawVariantError):
    def __init__(self, intron):
        self.intron = intron


def _is_coding_intronic(loc):
    """
    Check whether a location is an intronic c. position.

    @arg loc: A location from the Parser module.
    @type loc: pyparsing.ParseResults

    @return: True if the location is c. intronic, False otherwise.
    @rtype: boolean
    """
    if not loc:
        return False
    if not loc.PtLoc:
        return False
    if not loc.PtLoc.Offset:
        return False
    return True
#_is_coding_intronic


def _check_intronic_position(main, offset, transcript):
    """
    Check whether a c. position is really in an intron: The main coordinate
    must be a splice site and the offset coordinate must have the correct
    sign. Raise _RawVariantError exception if this check fails.

    @arg main: Main coordinate of the position.
    @type main: int
    @arg offset: Offset coordinate of the position.
    @type offset: int
    @arg transcript: Transcript under scrutiny.
    @type transcript: object

    @raise _OffsetNotFromBoundary: An offset from a non-boundary position
                                   was used.
    @raise _OffsetSignError: Offset from exon boundary has the wrong sign.

    @todo: Check if the offset is really in the flanking intron.
    """
    main_g = transcript.CM.x2g(main, 0)
    sites = transcript.CM.RNA

    if offset:
        oriented_offset = offset * transcript.CM.orientation
        try:
            i = sites.index(main_g)
            if not i % 2:
                # Splice acceptor, so sign must be -.
                if oriented_offset > 0:
                    raise _OffsetSignError(
                        transcript.CM.int2main(main),
                        transcript.CM.int2offset((main, offset)),
                        True)
            else:
                # Splice donor, so sign must be +.
                if oriented_offset < 0:
                    raise _OffsetSignError(
                        transcript.CM.int2main(main),
                        transcript.CM.int2offset((main, offset)),
                        False)
        except ValueError:
            # The main coordinate is not a splice site.
            raise _OffsetNotFromBoundaryError(transcript.CM.int2main(main))
#_check_intronic_position


def _check_argument(argument, reference, first, last, output):
    """
    Do several checks for the optional argument of a variant. Raise a
    _RawVariantError exception if the checks fail.

    @arg reference: The reference sequence.
    @type reference: string
    @arg first: Start position of the variant.
    @type first: int
    @arg last: End position of the variant.
    @type last: int
    @arg argument: The optional argument.
    @type argument: string
    @arg output: The Output object.
    @type output: mutalyzer.Output.Output

    @raise _LengthMismatchError: The argument is a length, but it does not
                                 match the given range length.
    @raise _NotDNAError: The argument should be DNA, but it is not.
    @raise _ReferenceMismatchError: The argument is DNA, but it does not
                                    match the given reference.
    """
    if not argument:
        # The argument is optional, if it is not present, it is correct.
        return

    if argument.isdigit():
        # If it is a digit (3_9del7 for example), the digit must be equal to
        # the length of the given range.
        length = int(argument)
        interval = first - last + 1
        if length != interval:
            output.addMessage(__file__, 3, 'EARGLEN',
                              'The length (%i) differed from that of the ' \
                              'range (%i).' % (length, interval))
            raise _LengthMismatchError()
    else:
        # If it is not a digit, it muse be DNA.
        if not util.is_dna(argument):
            output.addMessage(__file__, 4, 'ENODNA',
                              'Invalid letters in argument.')
            raise _NotDNAError()
        # And the DNA must match the reference sequence.
        reference_slice = str(reference[first - 1:last])
        if reference_slice != str(argument):
            # Todo: Be more informative.
            output.addMessage(__file__, 3, 'EREF',
                              '%s not found at position %s, found %s ' \
                              'instead.' % (argument,
                                            util.format_range(first, last),
                                            reference_slice))
            raise _ReferenceMismatchError()
#_check_argument


def _add_batch_output(O):
    """
    Format the results to a batch output.

    Filter the mutalyzer output and reformat it for use in the batch system
    as output object 'batchDone'.

    @arg O: The Output object
    @type O: Modules.Output.Output

    @todo: More documentation.
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

    # Add effects on restriction sites as two columns (created and deleted).
    # The value for each column is a semicolon-separated list of
    # comma-separated lists: for each raw variant, a list of restriction
    # sites.
    sites_created = []
    sites_deleted = []
    for variant in O.getOutput("restrictionSites"):
        sites_created.append(','.join(variant[0]))
        sites_deleted.append(','.join(variant[1]))
    outputline += "%s\t%s\t" % (';'.join(sites_created), ';'.join(sites_deleted))

    #Link naar additional info:
    #outputline+="http://localhost/mutalyzer2/redirect?mutationName=%s" %\
    #        "todovariant"

    O.addOutput("batchDone", outputline)
#_add_batch_output


def apply_substitution(position, original, substitute, mutator, record, O):
    """
    Do a semantic check for a substitution, do the actual substitution and
    give it a name.

    @arg position: Genomic location of the substitution.
    @type position: int
    @arg original: Nucleotide in the reference sequence.
    @type original: string
    @arg substitute: Nucleotide in the mutated sequence.
    @type substitute: string
    @arg mutator: A Mutator instance.
    @type mutator: mutalyzer.mutator.Mutator
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg O: The Output object.
    @type O: Modules.Output.Output

    @raise _NotDNAError: Invalid (non-DNA) letter in input.
    """
    if not util.is_dna(substitute):
        O.addMessage(__file__, 4, 'ENODNA', 'Invalid letter in input')
        raise _NotDNAError()

    if original == substitute:
        O.addMessage(__file__, 2, 'WNOCHANGE',
                     'No mutation given (%c>%c) at position %i.' % \
                     (original, substitute, position))
        return

    mutator.subM(position, substitute)

    record.name(position, position, 'subst', mutator.orig[position - 1],
                substitute, None)
#apply_substitution


def apply_deletion_duplication(first, last, type, mutator, record, O,
                               first_fuzzy=False, last_fuzzy=False):
    """
    Do a semantic check for a deletion or duplication, do the actual
    deletion/duplication and give it a name.

    @arg first: Genomic start position of the del/dup.
    @type first: int
    @arg last: Genomic end position of the del/dup.
    @type last: int
    @arg type: The variant type (del or dup).
    @type type: string
    @arg mutator: A Mutator instance.
    @type mutator: mutalyzer.mutator.Mutator
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg O: The Output object.
    @type O: Modules.Output.Output

    @kwarg first_fuzzy: Denotes that the start position is fuzzy (e.g. in the
        case of an unknown offset in c. notation).
    @type first_fuzzy: bool
    @kwarg last_fuzzy: Denotes that the end position is fuzzy (e.g. in the
        case of an unknown offset in c. notation).
    @type last_fuzzy: bool
    """
    reverse_roll, forward_roll = util.roll(mutator.orig, first, last)

    # In the case of RNA, check if we roll over a splice site. If so, make
    # the roll shorter, just up to the splice site. (Effectively this always
    # means we roll over two splice sites, since they are adjacent.)
    # We only have to consider the forward roll, since RNA reference
    # sequences are always orientated in correspondence with the transcript.
    original_forward_roll = forward_roll
    if record.record.molType != 'g':
        # Todo: Do we assume .geneList[0].transcriptList[0] is the selected
        # transcript here?? Why not use record.current_transcript?
        splice_sites = record.record.geneList[0].transcriptList[0] \
                       .mRNA.positionList
        for acceptor, donor in util.grouper(splice_sites):
            # Note that acceptor and donor splice sites both point to the
            # first, respectively last, position of the exon, so they are
            # both at different sides of the boundary.
            if last < acceptor and last + forward_roll >= acceptor:
                forward_roll = acceptor - 1 - last
                break
            if last <= donor and last + forward_roll > donor:
                forward_roll = donor - last
                break

    # Did we select a transcript on the reverse strand?
    transcript = record.current_transcript()
    reverse_strand = transcript and transcript.CM.orientation == -1

    if forward_roll and not reverse_strand:
        new_first = first + forward_roll
        new_stop = last + forward_roll
        O.addMessage(__file__, 2, 'WROLLFORWARD',
            'Sequence "%s" at position %s was given, however, ' \
            'the HGVS notation prescribes that on the forward strand ' \
            'it should be "%s" at position %s.' % (
            mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
            util.format_range(first, last),
            mutator.visualiseLargeString(str(mutator.orig[new_first - 1:new_stop])),
            util.format_range(new_first, new_stop)))

    if forward_roll != original_forward_roll and not reverse_strand:
        # The original roll was decreased because it crossed a splice site.
        incorrect_first = first + original_forward_roll
        incorrect_stop = last + original_forward_roll
        O.addMessage(__file__, 1, 'IROLLBACK',
            'Sequence "%s" at position %s was not corrected to "%s" at ' \
            'position %s, since they reside in different exons.' % (
            mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
            util.format_range(first, last),
            mutator.visualiseLargeString(str(mutator.orig[incorrect_first - 1:incorrect_stop])),
            util.format_range(incorrect_first, incorrect_stop)))

    if reverse_roll and reverse_strand:
        new_first = first - reverse_roll
        new_stop = last - reverse_roll
        O.addMessage(__file__, 2, 'WROLLREVERSE',
            'Sequence "%s" at position %s was given, however, ' \
            'the HGVS notation prescribes that on the reverse strand ' \
            'it should be "%s" at position %s.' % (
            mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
            util.format_range(first, last),
            mutator.visualiseLargeString(str(mutator.orig[new_first - 1:new_stop])),
            util.format_range(new_first, new_stop)))

    # We don't go through the trouble of visualising the *corrected* variant
    # and are happy with visualising what the user gave us.
    if type == 'del':
        mutator.delM(first, last)
    else:
        mutator.dupM(first, last)

    record.name(first, last, type, '', '', (reverse_roll, forward_roll),
                start_fuzzy=first_fuzzy,
                stop_fuzzy=last_fuzzy)
#apply_deletion_duplication


def apply_inversion(first, last, mutator, record, O):
    """
    Do a semantic check for an inversion, do the actual inversion, and give
    it a name.

    @arg first: Genomic start position of the inversion.
    @type first: int
    @arg last: Genomic end position of the inversion.
    @type last: int
    @arg mutator: A Mutator instance.
    @type mutator: mutalyzer.mutator.Mutator
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg O: The Output object.
    @type O: Modules.Output.Output
    """
    snoop = util.palinsnoop(mutator.orig[first - 1:last])

    if snoop:
        # We have a reverse-complement-palindromic prefix.
        if snoop == -1 :
            # Actually, not just a prefix, but the entire selected sequence is
            # a 'palindrome'.
            O.addMessage(__file__, 2, 'WNOCHANGE',
                'Sequence "%s" at position %i_%i is a palindrome ' \
                '(its own reverse complement).' % (
                mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
                first, last))
            return
        else:
            O.addMessage(__file__, 2, 'WNOTMINIMAL',
                'Sequence "%s" at position %i_%i is a partial ' \
                'palindrome (the first %i nucleotide(s) are the reverse ' \
                'complement of the last one(s)), the HGVS notation ' \
                'prescribes that it should be "%s" at position %i_%i.' % (
                mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
                first, last, snoop,
                mutator.visualiseLargeString(
                    str(mutator.orig[first + snoop - 1: last - snoop])),
                first + snoop, last - snoop))
            first += snoop
            last -= snoop

    mutator.invM(first, last)

    if first == last:
        O.addMessage(__file__, 2, 'WWRONGTYPE', 'Inversion at position ' \
            '%i is actually a substitution.' % first)
        record.name(first, first, 'subst', mutator.orig[first - 1],
            Bio.Seq.reverse_complement(mutator.orig[first - 1]), None)
    else :
        record.name(first, last, 'inv', '', '', None)
#apply_inversion


def apply_insertion(before, after, s, mutator, record, O):
    """
    Do a semantic check for an insertion, do the actual insertion, and give
    it a name.

    @arg before: Genomic position before the insertion.
    @type before: int
    @arg after: Genomic position after the insertion.
    @type after: int
    @arg s: Nucleotides to be inserted.
    @type s: string
    @arg mutator: A Mutator instance.
    @type mutator: mutalyzer.mutator.Mutator
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg O: The Output object.
    @type O: Modules.Output.Output

    @raise _NotDNAError: Invalid (non-DNA) letter in input.
    @raise _PositionsNotConsecutiveError: Positions {before} and {after} are
                                          not consecutive.
    """
    if before + 1 != after:
        O.addMessage(__file__, 3, 'EINSRANGE',
            '%i and %i are not consecutive positions.' % (before, after))
        raise _PositionsNotConsecutiveError()

    if not s or not util.is_dna(s):
        O.addMessage(__file__, 3, 'EUNKVAR', 'Although the syntax of this ' \
            'variant is correct, the effect can not be analysed.')
        raise _NotDNAError()

    insertion_length = len(s)

    # We don't go through the trouble of visualising the *corrected* variant
    # and are happy with visualising what the user gave us.
    mutator.insM(before, s)

    new_before = mutator.shiftpos(before)
    new_stop = mutator.shiftpos(before) + insertion_length

    reverse_roll, forward_roll = util.roll(mutator.mutated, new_before + 1, new_stop)

    # In the case of RNA, check if we roll over a splice site. If so, make
    # the roll shorter, just up to the splice site. (Effectively this always
    # means we roll over two splice sites, since they are adjacent.)
    # We only have to consider the forward roll, since RNA reference
    # sequences are always orientated in correspondence with the transcript.
    original_forward_roll = forward_roll
    if record.record.molType != 'g' :
        splice_sites = record.record.geneList[0].transcriptList[0] \
                       .mRNA.positionList
        for acceptor, donor in util.grouper(splice_sites):
            # Note that acceptor and donor splice sites both point to the
            # first, respectively last, position of the exon, so they are
            # both at different sides of the boundary.
            if new_stop < acceptor and new_stop + forward_roll >= acceptor:
                forward_roll = acceptor - 1 - new_stop
                break
            if new_stop <= donor and new_stop + forward_roll > donor:
                forward_roll = donor - new_stop
                break

    if reverse_roll + forward_roll >= insertion_length:
        # Todo: Could there also be a IROLLBACK message in this case?
        O.addMessage(__file__, 2, 'WINSDUP',
            'Insertion of %s at position %i_%i was given, ' \
            'however, the HGVS notation prescribes that it should be a ' \
            'duplication of %s at position %i_%i.' % (
            s, before, before + 1,
            mutator.mutated[new_before + forward_roll:new_stop + forward_roll],
            before + forward_roll,
            before + forward_roll + insertion_length - 1))
        after += forward_roll - 1
        before = after - insertion_length + 1
        record.name(before, after, 'dup', '', '',
                    (reverse_roll + forward_roll - insertion_length, 0))
        return

    # Did we select a transcript on the reverse strand?
    transcript = record.current_transcript()
    reverse_strand = transcript and transcript.CM.orientation == -1

    if forward_roll and not reverse_strand:
        O.addMessage(__file__, 2, 'WROLLFORWARD', 'Insertion of %s at position ' \
            '%i_%i was given, however, the HGVS notation prescribes ' \
            'that on the forward strand it should be an insertion of %s ' \
            'at position %i_%i.' % (
            s, before, before + 1,
            mutator.mutated[new_before + forward_roll:new_stop + forward_roll],
            new_before + forward_roll, new_before + forward_roll + 1))

    if forward_roll != original_forward_roll and not reverse_strand:
        # The original roll was decreased because it crossed a splice site.
        O.addMessage(__file__, 1, 'IROLLBACK',
            'Insertion of %s at position %i_%i was not corrected to an ' \
            'insertion of %s at position %i_%i, since they reside in ' \
            'different exons.' % (
            s, before, before + 1,
            mutator.mutated[new_before + original_forward_roll:new_stop + original_forward_roll],
            new_before + original_forward_roll, new_before + original_forward_roll + 1))

    if reverse_roll and reverse_strand:
        O.addMessage(__file__, 2, 'WROLLREVERSE', 'Insertion of %s at position ' \
            '%i_%i was given, however, the HGVS notation prescribes ' \
            'that on the reverse strand it should be an insertion of %s ' \
            'at position %i_%i.' % (
            s, before, before + 1,
            mutator.mutated[new_before - reverse_roll:new_stop - reverse_roll],
            new_before - reverse_roll, (new_before - reverse_roll) + 1))

    record.name(before, before + 1, 'ins',
                mutator.mutated[new_before + forward_roll:new_stop + forward_roll],
                '', (reverse_roll, forward_roll),
                mutator.mutated[new_before - reverse_roll:new_stop - reverse_roll])
#apply_insertion


def apply_delins(first, last, delete, insert, mutator, record, output):
    """
    Do a semantic check for an delins, do the actual delins, and give
    it a name.

    @arg first: Genomic start position of the delins.
    @type first: int
    @arg last: Genomic end position of the delins.
    @type last: int
    @arg delete: Sequence to delete (may be None, in which case it will be
                 constructed from the reference sequence).
    @type delete: string
    @arg insert: Sequence to insert.
    @type insert: string
    @arg mutator: A Mutator instance.
    @type mutator: mutalyzer.mutator.Mutator
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg output: The Output object.
    @type output: Modules.Output.Output
    """
    if not delete:
        delete = mutator.orig[first - 1:last]

    if str(delete) == str(insert):
        output.addMessage(__file__, 2, 'WNOCHANGE',
                          'Sequence "%s" at position %i_%i is identical to ' \
                          'the variant.' % (
                mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
                              first, last))
        return

    delete_trimmed, insert_trimmed, lcp, lcs = util.trim_common(delete, insert)

    if not len(delete_trimmed):
        output.addMessage(__file__, 2, 'WWRONGTYPE', 'The given DelIns ' \
                          'is actually an insertion.')
        apply_insertion(first + lcp - 1, first + lcp, insert_trimmed, mutator,
                        record, output)
        return

    if len(delete_trimmed) == 1 and len(insert_trimmed) == 1:
            output.addMessage(__file__, 2, 'WWRONGTYPE', 'The given DelIns ' \
                              'is actually a substitution.')
            apply_substitution(first + lcp, delete_trimmed, insert_trimmed,
                               mutator, record, output)
            return

    if not len(insert_trimmed):
        output.addMessage(__file__, 2, 'WWRONGTYPE', 'The given DelIns ' \
                          'is actually a deletion.')
        apply_deletion_duplication(first + lcp, last - lcs, 'del',
                                   mutator, record, output)
        return

    if str(Bio.Seq.reverse_complement(delete_trimmed)) == insert_trimmed:
        output.addMessage(__file__, 2, 'WWRONGTYPE', 'The given DelIns ' \
                          'is actually an inversion.')
        apply_inversion(first + lcp, last - lcs, mutator,
                        record, output)
        return

    if len(insert) != len(insert_trimmed):
        output.addMessage(__file__, 2, 'WNOTMINIMAL',
                'Sequence "%s" at position %i_%i has the same prefix or ' \
                'suffix as the inserted sequence "%s". The HGVS notation ' \
                'prescribes that it should be "%s" at position %i_%i.' % (
                mutator.visualiseLargeString(str(mutator.orig[first - 1:last])),
                first, last, insert, insert_trimmed, first + lcp, last - lcs))

    mutator.delinsM(first + lcp, last - lcs, insert_trimmed)

    record.name(first + lcp, last - lcs, 'delins', insert_trimmed, '', None)
#apply_delins


def _get_offset(location, main_genomic, sites, output):
    """
    Convert the offset coordinate in a location (from the Parser) to an
    integer.

    @arg location: A location.
    @type location: pyparsing.ParseResults
    @arg main_genomic: Genomic main position to which the offset belongs.
    @type main_genomic: int
    @arg sites: List of splice sites.
    @type sites: list
    @arg output: The Output object.
    @type output: Modules.Output.Output

    @return: Integer representation of the offset coordinate.
    @rtype: int
    """
    if location.Offset :
        if location.Offset == '?' :
            try:
                # Todo: If it removes CDS start, don't do protein translation.
                # Todo: Wrt orientation, perhaps always go to splice site
                #     locations via the crossmapper...
                # Todo: Also check if +? and -? are correctly used.
                # Todo: Exactly centering might not be so nice, since the center
                # might be closer to a neighbouring exon, making a+xxx from b-?
                # and vice versa. This might not be fixed directly by doing a
                # center +/- 1 because there might be rolling. Ideally we
                # disable rolling entirely for these positions...
                #
                # Note that the code below might be a bit confusing, especially
                # considering reverse strand transcripts. Magically, it works
                # for both orientations.
                i = sites.index(main_genomic)
                if i == 0:
                    # Before first exon (or last on the reverse strand).
                    offset = main_genomic / 2
                elif i == len(sites) - 1:
                    # After last exon (or first on the reverse strand).
                    # Todo: Get length of reference, and calculate a sensible
                    # offset.
                    #
                    # We now use that 2000 is the default downstream length,
                    # but of course this is bogus on the reverse strand and
                    # just a hack anyway.
                    offset = 1000
                elif i % 2 == 0:
                    # Acceptor site (or donor on the reverse strand).
                    offset = abs(main_genomic - sites[i - 1]) / 2 - 1
                else:
                    # Donor site (or acceptor on the reverse strand).
                    offset = abs(sites[i + 1] - main_genomic) / 2 - 1
                # Todo: We would like to use the c. position in this message.
                output.addMessage(__file__, 1, "IUNKNOWNOFFSET", "Unknown offset " \
                                  "relative to %s interpreted as middle of " \
                                  "flanking intron." % main_genomic)
            except ValueError:
                # Todo: This means we don't get an error if the main position
                # was not on an exon boundary. We should return something else
                # than 0 I guess.
                #return 0  # This is highly debatable.
                # Any non-zero value will do.
                return 1
        else:
            offset = int(location.Offset)
        if location.OffSgn == '-' :
            return -offset
        return offset

    return 0
#_get_offset


def _intronic_to_genomic(location, transcript):
    """
    Get genomic location from IVS location.

    @arg location: A location.
    @type location: pyparsing.ParseResults
    @arg transcript: todo
    @type transcript: todo

    @return: Genomic location represented by given IVS location.
    @rtype: int

    @raise _InvalidIntronError: Intron does not exist.
    """
    ivs_number = int(location.IVSNumber)

    if ivs_number < 1 or ivs_number > transcript.CM.numberOfIntrons():
        raise _InvalidIntronError(ivs_number)

    if location.OffSgn == '+':
        return transcript.CM.getSpliceSite(ivs_number * 2 - 1) + \
               transcript.CM.orientation * int(location.Offset)
    else:
        return transcript.CM.getSpliceSite(ivs_number * 2) - \
               transcript.CM.orientation * int(location.Offset)
#_intronic_to_genomic


def _exonic_to_genomic(location, transcript) :
    """
    Get genomic range from EX location.

    @arg location: A location.
    @type location: pyparsing.ParseResults
    @arg transcript: todo
    @type transcript: todo

    @return: A tuple of:
        - first: Genomic start location represented by given EX location.
        - last:  Genomic end location represented by given EX location.
    @rtype: tuple(int, int)

    @raise _InvalidExonError: Exon does not exist.

    @todo: We probably want to treat this as a-?_b+?, so take the centers of
           flanking exons.
    """
    first_exon = int(location.EXNumberStart)
    if first_exon < 1 or first_exon > transcript.CM.numberOfExons():
        raise _InvalidExonError(first_exon)
    first = transcript.CM.getSpliceSite(first_exon * 2 - 2)

    if location.EXNumberStop:
        last_exon = int(location.EXNumberStop)
        if last_exon < 1 or last_exon > transcript.CM.numberOfExons():
            raise _InvalidExonError(last_exon)
        last = transcript.CM.getSpliceSite(last_exon * 2 - 1)
    else:
        last = transcript.CM.getSpliceSite(first_exon * 2 - 1)

    return first, last
#_exonic_to_genomic


def _genomic_to_genomic(first_location, last_location):
    """
    Get genomic range from parsed genomic location.

    @arg first_location: The start location (g.) of the variant.
    @type first_location: pyparsing.ParseResults
    @arg last_location: The start location (g.) of the variant.
    @type last_location: pyparsing.ParseResults

    @return: A tuple of:
        - first: Genomic start location represented by given location.
        - last:  Genomic end location represented by given location.
    @rtype: tuple(int, int)

    @raise _UnknownPositionError: Unknown positions were used.
    @raise _RawVariantError: Range cannot be intepreted.
    """
    if not first_location.Main or not last_location.Main:
        # Unknown positions are denoted by the '?' character.
        raise _UnknownPositionError()

    if not first_location.Main.isdigit() or not last_location.Main.isdigit():
        raise _RawVariantError()

    first = int(first_location.Main)
    last = int(last_location.Main)

    return first, last


def _coding_to_genomic(first_location, last_location, transcript, output):
    """
    Get genomic range from parsed c. location.

    @arg first_location: The start location (c.) of the variant.
    @type first_location: pyparsing.ParseResults
    @arg last_location: The start location (c.) of the variant.
    @type last_location: pyparsing.ParseResults
    @arg transcript: todo
    @type transcript: todo
    @arg output: The Output object.
    @type output: Modules.Output.Output

    @return: A tuple of:
        - first: Genomic start location represented by given location.
        - last:  Genomic end location represented by given location.
    @rtype: tuple(int, int)

    @raise _UnknownPositionError: Unknown positions were used.
    @raise _RawVariantError: Range cannot be interpreted.
    @raise _OffsetNotFromBoundary: An offset from a non-boundary position
                                    was used.
    @raise _OffsetSignError: Offset from exon boundary has the wrong sign.
    """
    if not first_location.Main or not last_location.Main:
        # Unknown positions are denoted by the '?' character.
        raise _UnknownPositionError()

    if not first_location.Main.isdigit() or not last_location.Main.isdigit():
        raise _RawVariantError()

    first_main = transcript.CM.main2int(first_location.MainSgn + \
                                        first_location.Main)
    first_main_genomic = transcript.CM.x2g(first_main, 0)
    first_offset = _get_offset(first_location, first_main_genomic,
                               transcript.CM.RNA, output)

    last_main = transcript.CM.main2int(last_location.MainSgn + \
                                       last_location.Main)
    last_main_genomic = transcript.CM.x2g(last_main, 0)
    last_offset = _get_offset(last_location, last_main_genomic,
                              transcript.CM.RNA, output)

    # These raise _RawVariantError exceptions on invalid positions.
    _check_intronic_position(first_main, first_offset, transcript)
    _check_intronic_position(last_main, last_offset, transcript)

    first = transcript.CM.x2g(first_main, first_offset)
    last = transcript.CM.x2g(last_main, last_offset)

    if transcript.CM.orientation == -1:
        first, last = last, first

    return first, last
#_coding_to_genomic


def process_raw_variant(mutator, variant, record, transcript, output):
    """
    Process a raw variant.

    We raise _RawVariantError if there was something seriously in error
    with the raw variant (but it is still okay to process other raw
    variants). We might (don't at the moment) also raise _VariantError,
    meaning to stop processing the entire variant.

    @arg mutator: A Mutator instance.
    @type mutator: mutalyzer.mutator.Mutator
    @arg variant: A parsed raw (simple, noncompound) variant.
    @type variant: pyparsing.ParseResults
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg transcript: A transcript object.
    @type transcript: Modules.GenRecord.Locus
    @arg output: The Output object.
    @type output: Modules.Output.Output

    @raise _RawVariantError: Cannot process this raw variant.
    @raise _VariantError: Cannot further process the entire variant.
    """
    variant, original_description = variant.RawVar, variant[-1]

    # {argument} may be a number, or a subsequence of the reference.
    # {sequence} is the variant subsequence.
    argument = variant.Arg1
    sequence = variant.Arg2

    # If we are on the reverse strand, subsequences must be in reverse
    # complement.
    if transcript and transcript.CM.orientation == -1:
        sequence = Bio.Seq.reverse_complement(sequence)
        if util.is_dna(argument):
            argument = Bio.Seq.reverse_complement(argument)

    # Get genomic first and last positions for this variant. Below we handle
    # the different ways of describing these positions.

    if variant.EXLoc:
        # EX positioning.
        try:
            first, last = _exonic_to_genomic(variant.EXLoc, transcript)
        except _InvalidExonError as e:
            output.addMessage(__file__, 4, 'EINVALIDEXON',
                              'Non-existing exon number %d given.' % e.exon)
            raise
        if last < first:
            # Todo: Why could this ever happen?
            first, last = last, first

    elif not variant.StartLoc:
        # All non-EX positioning ways need a start location.
        # Todo: Better message.
        output.addMessage(__file__, 4, 'EUNKNOWN',
                          'An unknown error occurred.')
        raise _RawVariantError()

    elif variant.StartLoc.IVSLoc:
        # IVS positioning.
        if record.record.molType != 'g':
            output.addMessage(__file__, 3, 'ENOINTRON', 'Intronic ' \
                'position given for a non-genomic reference sequence.')
            raise _RawVariantError()
        try:
            first = last = _intronic_to_genomic(variant.StartLoc.IVSLoc,
                                                transcript)
        except _InvalidIntronError as e:
            output.addMessage(__file__, 4, 'EINVALIDINTRON',
                              'Non-existing intron number %d given.' % \
                              e.intron)
            raise
        if variant.EndLoc and variant.EndLoc.IVSLoc:
            try:
                last = _intronic_to_genomic(variant.EndLoc.IVSLoc, transcript)
            except _InvalidIntronError as e:
                output.addMessage(__file__, 4, 'EINVALIDINTRON',
                                  'Non-existing intron number %d given.' % \
                                  e.intron)
                raise
            if last < first:
                # Todo: Why could this ever happen?
                first, last = last, first

    else:
        # Genomic or coding positioning.
        if record.record.molType != 'g' and \
               (_is_coding_intronic(variant.StartLoc) or
                _is_coding_intronic(variant.EndLoc)):
            output.addMessage(__file__, 3, 'ENOINTRON', 'Intronic ' \
                'position given for a non-genomic reference sequence.')
            raise _RawVariantError()

        first_location = last_location = variant.StartLoc.PtLoc
        if variant.EndLoc:
            last_location = variant.EndLoc.PtLoc

        # Todo: Check these error messages for minus strand variants etc.
        try:
            if transcript:
                # Coding positioning.
                first, last = _coding_to_genomic(first_location, last_location,
                                                 transcript, output)
                output.addOutput('rawVariantsCoding',
                                 (original_description, first_location, last_location))
            else:
                # Genomic positioning.
                first, last = _genomic_to_genomic(first_location, last_location)
        except _UnknownPositionError:
            output.addMessage(__file__, 4, 'EUNKNOWN',
                              'Unknown positions (denoted by "?") are ' \
                              'not supported.')
            raise
        except _OffsetSignError as e:
            output.addMessage(__file__, 3, 'EOFFSETSIGN', 'Offset %s from ' \
                              'position %s is in %s direction but should ' \
                              'be in %s direction.' % \
                              (e.offset, e.main,
                               'downstream' if e.acceptor else 'upstream',
                               'upstream' if e.acceptor else 'downstream'))
            raise
        except _OffsetNotFromBoundaryError as e:
            output.addMessage(__file__, 3, 'EOFFSETFROMBOUNDARY',
                              'Offset may not be from position %s because ' \
                              ' this is not an exon boundary.' % e.main)
            raise
        except _RawVariantError:
            # Todo: What does this situation really mean? I don't think
            # this is the right message.
            #output.addMessage(__file__, 3, 'ESPLICE', 'Invalid intronic ' \
            #                  'position given.')
            output.addMessage(__file__, 4, 'EPOSITION',
                              'Could not determine position.')
            raise

    if last < first:
        output.addMessage(__file__, 3, 'ERANGE', 'End position is smaller than ' \
                          'the begin position.')
        raise _RawVariantError()

    if first < 1:
        output.addMessage(__file__, 4, 'ERANGE', 'Position %i is out of range.' %
                          first)
        raise _RawVariantError()

    if last > len(mutator.orig):
        output.addMessage(__file__, 4, 'ERANGE', 'Position %s is out of range.' %
                          last)
        raise _RawVariantError()

    splice_abort = False

    # If we hit a splice site, issue a warning. Later on we decide if we
    # can still process this variant in any way (e.g. if it deletes an
    # entire exon).
    if transcript and util.over_splice_site(first, last, transcript.CM.RNA):
        splice_abort = True
        output.addMessage(__file__, 2, 'WOVERSPLICE',
                          'Variant hits one or more splice sites in ' \
                          'selected transcript.')

    # If we have a deletion, and it covers exactly an even number of splice
    # sites, remove these splice sites.
    # Note, this is not the same as util.over_splice_site(). Here we collect
    # sites where the deletion borders the exon/intron boundary.
    # Todo: Special cases for first/last exon? Upstream/downstream exons?
    # Todo: This still goes horribly wrong in some cases, example:
    # NM_000088.3(COL1A1_v001):c.588del
    if transcript and variant.MutationType == 'del':
        removed_sites = []
        for acceptor, donor in util.grouper(transcript.CM.RNA):

            # If we have introns, we match splice sites in a fuzzy way. This
            # Means that in the case of
            #
            #               a            b
            #     ===========------------=============
            #
            # with splice sites a and b, a deletion a+1_b-1 of the entire
            # intron gets treated as a deletion of both splice sites.
            #
            # We don't want this behaviour on e.g. RNA, where we only have
            # exons. In the case of
            #
            #              a b           c d
            #     ========== ============= ===========
            #
            # with splice sites a b c d, a deletion b_c of the middle exon
            # should only remove splice sites b and c, not a and d.
            if record.record.molType == 'g':
                fuzzy = 1
            else:
                fuzzy = 0

            if first <= acceptor <= last + fuzzy:
                removed_sites.append(acceptor)
            if first - fuzzy <= donor <= last:
                removed_sites.append(donor)

        if len(removed_sites) and not len(removed_sites) % 2:
            # An even number of splice sites was removed. We can deal with
            # this, but issue a warning.
            # However, don't do this trick if we end up removing an odd number
            # of sites from the CDS.
            # Todo: We might cripple the start codon, fix the translation code
            # (further down) to deal with this.
            # Todo: Bit unrelated, but find out the difference between
            # - transcript.CM.RNA
            # - transcript.mRNA.positionList
            # and what we should use (likewise for CDS).
            removed_cds_sites = filter(lambda s: s in transcript.CDS.positionList,
                                       removed_sites) if transcript.CDS else []
            if not len(removed_cds_sites) % 2:
                # Todo: splice_abort=False undoes the warning (sort of), but
                # the warning might (also) be about other sites...
                splice_abort = False
                mutator.add_removed_sites(removed_sites)
                output.addMessage(__file__, 1, 'IDELSPLICE',
                                  'Removed %i splice sites from selected ' \
                                  'transcript.' % len(removed_sites))
                # This is primarily for use in unittests.
                output.addOutput('removedSpliceSites', len(removed_sites))

    # If splice_abort is set, this basically means WOVERSPLICE was called and
    # IDELSPLICE was not called.
    # I guess in that case we do want to generate the visualisation, the
    # genomic description, and affected transcripts. But NOT the predicted
    # protein.
    # The following solution is a bit of a hack. By setting the .translate
    # field of the transcript to False, we force that no protein is predicted.
    if splice_abort:
        transcript.translate = False

    # The following functions can raise _RawVariantError exceptions, but we
    # just let them flow through to the caller.

    # Check if the (optional) argument is valid.
    if variant.MutationType in ['del', 'dup', 'subst', 'delins']:
        _check_argument(argument, mutator.orig, first, last, output)

    # Substitution.
    if variant.MutationType == 'subst':
        apply_substitution(first, argument, sequence, mutator, record, output)

    # Deletion or duplication.
    if variant.MutationType in ['del', 'dup']:
        # The fuzzy flags are to support deletions of the form c.a-?_b+?del.
        first_fuzzy = variant.StartLoc and variant.StartLoc.PtLoc.Offset == '?'
        last_fuzzy = variant.EndLoc and variant.EndLoc.PtLoc.Offset == '?'
        apply_deletion_duplication(first, last, variant.MutationType, mutator,
                                   record, output, first_fuzzy=first_fuzzy,
                                   last_fuzzy=last_fuzzy)

    # Inversion.
    if variant.MutationType == 'inv':
        apply_inversion(first, last, mutator, record, output)

    # Insertion.
    if variant.MutationType == 'ins':
        # Check if the inserted sequence is not a range.
        # Todo: Implement this feature.
        if not argument:
            output.addMessage(__file__, 4, 'ENOTIMPLEMENTED',
                              'Insertion of a range is not implemented yet.')
            raise _RangeInsertionError()
        apply_insertion(first, last, argument, mutator, record, output)

    # DelIns.
    if variant.MutationType == 'delins':
        # Check if the inserted sequence is not a range.
        # Todo: Implement this feature.
        if not sequence:
            output.addMessage(__file__, 4, 'ENOTIMPLEMENTED',
                              'Insertion of a range is not implemented yet.')
            raise _RangeInsertionError()
        apply_delins(first, last, argument, sequence, mutator, record, output)
#process_raw_variant


def _add_static_transcript_info(transcript, output):
    """
    Add static (unrelated to the variant) transcript-specific information to
    the {output} object.

    @arg transcript: A transcript object.
    @type transcript: Modules.GenRecord.Locus
    @arg output: The Output object.
    @type output: Modules.Output.Output
    """
    output.addOutput('hasTranscriptInfo', True)

    # Add exon table to output.
    for i in range(0, transcript.CM.numberOfExons() * 2, 2):
        acceptor = transcript.CM.getSpliceSite(i)
        donor = transcript.CM.getSpliceSite(i + 1)
        output.addOutput('exonInfo', [acceptor, donor,
                                      transcript.CM.g2c(acceptor),
                                      transcript.CM.g2c(donor)])

    # Add CDS info to output.
    cds_stop = transcript.CM.info()[2]
    output.addOutput('cdsStart_g', transcript.CM.x2g(1, 0))
    output.addOutput('cdsStart_c', 1)
    output.addOutput('cdsStop_g', transcript.CM.x2g(cds_stop, 0))
    output.addOutput('cdsStop_c', cds_stop)

    # Is this transcript coding?
    # Example of non-coding transcript variant:
    # AL449423.14(CDKN2A_v004):n.42_437del
    output.addOutput('transcriptCoding', bool(transcript.CM.CDS))

    # Is this transcript on the reverse strand?
    output.addOutput('transcriptReverse', transcript.CM.orientation == -1)


def _add_transcript_info(mutator, transcript, output):
    """
    Add transcript-specific information (including protein prediction) to
    the {output} object.

    @arg mutator: A Mutator instance.
    @type mutator: mutalyzer.mutator.Mutator
    @arg transcript: A transcript object.
    @type transcript: Modules.GenRecord.Locus
    @arg output: The Output object.
    @type output: Modules.Output.Output

    @todo: Documentation.
    @todo: Don't generate the fancy HTML protein descriptions here.
    @todo: Add mutated transcript and CDS info.
    """
    # Add transcript info to output.
    if transcript.transcribe:
        output.addOutput('myTranscriptDescription', transcript.description or '=')
        output.addOutput('origMRNA',
            str(util.splice(mutator.orig, transcript.mRNA.positionList)))
        output.addOutput('mutatedMRNA',
            str(util.splice(mutator.mutated,
                        mutator.newSplice(transcript.mRNA.positionList))))

    # Add protein prediction to output.
    if transcript.translate:
        cds_original = Seq(str(util.splice(mutator.orig, transcript.CDS.positionList)),
                           IUPAC.unambiguous_dna)
        cds_variant = Seq(str(util.__nsplice(mutator.mutated,
                                        mutator.newSplice(transcript.mRNA.positionList),
                                        mutator.newSplice(transcript.CDS.location),
                                        transcript.CM.orientation)),
                          IUPAC.unambiguous_dna)

        #output.addOutput('origCDS', cds_original)

        if transcript.CM.orientation == -1:
            cds_original = Bio.Seq.reverse_complement(cds_original)
            cds_variant = Bio.Seq.reverse_complement(cds_variant)

        if not util.is_dna(cds_original):
            output.addMessage(__file__, 4, 'ENODNA',
                              'Invalid letters in reference sequence.')
            return

        if '*' in cds_original.translate(table=transcript.txTable)[:-1]:
            output.addMessage(__file__, 3, 'ESTOP',
                              'In frame stop codon found.')
            return

        protein_original = cds_original.translate(table=transcript.txTable,
                                                  to_stop=True)
        protein_variant = cds_variant.translate(table=transcript.txTable,
                                                to_stop=True)

        # Note: addOutput('origCDS', ...) was first before the possible
        #       reverse complement operation above.
        output.addOutput('origCDS', cds_original)
        output.addOutput("newCDS", cds_variant[:(len(str(protein_variant)) + 1) * 3])

        output.addOutput('oldprotein', protein_original + '*')

        # Todo: Don't generate the fancy HTML protein views here, do this in
        # website.py.
        # I think it would also be nice to include the mutated list of splice
        # sites.
        if not protein_variant or protein_variant[0] != 'M':
            # Todo: Protein differences are not color-coded,
            # use something like below in protein_description().
            util.print_protein_html(protein_original + '*', 0, 0, output,
                                    'oldProteinFancy')
            if str(cds_variant[0:3]) in \
                   Bio.Data.CodonTable.unambiguous_dna_by_id \
                   [transcript.txTable].start_codons:
                output.addOutput('newprotein', '?')
                util.print_protein_html('?', 0, 0, output, 'newProteinFancy')
                output.addOutput('altStart', str(cds_variant[0:3]))
                if str(protein_original[1:]) != str(protein_variant[1:]):
                    output.addOutput('altProtein',
                                     'M' + protein_variant[1:] + '*')
                    util.print_protein_html('M' + protein_variant[1:] + '*', 0, 0,
                                            output, 'altProteinFancy')
            else :
                output.addOutput('newprotein', '?')
                util.print_protein_html('?', 0, 0, output, 'newProteinFancy')

        else:
            cds_length = util.cds_length(
                mutator.newSplice(transcript.CDS.positionList))
            descr, first, last_original, last_variant = \
                   util.protein_description(cds_length, protein_original,
                                            protein_variant)

            # This is never used.
            output.addOutput('myProteinDescription', descr)

            util.print_protein_html(protein_original + '*', first, last_original,
                                    output, 'oldProteinFancy')
            if str(protein_original) != str(protein_variant):
                output.addOutput('newprotein', protein_variant + '*')
                util.print_protein_html(protein_variant + '*', first, last_variant,
                                        output, 'newProteinFancy')
#_add_transcript_info


def process_variant(mutator, description, record, output):
    """
    @arg mutator: A Mutator instance.
    @type mutator: mutalyzer.mutator.Mutator
    @arg description: Parsed HGVS variant description.
    @type description: pyparsing.ParseResults
    @arg record: A GenRecord object.
    @type record: Modules.GenRecord.GenRecord
    @arg output: The Output object.
    @type output: Modules.Output.Output

    @raise _VariantError: Cannot process this variant.

    @todo: Documentation.
    """
    if not description.RawVar and not description.SingleAlleleVarSet:
        output.addMessage(__file__, 4, 'ENOVARIANT',
                          'Variant description contains no mutation.')
        raise _VariantError()

    if description.RefType == 'r':
        output.addMessage(__file__, 4, 'ERNA',
                          'Descriptions on RNA level are not supported.')
        raise _VariantError()

    transcript = None

    if description.RefType in ['c', 'n']:

        gene = None
        gene_symbol, transcript_id = output.getOutput('geneSymbol')[-1]

        if description.LrgAcc:
            # LRG case, pick the top gene.
            gene = record.record.geneList[0]
            if transcript_id:
                transcript = gene.findLocus(transcript_id)
                if not transcript:
                    # Todo: Incorrect error message, it might also be that
                    # there are no transcripts at all (e.g. N4BP2L1 on
                    # NG_012772.1).
                    output.addMessage(__file__, 4, "ENOTRANSCRIPT",
                        "Multiple transcripts found for gene %s. Please " \
                        "choose from: %s" %(gene.name,
                            ", ".join(gene.listLoci())))
            else:
                # No transcript id given.
                if len(gene.transcriptList) == 1:
                    # No transcript given, only 1 found, pick that.
                    transcript = gene.transcriptList[0]
                else:
                    output.addMessage(__file__, 4, "ENOTRANSCRIPT",
                        "No transcript given for gene %s. Please " \
                        "choose from: %s" %(gene.name,
                            ", ".join(gene.listLoci())))

        else:
            # Not an LRG, find our gene manually.
            genes = record.record.listGenes()
            transcript_id = transcript_id and "%.3i" % int(transcript_id)

            if gene_symbol in genes:
                # We found our gene.
                gene = record.record.findGene(gene_symbol)
            elif (len(genes) == 1) and not(gene_symbol):
                # No gene given and there is only one gene in the record.
                # Todo: message?
                gene = record.record.geneList[0]
            else:
                output.addMessage(__file__, 4, "EINVALIDGENE",
                    "Gene %s not found. Please choose from: %s" % (
                    gene_symbol, ", ".join(genes)))

            if gene:
                # Find transcript.
                transcripts = gene.listLoci()
                if transcript_id in transcripts:
                    # Found our transcript.
                    transcript = gene.findLocus(transcript_id)
                elif (len(transcripts) == 1) and not(transcript_id):
                    # No transcript given and there is only one transcript for
                    # this gene.
                    transcript = gene.transcriptList[0]
                else:
                    # Todo: Incorrect error message, it might also be that
                    # there are no transcripts at all (e.g. N4BP2L1 on
                    # NG_012772.1).
                    output.addMessage(__file__, 4, "ENOTRANSCRIPT",
                        "Multiple transcripts found for gene %s. Please " \
                        "choose from: %s" %(gene.name,
                        ", ".join(gene.listLoci())))

        # Add selected gene symbol to output
        output.addOutput('geneSymbol', (gene and gene.name or '',
                                        transcript and transcript.name or ''))

        # Return if no transcript is selected
        if not transcript:
            # Skip all BatchJobs with the same preColon data.
            output.addOutput('BatchFlags',
                             ('S2', output.getOutput('preColon')[-1]))
            raise _VariantError()
        elif not transcript.transcribe:
            # Todo: Shouldn't we add some message here?
            raise _VariantError()

        # Mark this as the current transcript we work with.
        transcript.current = True

    # Add static transcript-specific information.
    if transcript and record.record.geneList:
        _add_static_transcript_info(transcript, output)

    # Now process all raw variants (or just the only one). The function
    # process_raw_variant might raise a _VariantError exception.
    if description.SingleAlleleVarSet:
        for var in description.SingleAlleleVarSet:
            try:
                process_raw_variant(mutator, var, record, transcript,
                                    output)
            except _RawVariantError:
                #output.addMessage(__file__, 2, 'WSKIPRAWVARIANT',
                #                  'Ignoring raw variant "%s".' % var[0])
                output.addMessage(__file__, 1, 'IRAWVARIANTERROR',
                                  'Aborted variant check due to error in ' \
                                  'raw variant "%s".' % var[-1])
                raise
    else:
        process_raw_variant(mutator, description, record,
                            transcript, output)

    # Add transcript-specific variant information.
    if transcript and record.record.geneList:
        _add_transcript_info(mutator, transcript, output)
#process_variant


def check_variant(description, output):
    """
    Check the variant described by {description} according to the HGVS variant
    nomenclature and populate the {output} object with various information
    about the variant and its reference sequence.

    @arg description: Variant description in HGVS notation.
    @type description: string
    @arg output: An output object.
    @type output: Modules.Output.Output

    @todo: Documentation.
    @todo: Raise exceptions on failure instead of just return.
    """
    output.addOutput('inputvariant', description)

    grammar = Grammar(output)
    parsed_description = grammar.parse(description)

    if not parsed_description:
        # Parsing went wrong.
        return

    # Add GeneSymbol and Transcript Var to the Output object for batch.
    if parsed_description.Gene:
        output.addOutput('geneOfInterest',
                         dict(parsed_description.Gene.items()))
    else:
        output.addOutput('geneOfInterest', dict())

    if parsed_description.LrgAcc:
        record_id = parsed_description.LrgAcc
    elif parsed_description.Version:
        record_id = parsed_description.RefSeqAcc + '.' + parsed_description.Version
    else:
        record_id = parsed_description.RefSeqAcc

    if not record_id:
        output.addMessage(__file__, 4, 'ENOREF', 'No reference sequence given.')
        return

    gene_symbol = transcript_id = ''

    database = Db.Cache()
    if parsed_description.LrgAcc:
        filetype = 'LRG'
        transcript_id = parsed_description.LRGTranscriptID
        retriever = Retriever.LRGRetriever(output, database)
    else:
        filetype = 'GB'
        if parsed_description.Gene:
            gene_symbol = parsed_description.Gene.GeneSymbol or ''
            transcript_id = parsed_description.Gene.TransVar or ''
            if parsed_description.Gene.ProtIso:
                output.addMessage(__file__, 4, 'EPROT',
                    'Indexing by protein isoform is not supported.')
        retriever = Retriever.GenBankRetriever(output, database)

    retrieved_record = retriever.loadrecord(record_id)

    if not retrieved_record:
        return

    # Add recordType to output for output formatting.
    output.addOutput('recordType', filetype)
    output.addOutput('organism', retrieved_record.organism)
    output.addOutput('reference', record_id)

    # Add some more reference info.
    # Todo: Add selected transcript (LRGTranscriptID for LRGs, Gene for
    #     genbank files), but I think this is part of the variant and not
    #     part of the reference.
    output.addOutput('reference_id', retrieved_record.id)
    output.addOutput('source_id', retrieved_record.source_id)
    if filetype == 'GB':
        output.addOutput('source_accession', retrieved_record.source_accession)
        output.addOutput('source_version', retrieved_record.source_version)
        output.addOutput('source_gi', retrieved_record.source_gi)
    output.addOutput('molecule', retrieved_record.molType)

    # Note: geneSymbol[0] is used as a filter for batch runs.
    output.addOutput('geneSymbol', (gene_symbol, transcript_id))

    # Note: preColon is used to filter out Batch entries that will result in
    # identical errors.
    output.addOutput('preColon', description.split(':')[0])
    output.addOutput('variant', description.split(':')[-1])

    record = GenRecord.GenRecord(output)
    record.record = retrieved_record
    record.checkRecord()

    # Create the legend.
    for gene in record.record.geneList:
        for transcript in sorted(gene.transcriptList, key=attrgetter('name')):
            if not transcript.name:
                continue
            output.addOutput('legends',
                             ['%s_v%s' % (gene.name, transcript.name),
                              transcript.transcriptID, transcript.locusTag,
                              transcript.transcriptProduct,
                              transcript.linkMethod])
            if transcript.translate:
                output.addOutput('legends',
                                 ['%s_i%s' % (gene.name, transcript.name),
                                  transcript.proteinID, transcript.locusTag,
                                  transcript.proteinProduct,
                                  transcript.linkMethod])

    # Note: The GenRecord instance is carrying the sequence in .record.seq.
    #       So is the Mutator instance in .mutator.orig.

    mutator = Mutator(record.record.seq, output)

    # Todo: If processing of the variant fails, we might still want to show
    # information about the record, gene, transcript.

    try:
        process_variant(mutator, parsed_description, record, output)
    except _VariantError:
        return

    output.addOutput('original', str(mutator.orig))
    output.addOutput('mutated', str(mutator.mutated))

    # Chromosomal region (only for GenBank human transcript references).
    # This is still quite ugly code, and should be cleaned up once we have
    # a refactored mapping module.
    if retrieved_record.organism == 'Homo sapiens' \
           and parsed_description.RefSeqAcc \
           and parsed_description.RefType == 'c':
        raw_variants = output.getOutput('rawVariantsCoding')
        # Example value: [('29+4T>C', 29+4, 29+4), ('230_233del', 230, 233)]
        if raw_variants:
            locations = [pos
                         for descr, first, last in raw_variants
                         for pos in (first, last)]
            converter = Converter('hg19', output)
            chromosomal_positions = converter.chromosomal_positions(
                locations, parsed_description.RefSeqAcc, parsed_description.Version)
            if chromosomal_positions:
                output.addOutput('rawVariantsChromosomal',
                                 (chromosomal_positions[0], chromosomal_positions[1],
                                  zip([descr for descr, first, last in raw_variants],
                                      util.grouper(chromosomal_positions[2]))))
                # Example value: ('chr12', [('29+4T>C', (2323, 2323)), ('230_233del', (5342, 5345))])

    # Protein.
    for gene in record.record.geneList:
        for transcript in gene.transcriptList:

            if not (transcript.CDS and transcript.translate) \
                   or ';' in transcript.description \
                   or transcript.description == '?':
                # Default value is '?', but later on we don't prefix a 'p.'
                # string, so we include it here. If there's no good reason
                # for this, I think we should only add the 'p.' later (so
                # __toProtDescr should also not add it).
                transcript.proteinDescription = 'p.?'
                continue

            cds_original = Seq(str(util.splice(mutator.orig, transcript.CDS.positionList)),
                               IUPAC.unambiguous_dna)
            cds_variant = Seq(str(util.__nsplice(mutator.mutated,
                                            mutator.newSplice(transcript.mRNA.positionList),
                                            mutator.newSplice(transcript.CDS.location),
                                            transcript.CM.orientation)),
                              IUPAC.unambiguous_dna)

            if transcript.CM.orientation == -1:
                cds_original = Bio.Seq.reverse_complement(cds_original)
                cds_variant = Bio.Seq.reverse_complement(cds_variant)

            #if '*' in cds_original.translate()[:-1]:
            #    output.addMessage(__file__, 3, "ESTOP",
            #                      "In frame stop codon found.")
            #    return
            ##if

            # Todo: Figure out if this is all ok, even if the CDS stop is
            # somehow removed, if the sequence is really short, etc.

            if not len(cds_original) % 3:
                try:
                    # FIXME this is a bit of a rancid fix.
                    protein_original = cds_original.translate(table=transcript.txTable,
                                                              cds=True,
                                                              to_stop=True)
                except Bio.Data.CodonTable.TranslationError:
                    output.addMessage(__file__, 4, "ETRANS", "Original " \
                                      "CDS could not be translated.")
                    return
                protein_variant = cds_variant.translate(table=transcript.txTable,
                                                        to_stop=True)
                try:
                    cds_length = util.cds_length(
                        mutator.newSplice(transcript.CDS.positionList))
                    transcript.proteinDescription = util.protein_description(
                        cds_length, protein_original, protein_variant)[0]
                except IndexError:
                    # Todo: Probably CDS start was hit by removal of exon...
                    transcript.proteinDescription = 'p.?'

            else:
                output.addMessage(__file__, 2, "ECDS", "CDS length is " \
                    "not a multiple of three in gene %s, transcript " \
                    "variant %s." % (gene.name, transcript.name))
                transcript.proteinDescription = 'p.?'

    reference = output.getOutput('reference')[-1]
    if ';' in record.record.description:
        generated_description = '[' + record.record.description + ']'
    else:
        generated_description = record.record.description or '='

    output.addOutput('genomicDescription', '%s:%c.%s' % \
                     (reference, record.record.molType, generated_description))
    output.addOutput('gDescription', '%c.%s' % \
                     (record.record.molType, generated_description))
    output.addOutput('molType', record.record.molType)

    if record.record.chromOffset:
        if ';' in record.record.chromDescription:
            chromosomal_description = '[' + record.record.chromDescription + ']'
        else:
            chromosomal_description = record.record.chromDescription or '='
        output.addOutput('genomicChromDescription', '%s:%c.%s' % \
                         (record.record.source_id,
                          record.record.molType, chromosomal_description))

    # Now we add variant descriptions for all transcripts, including protein
    # level descriptions.
    for gene in record.record.geneList:
        for transcript in sorted(gene.transcriptList, key=attrgetter('name')):

            # Note: I don't think genomic_id is ever used, because it is
            # always ''.
            coding_description = ''
            protein_description = ''
            full_description = ''
            full_protein_description = ''
            genomic_id = coding_id = protein_id = ''

            if ';' in transcript.description:
                generated_description = '[' + transcript.description + ']'
            else:
                generated_description = transcript.description or '='

            if record.record._sourcetype == 'LRG':
                if transcript.name:
                    full_description = '%st%s:%c.%s' % \
                                       (reference, transcript.name,
                                        transcript.molType,
                                        generated_description)
                    output.addOutput('descriptions', full_description)
                else:
                    output.addOutput('descriptions', gene.name)
            else:
                full_description = '%s(%s_v%s):%c.%s' % \
                                   (reference, gene.name, transcript.name,
                                    transcript.molType,
                                    generated_description)
                output.addOutput('descriptions', full_description)

            if transcript.molType == 'c':
                coding_description = 'c.%s' % generated_description
                protein_description = transcript.proteinDescription
                if record.record._sourcetype == 'LRG':
                    full_protein_description = '%sp%s:%s' % \
                                               (reference, transcript.name,
                                                protein_description)
                else:
                    full_protein_description = '%s(%s_i%s):%s' % \
                                               (reference, gene.name,
                                                transcript.name,
                                                protein_description)

                coding_id, protein_id = \
                           transcript.transcriptID, transcript.proteinID
                output.addOutput('protDescriptions',
                                 full_protein_description)

            # The 'NewDescriptions' field is used in _add_batch_output.
            output.addOutput('NewDescriptions',
                             (gene.name, transcript.name,
                              transcript.molType, coding_description,
                              protein_description, genomic_id, coding_id,
                              protein_id, full_description,
                              full_protein_description))

    _add_batch_output(output)
#check_variant
