"""
General utility functions.

@todo: All these functions come from the old Mutalyzer.py file. Try to find
       general utility functions in other modules too.
@todo: Use exceptions for failure handling.
@todo: End vs stop. I guess we should use start/stop (end goes with beginning).
       Or first/last, or acceptor/donor. Anyway, CDS is always denoted with
       start/stop. Important thing is that the semantics should be clear.
       Idea:
       * CDS -> use start/stop
       * splice sites or exons -> acceptor/donor
       * translation -> begin/end
       * any range of bases -> first/last
       * interbase position (if two numbers are used) -> before/after
@todo: We can also group this in separate files in a util/ directory, according
       to function (e.g. util/sequences.py, util/positioning.py, etc).
@todo: Unit tests (some can directly be extracted from the docstring).
"""


import os
import math
import time
from itertools import izip_longest

import Bio.Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import seq3


def grouper(iterable, n=2, fillvalue=None):
    """
    Make an iterator that takes {n} elements at a time from {iterable}, using
    {fillvalue} as default values where we don't have a multiple of {n}.

        >>> for g in grouper('ABCDEFG', 3, 'x'):
        ...     print g
        ('A', 'B', 'C')
        ('D', 'E', 'F')
        ('G', 'x', 'x')

        >>> splice_sites = [1, 4, 5, 12, 13, 18]
        >>> for acceptor, donor in grouper(splice_sites):
        ...     print 'Exon of length %d' % (donor - acceptor + 1)
        Exon of length 4
        Exon of length 8
        Exon of length 6

    Modified from the example at [1].

    @arg iterable: Iterable to take groups of elements from.
    @type iterable: any iterable type
    @kwarg n: Number of elements to take at a time (default 2).
    @type n: int
    @kwarg fillvalue: Default value to use as padding if length of {iterable}
                      is not a multiple of {n} (default None).
    @return: Iterator that gives elements of {iterable} as groups of {n}.
    @rtype: tuple

    [1] http://docs.python.org/library/itertools.html#recipes
    """
    args = [iter(iterable)] * n
    return izip_longest(*args, fillvalue=fillvalue)
#grouper


def over_splice_site(first, last, splice_sites):
    """
    Check wheter a genomic range {first}_{last} hits a splice site. Hitting
    a splice site means that the range covers the exon boundary.

        >>> splice_sites = [1, 4, 8, 12, 19, 28]
        >>> over_splice_site(7, 8, splice_sites)
        True
        >>> over_splice_site(12, 13, splice_sites)
        True
        >>> over_splice_site(8, 9, splice_sites)
        False
        >>> over_splice_site(8, 8, splice_sites)
        False

    @arg first: The first coordinate of the range in g. notation.
    @type first: int
    @arg last: The last coordinate of the range in g. notation.
    @type last: int
    @arg sites: A list of splice sites in g. notation.
    @type sites: list(int)

    @return: True if one or more splice sites are hit, False otherwise.
    @rtype: boolean

    @todo: Assert number of splice sites is even.
    """
    for acceptor, donor in grouper(splice_sites):
        if first < acceptor and last >= acceptor:
            return True
        if donor and first <= donor and last > donor:
            return True

    return False
#over_splice_site


def splice(s, splice_sites):
    """
    Construct the transcript or the coding sequence from a record and a list
    of splice sites.

        >>> splice('abcdefghijklmnopqrstuvwxyz', [2, 4, 7, 16, 20, 23])
        'bcdghijklmnoptuvw'

    @arg s: A DNA sequence.
    @type s: string
    @arg splice_sites: A list of even length of integers.
    @type splice_sites: list

    @return: The concatenation of slices from the sequence that is present
             in the GenBank record.
    @rtype: string

    @todo: Assert length of splice_sites is even.
    """
    transcript = ''

    for acceptor, donor in grouper(splice_sites):
        transcript += s[acceptor - 1:donor]

    return transcript
#splice


# Todo: refactor
def __nsplice(string, splice_sites, CDS, orientation) :
    """
    Just like _splice(), but it only keeps the parts between CDS[0] and
    CDS[1] (in the right orientation).

    I guess we could easily do this as a separate step after _splice()?

    @todo: keep this function?
    @todo: documentation
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
                transcript += string[splice_sites[i] - 1:CDS[1]]
            else :
                if splice_sites[i] < CDS[1] :
                    transcript += \
                        string[splice_sites[i] - 1:splice_sites[i + 1]]
        #for
    #else

    return transcript
#__nsplice


def cds_length(splice_sites):
    """
    Calculate the length of a CDS.

        >>> cds_length([2, 4, 7, 16, 20, 23])
        17

    @arg splice_sites: The coordinates of the CDS including internal splice
                       sites.
    @type splice_sites: list

    @return: Length of the CDS.
    @rtype: int

    @todo: Assert length of splice_sites is even.
    """
    l = 0

    for acceptor, donor in grouper(splice_sites):
        l += donor - acceptor + 1

    return l
#cds_length


def format_range(first, last):
    """
    Simplify a range to one position when applicable.

        >>> format_range(3, 5)
        '3_5'
        >>> format_range(3, 3)
        '3'

    @arg first: First coordinate of a range.
    @type first: integer
    @arg last: Second coordinate of a range.
    @type last: integer

    @return: {first}_{last} in case of a real range, {first} otherwise.
    @rtype: string
    """
    if first == last:
        return str(first)

    return '%i_%i' % (first, last)
#format_range


def roll(s, first, last):
    """
    Determine the variability of a variant by looking at cyclic
    permutations. Not all cyclic permutations are tested at each time, it
    is sufficient to check ``aW'' if ``Wa'' matches (with ``a'' a letter,
    ``W'' a word) when rolling to the left for example.

        >>> roll('abbabbabbabb', 4, 6)
        (3, 6)
        >>> roll('abbabbabbabb', 5, 5)
        (0, 1)
        >>> roll('abcccccde', 4, 4)
        (1, 3)

    @arg s: A reference sequence.
    @type s: string
    @arg first: First position of the pattern in the reference sequence.
    @type first: int
    @arg last: Last position of the pattern in the reference sequence.
    @type last: int

    @return: tuple:
        - left  ; Amount of positions that the pattern can be shifted to
                  the left.
        - right ; Amount of positions that the pattern can be shifted to
                  the right.
    @rtype: tuple(int, int)
    """
    pattern = s[first - 1:last]   # Extract the pattern
    pattern_length = len(pattern)

    # Keep rolling to the left as long as a cyclic permutation matches.
    minimum = first - 2
    j = pattern_length - 1
    while minimum > -1 and s[minimum] == pattern[j % pattern_length]:
        j -= 1
        minimum -= 1

    # Keep rolling to the right as long as a cyclic permutation matches.
    maximum = last
    j = 0
    while maximum < len(s) and s[maximum] == pattern[j % pattern_length]:
        j += 1
        maximum += 1

    return first - minimum - 2, maximum - last
#roll


def palinsnoop(s):
    """
    Check a sequence for a reverse-complement-palindromic prefix (and
    suffix). If one is detected, return the length of this prefix. If the
    string equals its reverse complement, return -1.

        >>> palinsnoop('TACGCTA')
        2
        >>> palinsnoop('TACGTA')
        -1
        >>> palinsnoop('TACGCTT')
        0

    @arg s: A nucleotide sequence.
    @type s: string

    @return: The number of elements that are palindromic or -1 if the string
             is a 'palindrome'.
    @rtype: string
    """
    s_revcomp = Bio.Seq.reverse_complement(s)

    for i in range(int(math.ceil(len(s) / 2.0))):
        if s[i] != s_revcomp[i]:
            # The first i elements are 'palindromic'.
            return i

    # Perfect 'palindrome'.
    return -1
#palinsnoop


def longest_common_prefix(s1, s2):
    """
    Calculate the longest common prefix of two strings.

        >>> longest_common_prefix('abcdefg', 'abcabcdefg')
        'abc'
        >>> longest_common_prefix('abcdefg', 'abcdefg')
        'abcdefg'

    @arg s1: The first string.
    @type s1: string
    @arg s2: The second string.
    @type s2: string

    @return: The longest common prefix of s1 and s2.
    @rtype: string

    @todo: This is mostly used just for the length of the returned string,
           and we could also return that directly.
    """
    pos = 0

    while pos < min(len(s1), len(s2)) and s1[pos] == s2[pos]:
        pos += 1

    return s1[:pos]
#longest_common_prefix


def longest_common_suffix(s1, s2):
    """
    Calculate the longest common suffix of two strings.

        >>> longest_common_suffix('abcdefg', 'abcabcdefg')
        'abcdefg'
        >>> longest_common_suffix('abcdefg', 'abcefg')
        'efg'

    @arg s1: The first string.
    @type s1: string
    @arg s2: The second string.
    @type s2: string

    @return: The longest common suffix of s1 and s2.
    @rtype: string
    """
    return longest_common_prefix(s1[::-1], s2[::-1])[::-1]
#longest_common_suffix


def trim_common(s1, s2):
    """
    Given two strings, trim their longest common prefix and suffix.

        >>> trim_common('abcdefg', 'abcabcdefg')
        ('', 'abc', 3, 4)
        >>> trim_common('abcxyzefg', 'abcabcg')
        ('xyzef', 'abc', 3, 1)

    @arg s1: A string.
    @type s1: string
    @arg s2: Another string.
    @type s2: string

    @return: A tuple of:
        - string: Trimmed version of s1.
        - string: Trimmed version of s2.
        - int:    Length of longest common prefix.
        - int:    Length of longest common suffix.

    @todo: More intelligently handle longest_common_prefix().
    """
    lcp = len(longest_common_prefix(s1, s2))
    lcs = len(longest_common_suffix(s1[lcp:], s2[lcp:]))
    return s1[lcp:len(s1) - lcs], s2[lcp:len(s2) - lcs], lcp, lcs
#trim_common


def is_dna(s):
    """
    Check whether a string is a DNA string.

        >>> is_dna('TACTGT')
        True
        >>> is_dna('TACUGT')
        False

    @arg s: Any string or Bio.Seq.Seq instance.
    @type s: string

    @return: True if the string is a DNA string, False otherwise.
    @rtype: boolean
    """
    for i in str(s):
        if not i in IUPAC.unambiguous_dna.letters:
            return False

    return True
#is_dna


def in_frame_description(s1, s2) :
    """
    Give a description of an inframe difference of two proteins. Also give
    the position at which the proteins start to differ and the positions at
    which they are the same again.

        >>> in_frame_description('MTAPQQMT', 'MTAQQMT')
        ('p.(Pro4del)', 3, 4, 3)
        >>> in_frame_description('MTAPQQMT', 'MTAQMT')
        ('p.(Pro4_Gln5del)', 3, 5, 3)
        >>> in_frame_description('MTAPQQT', 'MTAQQMT')
        ('p.(Pro4_Gln6delinsGlnGlnMet)', 3, 6, 6)

    @arg s1: The original protein.
    @type s1: string
    @arg s2: The mutated protein.
    @type s2: string

    @return: A tuple of:
        - string ; Protein description of the change.
        - int    ; First position of the change.
        - int    ; Last position of the change in the first protein.
        - int    ; Last position of the change in the second protein.
    @rtype: tuple(string, int, int, int)

    @todo: More intelligently handle longest_common_prefix().
    """
    if s1 == s2:
        # Nothing happened.
        return ('p.(=)', 0, 0, 0)

    lcp = len(longest_common_prefix(s1, s2))
    lcs = len(longest_common_suffix(s1[lcp:], s2[lcp:]))
    s1_end = len(s1) - lcs
    s2_end = len(s2) - lcs

    # Insertion / Duplication / Extention.
    if not s1_end - lcp:
        if len(s1) == lcp:
            return ('p.(*%i%sext*%i)' % \
                    (len(s1) + 1, seq3(s2[len(s1)]), abs(len(s1) - len(s2))),
                    len(s1), len(s1), len(s2))
        ins_length = s2_end - lcp

        if lcp - ins_length >= 0 and s1[lcp - ins_length:lcp] == s2[lcp:s2_end]:
            if ins_length == 1:
                return ('p.(%s%idup)' % \
                        (seq3(s1[lcp - ins_length]), lcp - ins_length + 1),
                        lcp, lcp, lcp + 1)
            return ('p.(%s%i_%s%idup)' % \
                    (seq3(s1[lcp - ins_length]),
                     lcp - ins_length + 1, seq3(s1[lcp - 1]), lcp),
                    lcp, lcp, lcp + ins_length)
        #if
        return ('p.(%s%i_%s%iins%s)' % \
                (seq3(s1[lcp - 1]), lcp, seq3(s1[lcp]),
                 lcp + 1, seq3(s2[lcp:s2_end])),
                lcp, lcp, s2_end)
    #if

    # Deletion / Inframe stop.
    if not s2_end - lcp:
        if len(s2) == lcp:
            return ('p.(%s%i*)' % (seq3(s1[len(s2)]), len(s2) + 1),
                    0, 0, 0)

        if lcp + 1 == s1_end:
            return ('p.(%s%idel)' % (seq3(s1[lcp]), lcp + 1),
                    lcp, lcp + 1, lcp)
        return ('p.(%s%i_%s%idel)' % \
                (seq3(s1[lcp]), lcp + 1, seq3(s1[s1_end - 1]), s1_end),
                lcp, s1_end, lcp)
    #if

    # Substitution.
    if s1_end == s2_end and s1_end == lcp + 1:
        return ('p.(%s%i%s)' % (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp])),
                lcp, lcp + 1, lcp + 1)

    # InDel.
    if lcp + 1 == s1_end:
        return ('p.(%s%idelins%s)' % \
                (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp:s2_end])),
                lcp, lcp + 1, s2_end)
    return ('p.(%s%i_%s%idelins%s)' % \
            (seq3(s1[lcp]), lcp + 1, seq3(s1[s1_end - 1]), s1_end,
             seq3(s2[lcp:s2_end])),
            lcp, s1_end, s2_end)
#in_frame_description


def out_of_frame_description(s1, s2):
    """
    Give the description of an out of frame difference between two
    proteins. Give a description of an inframe difference of two proteins.
    Also give the position at which the proteins start to differ and the
    end positions (to be compatible with the in_frame_description function).

        >>> out_of_frame_description('MTAPQQMT', 'MTAQQMT')
        ('p.(Pro4Glnfs*5)', 3, 8, 7)
        >>> out_of_frame_description('MTAPQQMT', 'MTAQMT')
        ('p.(Pro4Glnfs*4)', 3, 8, 6)
        >>> out_of_frame_description('MTAPQQT', 'MTAQQMT')
        ('p.(Pro4Glnfs*5)', 3, 7, 7)

    @arg s1: The original protein.
    @type s1: string
    @arg s2: The mutated protein.
    @type s2: string

    @return: A tuple of:
        - string ; Protein description of the change.
        - int    ; First position of the change.
        - int    ; Last position of the first protein.
        - int    ; Last position of the second protein.
    @rtype: tuple(string, int, int, int)

    @todo: More intelligently handle longest_common_prefix().
    """
    lcp = len(longest_common_prefix(s1, s2))

    if lcp == len(s2): # NonSense mutation.
        if lcp == len(s1): # Is this correct?
            return ('p.(=)', 0, 0, 0)
        return ('p.(%s%i*)' % (seq3(s1[lcp]), lcp + 1), lcp, len(s1), lcp)
    if lcp == len(s1) :
        return ('p.(*%i%sext*%i)' % \
                (len(s1) + 1, seq3(s2[len(s1)]), abs(len(s1) - len(s2))),
                len(s1), len(s1), len(s2))
    return ('p.(%s%i%sfs*%i)' % \
            (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp]), len(s2) - lcp + 1),
            lcp, len(s1), len(s2))
#out_of_frame_description


def protein_description(cds_stop, s1, s2) :
    """
    Wrapper function for the in_frame_description() and
    out_of_frame_description() functions. It uses the value cds_stop to
    decide which one to call.

        >>> protein_description(34, 'MTAPQQMT', 'MTAQQMT')
        ('p.(Pro4Glnfs*5)', 3, 8, 7)
        >>> protein_description(33, 'MTAPQQMT', 'MTAQQMT')
        ('p.(Pro4del)', 3, 4, 3)
        >>> protein_description(33, 'MTAPQQMT', 'TTAQQMT')
        ('p.?', 0, 4, 3)

    @arg cds_stop: Position of the stop codon in c. notation (CDS length).
    @type cds_stop: int
    @arg s1: The original protein.
    @type s1: string
    @arg s2: The mutated protein.
    @type s2: string

    @return: A tuple of:
        - string ; Protein description of the change.
        - int    ; First position of the change.
        - int    ; Last position of the change in the first protein.
        - int    ; Last position of the change in the second protein.
    @rtype: tuple(string, int, int, int)
    """
    if cds_stop % 3:
        description = out_of_frame_description(str(s1), str(s2))
    else:
        description = in_frame_description(str(s1), str(s2))

    if not s2 or str(s1[0]) != str(s2[0]):
        # Mutation in start codon.
        return 'p.?', description[1], description[2], description[3]

    return description
#protein_description


# Todo: cleanup
def _insert_tag(s, pos1, pos2, tag1, tag2):
    """
    Insert two tags (tag1 and tag2) in string s at positions pos1 and pos2
    respectively if the positions are within the length of s. If not,
    either insert one tag or do nothing. If pos1 equals pos2, don't do
    anything either.

    @arg s: A sequence.
    @type s:
    @arg pos1: Position of tag1.
    @type pos1: int
    @arg pos2: Position of tag2.
    @type pos2: int
    @arg tag1: Content of tag1.
    @type tag1: string
    @arg tag2: Content of tag2.
    @type tag2: string

    @return: The original sequence, or a sequence with eiter tag1, tag2 or
             both tags inserted.
    @rtype: string

    @todo: Cleanup (note: only used in print_protein_html).
    """
    output = s
    block = len(s)

    # Only do something if pos1 != pos2.
    if pos1 != pos2:
        if 0 <= pos1 < block:
            # Insert tag1.
            output = output[:pos1] + tag1 + output[pos1:]
        if 0 <= pos2 < block:
            # Insert tag2.
            output = output[:-(block - pos2)] + tag2 \
                     + output[-(block - pos2):]

    return output
#_insert_tag


# Todo: cleanup
def print_protein_html(s, first, last, O, where):
    """
    Make a fancy representation of a protein and put it in the Output
    object under the name 'where'. The representation contains HTML tags
    and is suitable for viewing in a monospaced font.

    @arg s: A protein sequence.
    @type s: string
    @arg first: First position to highlight.
    @type first: int
    @arg last: Last position to highlight.
    @type last: int
    @arg O: The Output object.
    @type O: Modules.Output.Output
    @arg where: Location in the {O} object to store the representation.
    @type where: string

    @todo: Cleanup.
    """
    if not s: return

    block = 10        # Each block consists of 10 amino acids.
    line = 6 * block  # Each line consists of 6 blocks.

    tag1 = '<b style="color:#FF0000">'  # Use this tag for highlighting.
    tag2 = '</b>'                       # And this one to end highlighting.

    # The maximum length for positions is the 10_log of the length of the
    # protein.
    m = int(math.floor(math.log(len(s), 10)) + 1)
    o = 1

    # Add the first position.
    output = '%s ' % str(o).rjust(m)

    for i in range(0, len(s), block):
        # Add the blocks.
        output += ' ' + _insert_tag(s[i:i + block], first - i, last - i,
                                    tag1, tag2)
        if not (i + block) % line and i + block < len(s):
            # One line done.
            o += line
            O.addOutput(where, output)
            # Add the position (while escaping any potential highlighting).
            output = '<tt style="color:000000;font-weight:normal">%s</tt> ' \
                     % str(o).rjust(m)

    # Add last line.
    O.addOutput(where, output)
#print_protein_html


def generate_id():
    """
    Generates a (somewhat) unique number, using time().

    Note: Don't use this in very high frequencies, because it utilizes a
    short time.sleep() call to get a higher uniqueness.

    @return: A somewhat unique number.
    @rtype: int
    """
    unique_per_second = 100
    time.sleep(1.0 / unique_per_second)
    return int(time.time() * unique_per_second)
#generate_id


def nice_filename(filename):
    """
    Strip the path and the extention from a filename.

    @arg filename: A complete path plus extention.
    @type filename: string

    @return: The bare filename without a path and extention.
    @rtype: string
    """
    return filename.split('/')[-1].split('.')[0]
#nice_filename


def message_info(message):
    """
    Construct a dictionary with information about {message}.

    @arg message: A message instance.
    @type message: output.Message

    @return: A dictionary with fields of {message}.
    @rtype: dictionary
    """
    classes = {0: 'debug',
               1: 'information',
               2: 'warning',
               3: 'error',
               4: 'error'}

    return {'level':       message.named_level(),
            'origin':      message.origin,
            'class':       classes[message.level],
            'description': message.description}
#message_info


def slow(f):
    """
    Decorator for slow tests. This makes them to pass immediately, without
    running them. But only if the environment variable MUTALYZER_QUICK_TEST
    is 1.

    @todo: I don't think this actually belongs here (a separate util module
      for the unit tests?).
    """
    def slow_f(*args, **kwargs):
        if 'MUTALYZER_QUICK_TEST' in os.environ \
               and os.environ['MUTALYZER_QUICK_TEST'] == '1':
            return
        else:
            f(*args, **kwargs)
    return slow_f
#slow


def monkey_patch_suds():
    """
    Apply our monkey-patch for the suds package.

    For some weird reason the location http://www.w3.org/2001/xml.xsd is used
    for the XML namespace, but the W3C seems to respond too slow on that url.
    We therefore use http://www.w3.org/2009/01/xml.xsd which fixes this.

    Call this function before importing anything from the suds package.
    """
    from suds.xsd.sxbasic import Import
    _import_open = Import.open

    # Only apply the patch once.
    if getattr(Import, 'MUTALYZER_MONKEY_PATCHED', False):
        return

    def _import_open_patched(self, *args, **kwargs):
        if self.location == 'http://www.w3.org/2001/xml.xsd':
            self.location = 'http://www.w3.org/2009/01/xml.xsd'
        return _import_open(self, *args, **kwargs)

    Import.open = _import_open_patched
    Import.MUTALYZER_MONKEY_PATCHED = True
#monkey_patch_suds
