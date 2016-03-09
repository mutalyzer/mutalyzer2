"""
General utility functions.

@todo: Most of these functions come from the old Mutalyzer.py file. Try to
    find general utility functions in other modules too.
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


from __future__ import unicode_literals

from functools import wraps
import inspect
from itertools import izip_longest
import math
import operator
import sys
import time

from Bio import Seq
from Bio import SeqIO
from Bio.SeqUtils import seq3

# NOTE: This is a temporary fix.
from extractor.describe import palinsnoop, roll


def reverse_complement(sequence):
    """
    Reverse complement of a sequence represented as unicode string.

    Unfortunately, BioPython's reverse_complement doesn't work on unicode
    strings. We work almost exclusively with unicode strings, so this is a
    convenience wrapper.
    """
    return unicode(Seq.reverse_complement(str(sequence)))


def is_utf8_alias(encoding):
    """
    Returns `True` if the given encoding is recognized as UTF-8.
    """
    aliases = ('utf_8', 'u8', 'utf', 'utf8')
    return encoding.lower().replace('-', '_') in aliases


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
    @type s: any sequence type
    @arg splice_sites: A list of even length of integers.
    @type splice_sites: list

    @return: The concatenation of slices from the sequence that is present
             in the GenBank record.
    @rtype: type(s)

    @todo: Assert length of splice_sites is even.
    """
    transcript = s[:0]

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

    transcript = string[:0]
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
    @rtype: unicode
    """
    if first == last:
        return unicode(first)

    return '%i_%i' % (first, last)
#format_range


def roll_(s, start, end) :
    """
    Different (and easier) way of finding the variability of a substring.
    """
    # TODO: Start counting at 1, testing, etc.

    l = len(s)

    i = 0
    while end + i + 1 < l and s[start + i] == s[end + i + 1] :
        i += 1

    j = 0
    while start - j and s[start - j - 1] == s[end - j] :
        j += 1

    return j, i
#roll


def longest_common_prefix(s1, s2):
    """
    Calculate the longest common prefix of two strings.

        >>> longest_common_prefix('abcdefg', 'abcabcdefg')
        'abc'
        >>> longest_common_prefix('abcdefg', 'abcdefg')
        'abcdefg'

    @arg s1: The first string.
    @type s1: unicode
    @arg s2: The second string.
    @type s2: unicode

    @return: The longest common prefix of s1 and s2.
    @rtype: unicode

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
    @type s1: unicode
    @arg s2: The second string.
    @type s2: unicode

    @return: The longest common suffix of s1 and s2.
    @rtype: unicode
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
    @type s1: unicode
    @arg s2: Another string.
    @type s2: unicode

    @return: A tuple of:
        - unicode: Trimmed version of s1.
        - unicode: Trimmed version of s2.
        - int:     Length of longest common prefix.
        - int:     Length of longest common suffix.

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

    @arg s: Any string.
    @type s: unicode

    @return: True if the string is a DNA string, False otherwise.
    @rtype: boolean
    """
    for i in s:
        if i not in 'ATCG':
            return False

    return True
#is_dna


def guess_file_type(handle):
    """
    Guess the file type of an NGS data file.

    :arg file handle: Open readable handle to an NGS data file.

    :returns: Either 'fasta', 'fastq' or 'text'.
    :rtype: unicode
    """
    try:
        extension = getattr(handle, 'name').split('.')[-1]
    except AttributeError:
        pass
    else:
        if extension in ('fastq', 'fq'):
            return 'fastq'
        elif extension in ('fasta', 'fa'):
            return 'fasta'

    try:
        position = handle.tell()
        handle.seek(0)
    except IOError:
        return 'text'

    token = handle.read(1)
    handle.seek(position)

    if token == '>':
        return 'fasta'
    elif token == '@':
        return 'fastq'
    return 'text'


def read_dna(handle):
    """
    Read the first record in an NGS data file.

    If the format is not recognised as FASTA or FASTQ, we assume that the input
    is in plain text. In this case, DNA is converted to uppercase and all
    non-DNA characters are removed.

    :arg stream handle: Open readable handle to an NGS data file.

    :returns: Content of the first record in the file.
    :rtype: unicode
    """
    file_format = guess_file_type(handle)

    if file_format != 'text':
        return unicode(SeqIO.parse(handle, file_format).next().seq)

    return ''.join(x for x in unicode(handle.read()).upper() if x in 'ATCG')


def in_frame_description(s1, s2):
    """
    Give a description of an inframe difference of two proteins. Also give
    the position at which the proteins start to differ and the positions at
    which they are the same again.

        >>> in_frame_description('MTAPQQMT*', 'MTAQQMT*')
        ('p.(Pro4del)', 3, 4, 3)
        >>> in_frame_description('MTAPQQMT*', 'MTAQMT*')
        ('p.(Pro4_Gln5del)', 3, 5, 3)
        >>> in_frame_description('MTAPQQT*', 'MTAQQMT*')
        ('p.(Pro4_Gln6delinsGlnGlnMet)', 3, 6, 6)
        >>> in_frame_description('MTAPQQMT*', 'MTAPQQMTMQ*')
        ('p.(*9Metext*2)', 8, 9, 11)
        >>> in_frame_description('MTAPQQMT*', 'MTAPQQMTMQ')
        ('p.(*9Metext*?)', 8, 9, 10)

    @arg s1: The original protein.
    @type s1: unicode
    @arg s2: The mutated protein.
    @type s2: unicode

    @return: A tuple of:
        - unicode ; Protein description of the change.
        - int     ; First position of the change.
        - int     ; Last position of the change in the first protein.
        - int     ; Last position of the change in the second protein.
    @rtype: tuple(unicode, int, int, int)

    @todo: More intelligently handle longest_common_prefix().
    @todo: Refactor this code (too many return statements).
    """
    s2_stop = '*' in s2
    s1 = s1.rstrip('*')
    s2 = s2.rstrip('*')

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
            # http://www.hgvs.org/mutnomen/FAQ.html#nostop
            stop = unicode(abs(len(s1) - len(s2))) if s2_stop else '?'

            return ('p.(*%i%sext*%s)' % \
                    (len(s1) + 1, seq3(s2[len(s1)]), stop),
                    len(s1), len(s1) + 1, len(s2) + (1 if s2_stop else 0))

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
                    lcp, len(s1) + 1, len(s2) + 1)

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

        >>> out_of_frame_description('MTAPQQMT*', 'MTAQQMT*')
        ('p.(Pro4Glnfs*5)', 3, 9, 8)
        >>> out_of_frame_description('MTAPQQMT*', 'MTAQMT*')
        ('p.(Pro4Glnfs*4)', 3, 9, 7)
        >>> out_of_frame_description('MTAPQQT*', 'MTAQQMT*')
        ('p.(Pro4Glnfs*5)', 3, 8, 8)
        >>> out_of_frame_description('MTAPQQT*', 'MTAQQMT')
        ('p.(Pro4Glnfs*?)', 3, 8, 7)

    @arg s1: The original protein.
    @type s1: unicode
    @arg s2: The mutated protein.
    @type s2: unicode

    @return: A tuple of:
        - unicode ; Protein description of the change.
        - int     ; First position of the change.
        - int     ; Last position of the first protein.
        - int     ; Last position of the second protein.
    @rtype: tuple(unicode, int, int, int)

    @todo: More intelligently handle longest_common_prefix().
    """
    s1_seq = s1.rstrip('*')
    s2_seq = s2.rstrip('*')
    lcp = len(longest_common_prefix(s1_seq, s2_seq))

    if lcp == len(s2_seq): # NonSense mutation.
        if lcp == len(s1_seq): # Is this correct?
            return ('p.(=)', 0, 0, 0)
        return ('p.(%s%i*)' % (seq3(s1[lcp]), lcp + 1), lcp, len(s1), lcp)
    if lcp == len(s1_seq):
        # http://www.hgvs.org/mutnomen/FAQ.html#nostop
        stop = unicode(abs(len(s1_seq) - len(s2_seq))) if '*' in s2 else '?'

        return ('p.(*%i%sext*%s)' % \
                (len(s1_seq) + 1, seq3(s2[len(s1_seq)]), stop),
                len(s1_seq), len(s1), len(s2))

    # http://www.hgvs.org/mutnomen/FAQ.html#nostop
    stop = unicode(len(s2_seq) - lcp + 1) if '*' in s2 else '?'

    return ('p.(%s%i%sfs*%s)' % \
            (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp]), stop),
            lcp, len(s1), len(s2))
#out_of_frame_description


def protein_description(cds_stop, s1, s2):
    """
    Wrapper function for the in_frame_description() and
    out_of_frame_description() functions. It uses the value cds_stop to
    decide which one to call.

        >>> protein_description(34, 'MTAPQQMT*', 'MTAQQMT*')
        ('p.(Pro4Glnfs*5)', 3, 9, 8)
        >>> protein_description(34, 'MTAPQQMT*', 'MTAQQMT')
        ('p.(Pro4Glnfs*?)', 3, 9, 7)
        >>> protein_description(33, 'MTAPQQMT*', 'MTAQQMT*')
        ('p.(Pro4del)', 3, 4, 3)
        >>> protein_description(33, 'MTAPQQMT*', 'TTAQQMT*')
        ('p.?', 0, 4, 3)

    @arg cds_stop: Position of the stop codon in c. notation (CDS length).
    @type cds_stop: int
    @arg s1: The original protein.
    @type s1: unicode
    @arg s2: The mutated protein.
    @type s2: unicode

    @return: A tuple of:
        - unicode ; Protein description of the change.
        - int     ; First position of the change.
        - int     ; Last position of the change in the first protein.
        - int     ; Last position of the change in the second protein.
    @rtype: tuple(unicode, int, int, int)
    """
    if cds_stop % 3:
        description = out_of_frame_description(s1, s2)
    else:
        description = in_frame_description(s1, s2)

    return description
#protein_description


def visualise_sequence(sequence, max_length=25, flank_size=6):
    """
    If the length of a sequence is larger than a certain maxvissize, the
    string is clipped; otherwise the string is just returned.

    @arg sequence: DNA sequence.
    @type sequence: unicode
    @arg max_length: Maximum length of visualised sequence.
    @type max_length: int
    @arg flank_size: Length of the flanks in clipped visualised sequence.
    @type flank_size: int

    @return: Either the original sequence, or an abbreviation of it.
    @rtype: unicode
    """
    if len(sequence) > max_length:
        return '%s [%ibp] %s' % (sequence[:flank_size],
                                 len(sequence) - flank_size * 2,
                                 sequence[-flank_size:])
    return sequence
#visualise_sequence


# Todo: cleanup
def _insert_tag(s, pos1, pos2, tag1, tag2):
    """
    Insert two tags (tag1 and tag2) in string s at positions pos1 and pos2
    respectively if the positions are within the length of s. If not,
    either insert one tag or do nothing. If pos1 equals pos2, don't do
    anything either.

    @arg s: A sequence.
    @type s: unicode
    @arg pos1: Position of tag1.
    @type pos1: int
    @arg pos2: Position of tag2.
    @type pos2: int
    @arg tag1: Content of tag1.
    @type tag1: unicode
    @arg tag2: Content of tag2.
    @type tag2: unicode

    @return: The original sequence, or a sequence with eiter tag1, tag2 or
             both tags inserted.
    @rtype: unicode

    @todo: Cleanup (note: only used in print_protein_html).
    """
    output = s
    block = len(s)

    # Only do something if pos1 != pos2.
    if pos1 != pos2:
        if 0 <= pos1 < block:
            # Insert tag1.
            output = output[:pos1] + tag1 + output[pos1:]
        if 0 < pos2 < block:
            # Insert tag2.
            output = output[:-(block - pos2)] + tag2 \
                     + output[-(block - pos2):]
        if pos2 == block:
            # Insert tag2. Special case, since s[:-0] would yield the empty
            # string.
            output = output + tag2

    return output
#_insert_tag


# Todo: cleanup
def print_protein_html(s, first, last, O, where, text=False):
    """
    Make a fancy representation of a protein and put it in the Output
    object under the name 'where'. The representation contains HTML tags
    and is suitable for viewing in a monospaced font.

    @arg s: A protein sequence.
    @type s: unicode
    @arg first: First position to highlight.
    @type first: int
    @arg last: Last position to highlight.
    @type last: int
    @arg O: The Output object.
    @type O: Modules.Output.Output
    @arg where: Location in the {O} object to store the representation.
    @type where: unicode

    @todo: Cleanup.
    """
    if not s: return

    block = 10        # Each block consists of 10 amino acids.
    line = 6 * block  # Each line consists of 6 blocks.

    if text:
        tag1 = '\033[91m'  # Use this tag for highlighting.
        tag2 = '\033[0m'   # And this one to end highlighting.
    else:
        tag1 = '<b style="color:#FF0000">'  # Use this tag for highlighting.
        tag2 = '</b>'                       # And this one to end highlighting.
    #else

    # The maximum length for positions is the 10_log of the length of the
    # protein.
    m = int(math.floor(math.log(len(s), 10)) + 1)
    o = 1

    # Add the first position.
    output = '%s ' % unicode(o).rjust(m)

    for i in range(0, len(s), block):
        # Add the blocks.
        output += ' ' + _insert_tag(s[i:i + block], first - i, last - i,
                                    tag1, tag2)
        if not (i + block) % line and i + block < len(s):
            # One line done.
            o += line
            O.addOutput(where, output)
            # Add the position (while escaping any potential highlighting).
            if text:
                if first < o < last:
                    output = '%s%s%s ' % (tag2, unicode(o).rjust(m), tag1)
                else:
                    output = '%s ' % unicode(o).rjust(m)
            else:
                output = \
                    '<tt style="color:#333;font-weight:normal">%s</tt> ' % \
                    unicode(o).rjust(m)

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
    @type filename: unicode

    @return: The bare filename without a path and extention.
    @rtype: unicode
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


def format_usage(usage=None, keywords={}):
    """
    Format a usage string suitable for printing to the console. Some magic
    is employed so you can usually just call this function without arguments
    to have the calling module's docstring pretty-printed.

    @kwarg usage: The string to format. If omitted, the calling module's
        docstring is used.
    @type usage: unicode
    @kwarg keywords: A dictionary of (keyword, value) pairs used to format
        the usage string. If it does not contain the key 'command', it is
        added with the value of sys.argv[0].
    @type keywords: dictionary(unicode, unicode)

    @return: Formatted usage string. This is {usage} with any entries from
        {keywords} replaced and cut-off at the first occurence of two
        consecutive empty lines.
    @rtype: unicode
    """
    if not usage:
        caller = inspect.stack()[1]
        usage = inspect.getmodule(caller[0]).__doc__
    if not 'command' in keywords:
        keywords['command'] = sys.argv[0]

    return usage.split('\n\n\n')[0].strip().format(**keywords)
#format_usage


def singleton(cls):
    """
    Decorator to define a class with a singleton instance.

    Note that this decorator makes cls a function instead of a class and
    things like super() and classmethods won't work anymore. So be carefull
    with this and certainly don't use it with subclassing.

    By Shane Hathaway, taken from PEP318 [1].

    [1] http://www.python.org/dev/peps/pep-0318/#examples
    """
    instances = {}
    @wraps(cls)
    def getinstance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]
    return getinstance
#singleton


def monkey_patch_suds():
    """
    Apply our monkey-patch for the suds package.

    For some weird reason the location http://www.w3.org/2001/xml.xsd is used
    for the XML namespace, but the W3C seems to respond too slow on that url.
    We therefore use http://www.w3.org/2009/01/xml.xsd which fixes this.

    Call this function before importing anything from the suds package. For
    example, start your file with the following:

        from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()
        from suds.client import Client
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


class AttributeDictMixin(object):
    """
    Augment classes with a Mapping interface by adding attribute access.

    Taken from `Celery <http://www.celeryproject.org/>`_
    (`celery.datastructures.AttributeDictMixin`).
    """
    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError:
            raise AttributeError(
                '{0!r} object has no attribute {1!r}'.format(
                    type(self).__name__, k))

    def __setattr__(self, key, value):
        self[key] = value


# This is used in LazyObject to define the empty wrapper.
empty = object()


# Helper for LazyObject.
def _new_method_proxy(func):
    def inner(self, *args):
        if self._wrapped is empty:
            self._setup()
        return func(self._wrapped, *args)
    return inner


class LazyObject(object):
    """
    A wrapper for another class that can be used to delay instantiation of the
    wrapped class.

    Taken from `Django <https://www.djangoproject.com/>`_
    (`django.utils.functional.LazyObject`).
    """
    _wrapped = None

    def __init__(self):
        self._wrapped = empty

    __getattr__ = _new_method_proxy(getattr)

    def __setattr__(self, name, value):
        if name == '_wrapped':
            # Assign to __dict__ to avoid infinite __setattr__ loops.
            self.__dict__['_wrapped'] = value
        else:
            if self._wrapped is empty:
                self._setup()
            setattr(self._wrapped, name, value)

    def __delattr__(self, name):
        if name == '_wrapped':
            raise TypeError('can\'t delete _wrapped.')
        if self._wrapped is empty:
            self._setup()
        delattr(self._wrapped, name)

    def _setup(self):
        """
        Must be implemented by subclasses to initialize the wrapped object.
        """
        raise NotImplementedError('subclasses of LazyObject must provide a '
                                  '_setup() method')

    # Introspection support
    __dir__ = _new_method_proxy(dir)

    # Dictionary methods support
    __getitem__ = _new_method_proxy(operator.getitem)
    __setitem__ = _new_method_proxy(operator.setitem)
    __delitem__ = _new_method_proxy(operator.delitem)

    __len__ = _new_method_proxy(len)
    __contains__ = _new_method_proxy(operator.contains)


# We try to minimize non-trivial dependencies for non-critical features. The
# setproctitle package is implemented as a C extension and hence requires a C
# compiler and the Python development headers. Here we use it as an optional
# dependency.
try:
    import setproctitle as _setproctitle
    def set_process_name(name):
        _setproctitle.setproctitle(name)
except ImportError:
    def set_process_name(name):
        pass
