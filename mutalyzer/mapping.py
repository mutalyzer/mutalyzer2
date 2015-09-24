"""
Work with the mappings of transcripts to chromosomes.

Instances of the {Converter} class convert between transcript and chromosomal
locations, using the 'Mapping' table.

The {Updater} class is an abstract base class, subclassed by {NCBIUpdater}.
Instances of {NCBIUpdater} can load NCBI mapping information from a file and
update the database with this information.
"""


from __future__ import unicode_literals

from collections import defaultdict
from itertools import groupby
from operator import attrgetter, itemgetter

import MySQLdb

from mutalyzer.db import session
from mutalyzer.db.models import Chromosome, TranscriptMapping
from mutalyzer.grammar import Grammar
from mutalyzer.models import SoapMessage, Mapping, Transcript
from mutalyzer.output import Output
from mutalyzer import Crossmap
from mutalyzer import Retriever
from mutalyzer import util


class MapviewSortError(Exception):
    pass


def _construct_change(var, reverse=False):
    """
    Construct mutation description.

    @arg var: RawVar object.
    @type var: pyparsing.ParseResults
    @var reverse: Variant is on the reverse strand.
    @type reverse: bool

    @return: Description of mutation (without reference and positions).
    @rtype: unicode
    """
    # Note that the pyparsing parse tree yields `str('')` for nonexisting
    # attributes, so we wrap the optional attributes in `unicode()`.
    if reverse:
        try:
            arg1 = unicode(int(var.Arg1))
        except ValueError:
            arg1 = util.reverse_complement(unicode(var.Arg1))
        try:
            arg2 = unicode(int(var.Arg2))
        except ValueError:
            arg2 = util.reverse_complement(unicode(var.Arg2))
    else:
        arg1 = unicode(var.Arg1)
        arg2 = unicode(var.Arg2)

    def parse_sequence(seq):
        if not seq.Sequence:
            raise NotImplementedError('Only explicit sequences are supported '
                                      'for insertions.')
        if reverse:
            return util.reverse_complement(seq.Sequence)
        return seq.Sequence

    if var.MutationType == 'subst':
        change = '%s>%s' % (arg1, arg2)
    elif var.MutationType in ('ins', 'delins'):
        if var.SeqList:
            if reverse:
                seqs = reversed(var.SeqList)
            else:
                seqs = var.SeqList
            insertion = '[' + ';'.join(parse_sequence(seq)
                                       for seq in seqs) + ']'
        else:
            insertion = parse_sequence(var.Seq)
        change = '%s%s' % (var.MutationType, insertion)
    else:
        change = '%s%s' % (var.MutationType, arg1 or arg2 or '')

    return change
#_construct_change


class Converter(object) :
    """
    Convert between transcript and chromosomal locations.

    Search for an NM number in the MySQL database, if the version number
    matches, get the start and end positions in a variant. Translate these
    positions to I{g.} notation if the variant is in I{c.} notation or vice
    versa.

    - If no end position is present, the start position is assumed to be the
      end position.
    - If the version number is not found in the database, an error message is
      generated and a suggestion for an other version is given.
    - If the reference sequence is not found at all, an error is returned.
    - If no variant is present, the transcription start and end and CDS end in
      I{c.} notation is returned.
    - If the variant is not accepted by the nomenclature parser, a parse error
      will be printed.

    @todo: Refactor anything using {mutalyzer.models} into the {webservice}
    module.
    """
    def __init__(self, assembly, O) :
        """
        Initialise the class.

        @arg assembly: the genome build version of the organism (e.g. hg19 for
        human genome build version 19)
        @type assembly: string
        @arg O: output object
        @type O: object
        """
        self.assembly = assembly
        self.__output = O

        # Populated arguments
        self.parseTree = None
        self.crossmap = None
        self.mapping = None
    #__init__

    def _reset(self) :
        self.crossmap = None
        self.mapping = None
    #_reset

    def _parseInput(self, variant) :
        """
        Parse a variant.

        @arg variant: variant description
        @type variant: string

        @return: parsetree object
        @rtype: object
        """
        grammar = Grammar(self.__output)
        parseTree = grammar.parse(variant)
        if not parseTree :
            self.__output.addMessage(__file__, 4, "EPARSE",
                    "Could not parse the given variant")
            return None
        #if
        if not parseTree.RefSeqAcc: #In case of LRG for example
            self.__output.addMessage(__file__, 4, "EONLYGB",
                "Currently we only support GenBank Records")
            return None
        #if
        self.parseTree = parseTree
        return parseTree
    #_parseInput

    def _get_mapping(self, acc, version=None, selector=None, selector_version=None) :
        """
        Get data from database.

        @arg acc: NM_ accession number (without version)
        @type acc: unicode
        @arg version: version number
        @type version: integer
        @kwarg selector: Optional gene symbol selector.
        @type selector: unicode
        @kwarg selector_version: Optional transcript version selector.
        @type selector_version: int
        """
        versions = [m.version for m in TranscriptMapping.query.filter(
                      TranscriptMapping.accession == acc,
                      TranscriptMapping.version != None,
                      TranscriptMapping.chromosome.has(assembly=self.assembly))]

        if not versions:
            self.__output.addMessage(__file__, 4, "EACCNOTINDB",
                    "The accession number %s could not be "
                    "found in our database (or is not a transcript)." % acc)
            self.__output.addOutput("LOVDERR",
                    "Reference sequence not found.")
            return

        if version in versions:
            mappings = TranscriptMapping.query.join(Chromosome).filter(
                TranscriptMapping.accession == acc,
                TranscriptMapping.version == version,
                Chromosome.assembly == self.assembly)
            if selector:
                mappings = mappings.filter(TranscriptMapping.gene == selector)
            if selector_version:
                mappings = mappings.filter(TranscriptMapping.transcript == selector_version)

            # Todo: The 'order by chrom asc' is a quick hack to make sure we
            #   first get a primary assembly mapping instead of some haplotype
            #   mapping for genes in the HLA cluster.
            #   A better fix is to return the entire list of mappings, and/or
            #   remove all secondary mappings for the HLA cluster.
            #   See also test_converter.test_hla_cluster and bug #58.
            self.mapping = mappings.order_by(TranscriptMapping.version.desc(),
                                             Chromosome.name.asc()).first()

            if not self.mapping:
                self.__output.addMessage(__file__, 4, "EACCNOTINDB",
                                         "The accession number %s version %s "
                                         "with transcript %s version %s could not be found "
                                         "in our database." %
                                         (acc, version, selector, selector_version))
            return

        if not version:
            self.__output.addMessage(__file__, 4, "ENOVERSION",
                "Version number missing for %s" % acc)
        else :
            self.__output.addMessage(__file__, 4, "EACCNOTINDB",
                "The accession number %s version %s "
                "could not be found in our database (or is not a transcript)." %
                (acc, version))
        self.__output.addMessage(__file__, 2, "WDIFFFOUND",
            "We found these versions: %s" %
            (", ".join("%s.%s" % (acc, i) for i in sorted(versions))))

        #LOVD list of versions available
        self.__output.addOutput("LOVDERR",
                "Reference sequence version not found. "
                "Available: %s" %
            (", ".join("%s.%s" % (acc, i) for i in sorted(versions))))

        #LOVD Only newest
        #self.__output.addOutput("LOVDERR",
        #        "Reference sequence version not found. "
        #        "Available: %s.%s" % (acc, sorted(versions)[-1]))
    #_get_mapping

    def makeCrossmap(self) :
        """
        Build the crossmapper.

        @todo: ADD Error Messages

        @return: Cross ; A Crossmap object
        @rtype: object
        """
        #TODO: ADD Error Messages
        if not self.mapping:
            return None

        # Create Mutalyzer compatible exon list.
        mrna = []
        for exon in zip(self.mapping.exon_starts, self.mapping.exon_stops):
            mrna.extend(exon)

        cds = self.mapping.cds or []
        orientation = 1 if self.mapping.orientation == 'forward' else -1

        self.crossmap = Crossmap.Crossmap(mrna, cds, orientation)
        return self.crossmap
    #makeCrossmap

    @staticmethod
    def _getcoords(C, Loc, Type) :
        """
        Return main, offset and g positions given either a position in
        I{c.} or in I{g.} notation.

        @arg C: A crossmapper
        @type C: object
        @arg Loc: A location in either I{g.} or I{c.} notation
        @type Loc: object
        @arg Type: The reference type
        @type Type: unicode
        @returns: triple:
            0. Main coordinate in I{c.} notation
            1. Offset coordinate in I{c.} notation
            2. Position in I{g.} notation
        @rtype: triple (integer, integer, integer)
        """
        if Type in 'cn' :
            if Loc.IVSLoc:
                ivs_number = int(Loc.IVSLoc.IVSNumber)
                if ivs_number < 1 or ivs_number > C.numberOfIntrons():
                    # Todo: Error handling in this entire module is 'suboptimal'
                    raise Exception('Invalid intron')
                if Loc.IVSLoc.OffSgn == '+':
                    g = C.getSpliceSite(ivs_number * 2 - 1) + \
                        C.orientation * int(Loc.IVSLoc.Offset)
                else:
                    g = C.getSpliceSite(ivs_number * 2) - \
                        C.orientation * int(Loc.IVSLoc.Offset)
                main, offset = C.g2x(g)
            else:
                main = C.main2int(Loc.PtLoc.MainSgn +  Loc.PtLoc.Main)
                offset = C.offset2int(Loc.PtLoc.OffSgn +  Loc.PtLoc.Offset)
                g = C.x2g(main, offset)
                main, offset = C.g2x(g)
        else:
            g = int(Loc.PtLoc.Main)
            main, offset = C.g2x(g)

        return (main, offset, g)
    #_getcoords

    def _coreMapping(self) :
        """
        Build the Mapping ClassSerializer.

        @return: Mapping ; A ClassSerializer object
        @rtype: object
        """

        Cross = self.makeCrossmap()
        if not Cross :
            return None

        if self.parseTree.SingleAlleleVarSet:
            mutations = [v.RawVar for v in self.parseTree.SingleAlleleVarSet]
        else:
            mutations = [self.parseTree.RawVar]

        mappings = []

        for mutation in mutations:

            if not mutation.StartLoc :
                return None

            # Get the coordinates of the start position
            startmain, startoffset, start_g = \
                    self._getcoords(Cross, mutation.StartLoc,
                                    self.parseTree.RefType)

            # If there is an end position, calculate the coordinates.
            if mutation.EndLoc :
                endmain, endoffset, end_g = \
                    self._getcoords(Cross, mutation.EndLoc,
                                    self.parseTree.RefType)
            else :
                end_g, endmain, endoffset = start_g, startmain, startoffset

            # Assign these values to the Mapping ClassSerializer
            V = Mapping()
            V.startmain     = startmain
            V.startoffset   = startoffset
            V.endmain       = endmain
            V.endoffset     = endoffset
            V.start_g       = start_g
            V.end_g         = end_g
            V.mutationType  = mutation.MutationType

            mappings.append(V)

        return mappings
    #_coreMapping

    def giveInfo(self, accNo) :
        """
        Returns transcription start, transcription end and CDS stop, if
        available.

        @arg accNo: transcript (NM_) accession number (with or without version)
        @type accNo: unicode

        @return: transcription start, transcription end and CDS stop
        @rtype: triple
        """

        if '.' not in accNo :
            acc, ver = accNo, None
        else :
            acc, ver = accNo.split('.')
            ver = int(ver)
        self._get_mapping(acc, ver)
        CM = self.makeCrossmap()
        if CM :
            return CM.info()
    #giveInfo

    def mainTranscript(self, accNo) :
        """
        One of the entry points (called by the HTML publisher).

        @arg accNo: The full NM accession number (including version)
        @type accNo: unicode

        @return: T ; ClassSerializer object with the types trans_start,
        trans_stop and CDS_stop
        @rtype: object

        """

        # Initiate ClassSerializer object
        info = self.giveInfo(accNo)
        T = Transcript()
        if info :
            T.trans_start = info[0]
            T.trans_stop  = info[1]
            T.CDS_stop    = info[2]
        return T
    #mainTranscript

    def mainMapping(self, accNo, mutation) :
        """
        One of the entry points (called by the HTML publisher).

        @arg accNo: transcript (NM_) accession number (with version?)
        @type accNo: unicode
        @arg mutation: the 'mutation' (e.g. c.123C>T)
        @type mutation: unicode

        @return: ClassSerializer object
        @rtype: object
        """
        variant = "%s:%s" % (accNo, mutation)
        if self._parseInput(variant) :
            acc = self.parseTree.RefSeqAcc
            try:
                version = int(self.parseTree.Version)
            except ValueError:
                version = None
            self._get_mapping(acc, version)

        mappings = self._coreMapping()

        errors = []
        for message in self.__output.getMessages():
            soap_message = SoapMessage()
            soap_message.errorcode = message.code
            soap_message.message = message.description
            errors.append(soap_message)

        mapping = Mapping()

        mapping.messages = errors

        if not mappings:         # Something went wrong
            mapping.errorcode = len(errors)
            return mapping

        mapping.errorcode = 0

        if len(mappings) == 1:
            return mappings[0]

        mapping.mutationType = 'compound'

        # Now we have to do some tricks to combine the information from
        # possibly multiple variants into one Mapping object.
        # Todo: Maybe it is better to use self.mapping.orientation here.
        min_g = 9000000000
        max_g = 0
        forward = True
        for m in mappings:
            if m.start_g > m.end_g:
                forward = False
            if m.start_g < min_g:
                min_g = m.start_g
                min_main = m.startmain
                min_offset = m.startoffset
            if m.end_g < min_g:
                min_g = m.end_g
                min_main = m.endmain
                min_offset = m.endoffset
            if m.start_g > max_g:
                max_g = m.start_g
                max_main = m.startmain
                max_offset = m.startoffset
            if m.end_g > max_g:
                max_g = m.end_g
                max_main = m.endmain
                max_offset = m.endoffset

        if forward:
            mapping.startmain = min_main
            mapping.startoffset = min_offset
            mapping.endmain = max_main
            mapping.endoffset = max_offset
            mapping.start_g = min_g
            mapping.end_g = max_g
        else:
            mapping.startmain = max_main
            mapping.startoffset = max_offset
            mapping.endmain = min_main
            mapping.endoffset = min_offset
            mapping.start_g = max_g
            mapping.end_g = min_g

        return mapping
    #main_Mapping

    def c2chrom(self, variant) :
        """
        Converts a complete HGVS I{c.} notation into a chromosomal notation.

        @arg variant: The variant in HGVS I{c.} notation
        @type variant: unicode

        @return: var_in_g ; The variant in HGVS I{g.} notation
        @rtype: unicode
        """
        if self._parseInput(variant):
            acc = self.parseTree.RefSeqAcc
            try:
                version = int(self.parseTree.Version)
            except ValueError:
                version = None
            if self.parseTree.Gene:
                selector = self.parseTree.Gene.GeneSymbol
                selector_version = int(self.parseTree.Gene.TransVar or 1)
            else:
                selector = selector_version = None
            self._get_mapping(acc, version, selector, selector_version)

        mappings = self._coreMapping()
        if not mappings:
            return None

        if self.parseTree.SingleAlleleVarSet:
            variants = [v.RawVar for v in self.parseTree.SingleAlleleVarSet]
        else:
            variants = [self.parseTree.RawVar]

        # Construct the variant descriptions
        descriptions = []
        for variant, mapping in zip(variants, mappings):
            try:
                f_change = _construct_change(variant)
                r_change = _construct_change(variant, reverse=True)
            except NotImplementedError as e:
                self.__output.addMessage(__file__, 3, 'ENOTIMPLEMENTED',
                                         unicode(e))
                return None

            if self.mapping.orientation == 'forward':
                change = f_change
            else :
                change = r_change

            if mapping.start_g != mapping.end_g:
                if self.mapping.orientation == 'reverse':
                    last_g, first_g = mapping.start_g, mapping.end_g
                else:
                    first_g, last_g = mapping.start_g, mapping.end_g
                if last_g < first_g:
                    self.__output.addMessage(__file__, 3, 'ERANGE', 'End position '
                                             'is smaller than the begin position.')
                    return None
                descriptions.append('%s_%s%s' % (first_g, last_g, change))
            #if
            else :
                descriptions.append('%s%s' % (mapping.start_g, change))

        if len(descriptions) == 1:
            description = descriptions[0]
        else:
            description = '[' + ';'.join(descriptions) + ']'

        if self.mapping.chromosome.organelle == 'mitochondrion':
            return "%s:m.%s" % (self.mapping.chromosome.accession, description)
        else:
            return "%s:g.%s" % (self.mapping.chromosome.accession, description)
    #c2chrom

    def chromosomal_positions(self, positions, reference, version=None):
        """
        Convert c. positions to chromosomal positions.

        @arg positions: Positions in c. notation to convert.
        @type positions: list
        @arg reference: Transcript reference.
        @type reference: unicode
        @kwarg version: Transcript reference version. If omitted, '0' is
            assumed.
        @type version: unicode

        @return: Chromosome name, orientation (+ or -), and converted
            positions.
        @rtype: tuple(unicode, unicode, list)

        This only works for positions on transcript references in c. notation.
        """
        versions = [m.version for m in TranscriptMapping.query.filter(
                      TranscriptMapping.accession == reference,
                      TranscriptMapping.version != None,
                      TranscriptMapping.chromosome.has(assembly=self.assembly))]

        if version not in versions:
            return None

        self.mapping = TranscriptMapping.query \
            .join(Chromosome) \
            .filter(TranscriptMapping.accession == reference,
                    TranscriptMapping.version == version,
                    Chromosome.assembly == self.assembly) \
            .order_by(TranscriptMapping.version.desc(),
                      Chromosome.name.asc()).first()

        if not self.mapping:
            return

        mapper = self.makeCrossmap()
        if not mapper:
            return None

        chromosomal_positions = []

        for position in positions:
            main = mapper.main2int(position.MainSgn +  position.Main)
            offset = mapper.offset2int(position.OffSgn +  position.Offset)
            chromosomal_positions.append(mapper.x2g(main, offset))

        orientation = '+' if self.mapping.orientation == 'forward' else '-'

        return self.mapping.chromosome.name, orientation, chromosomal_positions
    #chromosomal_positions

    def correctChrVariant(self, variant) :
        """
        @arg variant:
        @type variant: unicode

        @return: variant ;
        @rtype: unicode
        """

        #Pre split check
        if ':' not in variant :
            self.__output.addMessage(__file__, 4, "EPARSE",
                "The variant needs a colon")
            return None

        #Remove whitespace
        variant = variant.replace(" ","")

        if variant.startswith('chr') and ':' in variant:
            preco, postco = variant.split(':', 1)

            chromosome = Chromosome.query.filter_by(assembly=self.assembly,
                                                    name=preco).first()
            if not chromosome:
                self.__output.addMessage(__file__, 4, "ENOTINDB",
                    "The accession number %s could not be found in our database (or is not a chromosome)." %
                    preco)
                return None

            variant = "%s:%s" % (chromosome.accession, postco)
        #if
        return variant
    #correctChrVariant

    def chrom2c(self, variant, rt, gene=None):
        """
        @arg variant: a variant description
        @type variant: unicode
        @arg rt: the return type
        @type rt: unicode
        @kwarg gene: Optional gene name. If given, return variant descriptions
            on all transcripts for this gene.
        @type gene: unicode

        @return: HGVS_notatations ;
        @rtype: dictionary or list
        """

        if not self._parseInput(variant) :
            return None

        acc = self.parseTree.RefSeqAcc
        version = self.parseTree.Version

        chromosome = Chromosome.query \
            .filter_by(assembly=self.assembly,
                       accession='%s.%s' % (acc, version)).first()
        if not chromosome :
            self.__output.addMessage(__file__, 4, "ENOTINDB",
                "The Accession number %s could not be found in our database (or is not a chromosome)." %
                acc)
            return None
        #if

        if self.parseTree.SingleAlleleVarSet:
            variants = [v.RawVar for v in self.parseTree.SingleAlleleVarSet]
        else:
            variants = [self.parseTree.RawVar]

        min_loc = 9000000000
        max_loc = 0
        for variant in variants:
            #FIXME This should be a proper conversion.
            loc = int(variant.StartLoc.PtLoc.Main)
            if variant.EndLoc :
                loc2 = int(variant.EndLoc.PtLoc.Main)
            else :
                loc2 = loc

            if loc2 < loc:
                self.__output.addMessage(__file__, 3, 'ERANGE', 'End position is '
                                         'smaller than the begin position.')
                return None

            min_loc = min(min_loc, loc)
            max_loc = max(max_loc, loc2)

        if gene:
            mappings = chromosome.transcript_mappings.filter_by(gene=gene)
        else:
            mappings = chromosome.transcript_mappings.filter(
                TranscriptMapping.start <= max_loc + 5000,
                TranscriptMapping.stop >= min_loc - 5000)

        HGVS_notatations = defaultdict(list)
        NM_list = []
        for mapping in mappings:
            self._reset()
            self.mapping = mapping
            core_mapping = self._coreMapping()
            if not core_mapping:
                #balen
                continue
            # construct the variant description
            accNo = "%s.%s" % (self.mapping.accession, self.mapping.version)
            if self.mapping.select_transcript:
                if self.mapping.transcript:
                    selector = '(%s_v%.3i)' % (self.mapping.gene, self.mapping.transcript)
                else:
                    selector = '(%s)' % self.mapping.gene
            else:
                selector = ''
            geneName = self.mapping.gene
            strand = self.mapping.orientation == 'forward'

            # Check if n or c type
            # Note: Originally, the below check using crossmap.info() was
            #     used (commented out now), but I do not understand this
            #     logic. Also, it breaks n. notation on non-coding mtDNA
            #     transcripts, so I replaced it with a simple .CDS check.
            #info = self.crossmap.info()
            #if info[0] == 1 and info[1] == info[2] :
            #    mtype = 'n'
            #else :
            #    mtype = 'c'
            if self.crossmap.CDS:
                mtype = 'c'
            else:
                mtype = 'n'

            mutations = []
            for variant, cmap in zip(variants, core_mapping):
                try:
                    f_change = _construct_change(variant)
                    r_change = _construct_change(variant, reverse=True)
                except NotImplementedError as e:
                    self.__output.addMessage(__file__, 4,
                                             "ENOTIMPLEMENTEDERROR", unicode(e))
                    return None

                startp = self.crossmap.tuple2string((cmap.startmain, cmap.startoffset))
                endp = self.crossmap.tuple2string((cmap.endmain, cmap.endoffset))

                if strand :
                    change = f_change
                else :
                    change = r_change
                    startp, endp = endp, startp

                if cmap.start_g != cmap.end_g :
                    loca = "%s_%s" % (startp, endp)
                else :
                    loca = "%s" % startp

                mutations.append('%s%s' % (loca, change))

            if len(mutations) == 1:
                mutation = mutations[0]
            else:
                mutation = '[' + ';'.join(mutations) + ']'

            description = "%s%s:%c.%s" % (accNo, selector, mtype, mutation)
            HGVS_notatations[geneName].append(description)
            NM_list.append(description)
        #for
        if rt == "list" :
            return NM_list
        return HGVS_notatations
    #chrom2c
#Converter


def import_from_ucsc_by_gene(assembly, gene):
    """
    Import transcript mappings for a gene from the UCSC.
    """
    connection = MySQLdb.connect(user='genome',
                                 host='genome-mysql.cse.ucsc.edu',
                                 db=assembly.alias,
                                 charset='utf8',
                                 use_unicode=True)

    query = """
        SELECT DISTINCT
          acc, version, txStart, txEnd, cdsStart, cdsEnd, exonStarts,
          exonEnds, name2 AS geneName, chrom, strand, protAcc
        FROM gbStatus, refGene, refLink
        WHERE type = "mRNA"
        AND refGene.name = acc
        AND acc = mrnaAcc
        AND name2 = %s
    """
    parameters = gene,

    cursor = connection.cursor()
    cursor.execute(query, parameters)
    result = cursor.fetchall()
    cursor.close()

    # All ranges in the UCSC tables are zero-based and open-ended. We convert
    # this to one-based, inclusive for our database.

    for (acc, version, txStart, txEnd, cdsStart, cdsEnd, exonStarts, exonEnds,
         geneName, chrom, strand, protAcc) in result:
        chromosome = assembly.chromosomes.filter_by(name=chrom).one()
        orientation = 'reverse' if strand == '-' else 'forward'
        exon_starts = [int(i) + 1 for i in exonStarts.split(',') if i]
        exon_stops = [int(i) for i in exonEnds.split(',') if i]
        if cdsStart and cdsEnd:
            cds = cdsStart + 1, cdsEnd
        else:
            cds = None
        mapping = TranscriptMapping.create_or_update(
            chromosome, 'refseq', acc, geneName, orientation, txStart + 1,
            txEnd, exon_starts, exon_stops, 'ucsc', cds=cds,
            version=int(version))
        session.add(mapping)

    session.commit()


def import_from_reference(assembly, reference):
    """
    Import transcript mappings from a genomic reference.

    .. todo: Also report how much was added/updated.

    .. note: Currently no exon locations are supported, this has only been
       tested on mtDNA.
    """
    chromosome = assembly.chromosomes.filter_by(name='chrM').one()

    output = Output(__file__)
    retriever = Retriever.GenBankRetriever(output)
    record = retriever.loadrecord(reference)

    select_transcript = len(record.geneList) > 1

    for gene in record.geneList:
        # We support exactly one transcript per gene.
        try:
            transcript = sorted(gene.transcriptList, key=attrgetter('name'))[0]
        except IndexError:
            continue

        # We use gene.location for now, it is always present and the same
        # for our purposes.
        #start, stop = transcript.mRNA.location[0], transcript.mRNA.location[1]
        start, stop = gene.location

        orientation = 'reverse' if gene.orientation == -1 else 'forward'

        try:
            cds = transcript.CDS.location
        except AttributeError:
            cds = None

        mapping = TranscriptMapping.create_or_update(
            chromosome, 'refseq', record.source_accession, gene.name,
            orientation, start, stop, [start], [stop], 'reference', cds=cds,
            select_transcript=select_transcript,
            version=int(record.source_version))
        session.add(mapping)

    session.commit()


def import_from_mapview_file(assembly, mapview_file, group_label):
    """
    Import transcript mappings from an NCBI mapview file.

    We require that this file is first sorted on the `feature_id` column
    (#11), which always contains the gene identifier, and then on the
    `chromosome` column (#2).

        sort -t $'\t' -k 11,11 -k 2,2 seq_gene.md > seq_gene.by_gene.md

    Raises :exc:`ValueError` if `mapview_file` is not sorted this way.

    The NCBI mapping file consists of entries, one per line, in order of
    their location in the genome (more specifically by start location).
    Every entry has a 'group_label' column, denoting the assembly it is
    from. We only use entries where this value is `group_label`.

    There are four types of entries (for our purposes):
    - Gene: Name, identifier, and location of a gene.
    - Transcript: Name, gene id, and location of a transcript.
    - UTR: Location and transcript of a non-coding exon (or part of it).
    - CDS: Location and transcript of a coding exon (or part of it).

    A bit troublesome for us is that exons are split in UTR exons and CDS
    exons, with exons overlapping the UTR/CDS border defined as two
    separate entries (one of type UTR and one of type CDS).

    Another minor annoyance is that some transcripts (~ 15) are split over
    two contigs (NT_*). In that case, they are defined by two entries in
    the file, where we should merge them by taking the start position of
    the first and the stop position of the second.

    To complicate this annoyance, some genes (e.g. in the PAR) are mapped
    on both the X and Y chromosomes, but stored in the file just like the
    transcripts split over two contigs. However, these ones should of
    course not be merged.

    Our strategy is too sort by gene and chromosome and process the file
    grouped by these two fields.

    For transcripts without any UTR and CDS entries (seems to happen for
    predicted genes), we generate one exon spanning the entire transcript.

    All positions are one-based, inclusive, and that is what we also use in
    our database.
    """
    columns = ['taxonomy', 'chromosome', 'start', 'stop', 'orientation',
               'contig', 'ctg_start', 'ctg_stop', 'ctg_orientation',
               'feature_name', 'feature_id', 'feature_type', 'group_label',
               'transcript', 'evidence_code']

    chromosomes = assembly.chromosomes.all()

    def read_records(mapview_file):
        for line in mapview_file:
            if line.startswith('#'):
                continue
            record = dict(zip(columns, line.rstrip().split('\t')))

            # Only use records from the given assembly.
            if record['group_label'] != group_label:
                continue

            # Only use records on chromosomes we know.
            try:
                record['chromosome'] = next(c for c in chromosomes if
                                            c.name == 'chr' + record['chromosome'])
            except StopIteration:
                continue

            record['start'] = int(record['start'])
            record['stop'] = int(record['stop'])

            yield record

    def build_mappings(records):
        # We structure the records per transcript and per record type. This is
        # generalized to a list of records for each type, but we expect only
        # one GENE record (with `-` as transcript value).
        # Note that there can be more than one RNA record per transcript if it
        # is split over different reference contigs.
        by_transcript = defaultdict(lambda: defaultdict(list))
        for r in records:
            by_transcript[r['transcript']][r['feature_type']].append(r)

        gene = by_transcript['-']['GENE'][0]['feature_name']

        for transcript, by_type in by_transcript.items():
            if transcript == '-':
                continue
            accession, version = transcript.split('.')
            version = int(version)
            chromosome = by_type['RNA'][0]['chromosome']
            orientation = 'reverse' if by_type['RNA'][0]['orientation'] == '-' else 'forward'
            start = min(t['start'] for t in by_type['RNA'])
            stop = max(t['stop'] for t in by_type['RNA'])

            exon_starts = []
            exon_stops = []
            cds_positions = []
            for exon in sorted(by_type['UTR'] + by_type['CDS'],
                               key=itemgetter('start')):
                if exon_stops and exon_stops[-1] > exon['start'] - 1:
                    # This exon starts before the end of the previous exon. We
                    # have no idea what to do in this case, so we ignore it.
                    # The number of transcripts affected is very small (e.g.,
                    # NM_031860.1 and NM_001184961.1 in the GRCh37 assembly).
                    continue
                if exon['feature_type'] == 'CDS':
                    cds_positions.extend([exon['start'], exon['stop']])
                if exon_stops and exon_stops[-1] == exon['start'] - 1:
                    # This exon must be merged with the previous one because
                    # it is split over two entries (a CDS part and a UTR part
                    # or split over different reference contigs).
                    exon_stops[-1] = exon['stop']
                else:
                    exon_starts.append(exon['start'])
                    exon_stops.append(exon['stop'])

            if cds_positions:
                cds = min(cds_positions), max(cds_positions)
            else:
                cds = None

            # If no exons are annotated, we create one spanning the entire
            # transcript.
            if not exon_starts:
                exon_starts = [start]
                exon_stops = [stop]

            yield TranscriptMapping.create_or_update(
                chromosome, 'refseq', accession, gene, orientation, start,
                stop, exon_starts, exon_stops, 'ncbi', cds=cds,
                version=version)

    processed_keys = set()

    for key, records in groupby(read_records(mapview_file),
                                itemgetter('feature_id', 'chromosome')):
        if key in processed_keys:
            raise MapviewSortError('Mapview file must be sorted by feature_id '
                                   'and chromosome (try `sort -k 11,11 -k '
                                   '2,2`)')
        processed_keys.add(key)

        for mapping in build_mappings(records):
            session.add(mapping)

    session.commit()
