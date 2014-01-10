"""
Work with the mappings of transcripts to chromosomes.

Instances of the {Converter} class convert between transcript and chromosomal
locations, using the 'Mapping' table.

The {Updater} class is an abstract base class, subclassed by {NCBIUpdater}.
Instances of {NCBIUpdater} can load NCBI mapping information from a file and
update the database with this information.
"""


from collections import defaultdict
from operator import attrgetter

from Bio.Seq import reverse_complement
import MySQLdb

from mutalyzer.db.models import Chromosome, TranscriptMapping
from mutalyzer.grammar import Grammar
from mutalyzer import Crossmap
from mutalyzer import Retriever
from mutalyzer.models import SoapMessage, Mapping, Transcript


def _construct_change(var, reverse=False):
    """
    Construct mutation description.

    @arg var: RawVar object.
    @type var: pyparsing.ParseResults
    @var reverse: Variant is on the reverse strand.
    @type reverse: bool

    @return: Description of mutation (without reference and positions).
    @rtype: string
    """
    if reverse:
        # todo: if var.Arg1 is unicode, this crashes
        try:
            arg1 = str(int(var.Arg1))
        except ValueError:
            arg1 = reverse_complement(str(var.Arg1) or '')
        try:
            arg2 = str(int(var.Arg2))
        except ValueError:
            arg2 = reverse_complement(str(var.Arg2) or '')
    else:
        arg1 = var.Arg1
        arg2 = var.Arg2

    if var.MutationType == 'subst':
        change = '%s>%s' % (arg1, arg2)
    elif var.MutationType == 'delins' and arg2:
        change = '%s%s' % (var.MutationType, arg2)
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
        @type acc: string
        @arg version: version number
        @type version: integer
        @kwarg selector: Optional gene symbol selector.
        @type selector: str
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
                TranscriptMapping.accession == acc, TranscriptMapping.version == version,
                Chromosome.assembly == self.assembly)
            if selector:
                mappings = mappings.filter_by(gene=selector)
            if selector_version:
                mappings = mappings.filter_by(transcript=selector_version)

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
        @type Type: string
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
        @type accNo: string

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
        @type accNo: string

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
        @type accNo: string
        @arg mutation: the 'mutation' (e.g. c.123C>T)
        @type mutation: string

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
        @type variant: string

        @return: var_in_g ; The variant in HGVS I{g.} notation
        @rtype: string
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
            f_change = _construct_change(variant)
            r_change = _construct_change(variant, reverse=True)
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

        if self.mapping.chromosome.organelle_type == 'mitochondrion':
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
        @type reference: string
        @kwarg version: Transcript reference version. If omitted, '0' is
            assumed.
        @type version: string

        @return: Chromosome name, orientation (+ or -), and converted
            positions.
        @rtype: tuple(string, string, list)

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
        @type variant: string

        @return: variant ;
        @rtype: string
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
        @type variant: string
        @arg rt: the return type
        @type rt: string
        @kwarg gene: Optional gene name. If given, return variant descriptions
            on all transcripts for this gene.
        @type gene: string

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
                f_change = _construct_change(variant)
                r_change = _construct_change(variant, reverse=True)

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
