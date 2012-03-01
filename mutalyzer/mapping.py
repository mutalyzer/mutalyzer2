"""
Work with the mappings of transcripts to chromosomes.

Instances of the {Converter} class convert between transcript and chromosomal
locations, using the 'Mapping' table.

The {Updater} class is an abstract base class, subclassed by {NCBIUpdater}.
Instances of {NCBIUpdater} can load NCBI mapping information from a file and
update the database with this information.
"""


from Bio.Seq import reverse_complement
from collections import defaultdict

from mutalyzer.grammar import Grammar
from mutalyzer import Db
from mutalyzer import Crossmap
from mutalyzer.models import SoapMessage, Mapping, Transcript


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
    def __init__(self, build, O) :
        """
        Initialise the class.

        @arg build: the genome build version of the organism (e.g. hg19 for
        human genome build version 19)
        @type build: string
        @arg O: output object
        @type O: object
        """
        self.build = None
        self.__output = O
        self.__database = Db.Mapping(build)

        # Populated arguments
        self.parseTree = None
        self.crossmap = None
        self.dbFields = {}
    #__init__

    def _reset(self) :
        self.crossmap = None
        self.dbFields = {}
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
        if parseTree.SingleAlleleVarSet :
            #Only simple mutations
            self.__output.addMessage(__file__, 4, "EPARSE",
                    "Can not process multiple mutation variant")
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

    def _populateFields(self, Fields) :
        """
        Create a Mutalyzer compatible exon list.

        @todo: ADD Error Messages, unlikely that CDS info is missing.

        @arg Fields: dictionary with exon start and end positions taken from the
        MySQL database
        @type Fields: dictionary

        @return: Exon list
        @rtype: list
        """
        Fields["exon_starts"] = map(int, Fields["exon_starts"].split(','))
        Fields["exon_stops"] = map(int, Fields["exon_stops"].split(','))
        assert(len(Fields["exon_starts"]) == len(Fields["exon_stops"]))

        if Fields['cds_start'] and Fields['cds_stop']:
            Fields["cds_start"] = int(Fields["cds_start"])
            Fields["cds_stop"]   = int(Fields["cds_stop"])

        # Create Mutalyzer compatible exon list
        Fields["exons"] = []
        for exon in zip(Fields["exon_starts"], Fields["exon_stops"]) :
            Fields["exons"].extend(exon)

        self.dbFields = Fields
        return Fields
    #_populateFields

    def _FieldsFromValues(self, values) :
        """
        Combines labels with the given values to a dictionary.
        (zip returns a list of tuples, where the i-th tuple contains the i-th
        element from each of the argument sequences or iterables.
        dict(arg) creates a new data dictionary, with items taken from arg.)

        @arg values: list of values take from the MySQL database
        @type values: list

        @return: dictionary with values taken from the MySQL database
        @rtype: dictionary
        """

        Fields = dict(zip(
            ("transcript", "start", "stop", "cds_start", "cds_stop",
             "exon_starts", "exon_stops", "gene",
             "chromosome", "orientation", "protein", "version"),
            values))
        return self._populateFields(Fields)
    #_FieldsFromValues

    def _FieldsFromDb(self, acc, version) :
        """
        Get data from database and populate dbFields dict.

        @arg acc: NM_ accession number (without version)
        @type acc: string
        @arg version: version number
        @type version: integer
        """

        if not version :
            version = 0
        version = int(version)
        versions = self.__database.get_NM_version(acc)
        if not versions :
            self.__output.addMessage(__file__, 4, "EACCNOTINDB",
                    "The accession number %s could not be "
                    "found in our database (or is not a transcript)." % acc)
            self.__output.addOutput("LOVDERR",
                    "Reference sequence not found.")
            return None     #Explicit return of None in case of error
        #if
        else :
            if version in versions :
                Values = self.__database.getAllFields(acc, version)
                return self._FieldsFromValues(Values)
            #if
            if not version :
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
            return None
        #else
    #_FieldsFromDb

    def makeCrossmap(self) :
        """
        Build the crossmapper.

        @todo: ADD Error Messages

        @return: Cross ; A Crossmap object
        @rtype: object
        """

        #TODO: ADD Error Messages
        if not self.dbFields: return None

        CDS = []
        if self.dbFields["cds_start"] and self.dbFields["cds_stop"]:
            CDS = [self.dbFields["cds_start"], self.dbFields["cds_stop"]]

        mRNA = self.dbFields["exons"]

        # Convert the strand information to orientation.
        orientation = 1
        if self.dbFields["orientation"] == '-' :
            orientation = -1

        # Build the crossmapper.
        self.crossmap = Crossmap.Crossmap(mRNA, CDS, orientation)
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
        if Type == 'c' :
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

        mutation = self.parseTree.RawVar

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

        return V
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
        self._FieldsFromDb(acc, ver)
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
            version = self.parseTree.Version
            self._FieldsFromDb(acc, version)

        mapping = self._coreMapping()

        errors = []
        for message in self.__output.getMessages():
            soap_message = SoapMessage()
            soap_message.errorcode = message.code
            soap_message.message = message.description
            errors.append(soap_message)

        if mapping is None :         # Something went wrong
            mapping = Mapping()
            mapping.errorcode = len(errors)
        else :
            mapping.errorcode = 0

        mapping.messages = errors

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

        if self._parseInput(variant) :
            acc = self.parseTree.RefSeqAcc
            version = self.parseTree.Version
            self._FieldsFromDb(acc, version)
        #if
        M = self._coreMapping()
        if M is None :
            return None

        # construct the variant description
        chromAcc = self.__database.chromAcc(self.dbFields["chromosome"])
        f_change = self._constructChange(False)
        r_change = self._constructChange(True)
        if self.dbFields["orientation"] == "+" :
            change = f_change
        else :
            change = r_change

        if M.start_g != M.end_g :
            if self.dbFields["orientation"] == '+' :
                var_in_g = "g.%s_%s%s" % (M.start_g, M.end_g, change)
            else :
                var_in_g = "g.%s_%s%s" % (M.end_g, M.start_g, change)
        #if
        else :
            var_in_g = "g.%s%s" % (M.start_g, change)

        return "%s:%s" % (chromAcc, var_in_g)
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
        if not version:
            version = 0
        version = int(version)

        versions = self.__database.get_NM_version(reference)

        if version not in versions:
            return None

        values = self.__database.getAllFields(reference, version)
        self._FieldsFromValues(values)

        mapper = self.makeCrossmap()
        if not mapper:
            return None

        chromosomal_positions = []

        for position in positions:
            main = mapper.main2int(position.MainSgn +  position.Main)
            offset = mapper.offset2int(position.OffSgn +  position.Offset)
            chromosomal_positions.append(mapper.x2g(main, offset))

        return self.dbFields['chromosome'], self.dbFields['orientation'], chromosomal_positions
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
            chrom = self.__database.chromAcc(preco)
            if chrom is None :
                self.__output.addMessage(__file__, 4, "ENOTINDB",
                    "The accession number %s could not be found in our database (or is not a chromosome)." %
                    preco)
                return None
            #if
            else :
                variant = "%s:%s" % (chrom, postco)
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
        chrom = self.__database.chromName("%s.%s" % (acc, version))
        if not chrom :
            self.__output.addMessage(__file__, 4, "ENOTINDB",
                "The Accession number %s could not be found in our database (or is not a chromosome)." %
                acc)
            return None
        #if
        f_change = self._constructChange(False)
        r_change = self._constructChange(True)

        #FIXME This should be a proper conversion.
        loc = int(self.parseTree.RawVar.StartLoc.PtLoc.Main)
        if self.parseTree.RawVar.EndLoc :
            loc2 = int(self.parseTree.RawVar.EndLoc.PtLoc.Main)
        else :
            loc2 = loc

        if gene:
            transcripts = self.__database.get_TranscriptsByGeneName(gene)
        else:
            transcripts = self.__database.get_Transcripts(chrom, loc-5000, loc2+5000, 1)

        HGVS_notatations = defaultdict(list)
        NM_list = []
        for transcript in transcripts :
            self._reset()
            self._FieldsFromValues(transcript)
            if self.dbFields['chromosome'] != chrom:
                # Could be the case if we got transcripts by gene name
                continue
            M = self._coreMapping()
            if M is None :
                #balen
                continue
            # construct the variant description
            accNo = "%s.%s" % (self.dbFields["transcript"],self.dbFields["version"])
            geneName = self.dbFields["gene"] or ""
            strand = self.dbFields["orientation"] == '+'
            startp = self.crossmap.tuple2string((M.startmain, M.startoffset))
            endp = self.crossmap.tuple2string((M.endmain, M.endoffset))

            if strand :
                change = f_change
            else :
                change = r_change
                startp, endp = endp, startp

            #Check if n or c type
            info = self.crossmap.info()
            if info[0] == '1' and info[1] == info[2] :
                mtype = 'n'
            else :
                mtype = 'c'

            if M.start_g != M.end_g :
                loca = "%s_%s" % (startp, endp)
            else :
                loca = "%s" % startp

            variant = "%s:%c.%s%s" % (accNo, mtype, loca, change)
            HGVS_notatations[geneName].append(variant)
            NM_list.append(variant)
        #for
        if rt == "list" :
            return NM_list
        return HGVS_notatations
    #chrom2c

    def _constructChange(self, revc = False) :
        """
        @todo document me

        @arg revc:
        @type revc:

        @return:
        @rtype: string
        """

        p = self.parseTree
        if not p or p.SingleAlleleVarSet :
            return None
        var = p.RawVar

        if revc :
            # todo: if var.Arg1 is unicode, this crashes
            try:
                arg1 = str(int(var.Arg1))
            except ValueError:
                arg1 = reverse_complement(var.Arg1 or "")
            try:
                arg2 = str(int(var.Arg2))
            except ValueError:
                arg2 = reverse_complement(var.Arg2 or "")
        #if
        else :
            arg1 = var.Arg1
            arg2 = var.Arg2
        #else

        if var.MutationType == "subst" :
            change = "%s>%s" % (arg1, arg2)
        else :
            change = "%s%s" % (var.MutationType, arg1 or arg2 or "")
        return change
    #_constructChange
#Converter


class Updater(object):
    """
    Abstract base class for updating the mapping information in the database.

    Subclasses should implement the {load} method, loading new mapping
    information into the 'MappingTemp' table. The {merge} method merges this
    table into the real 'Mapping' table.
    """
    def __init__(self, build):
        """
        @arg build: Human genome build (or database name), i.e. 'hg18' or
            'hg19'.
        @type build: string
        """
        self.build = build
        self.db = Db.Mapping(build)
    #__init__

    def load(self, *args, **kwargs):
        """
        The implementation of this method in subclasses should load mapping
        information in the 'MappingTemp' table.
        """
        raise NotImplementedError('Implement this method in subclasses')
    #load

    def merge(self):
        """
        Merge the 'Mapping' and 'MappingTemp' tables. The result is stored in
        the 'Mapping' table, of which a backup is created as 'MappingBackup'.

        @todo: Report how much was updated/added.
        """
        self.db.merge_update()
    #merge
#Updater


class NCBIUpdater(Updater):
    """
    Update the mapping information in the database with mapping information
    from the NCBI.

    Example usage:

        >>> updater = NCBIUpdater('hg19')
        >>> updater.load('/tmp/seq_gene.md', 'GRCh37.p2-Primary Assembly')
        >>> updater.merge()

    """
    COLUMNS = ['taxonomy', 'chromosome', 'start', 'stop', 'orientation',
               'contig', 'ctg_start', 'ctg_stop', 'ctg_orientation',
               'feature_name', 'feature_id', 'feature_type', 'group_label',
               'transcript', 'evidence_code']

    def __init__(self, build):
        """
        @arg build: Human genome build (or database name), i.e. 'hg18' or
            'hg19'.
        @type build: string
        """
        self.exon_backlog = {}
        super(NCBIUpdater, self).__init__(build)
    #__init__

    def load(self, mapping_file, assembly):
        """
        Load NCBI mapping information from {mapping_file} into the database.

        The NCBI mapping file consists of entries, one per line, in order of
        their location in the genome (more specifically by start location).
        Every entry has a 'group_name' column, denoting the assembly it is
        from. We only use entries where this value is {assembly}.

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

        Our strategy is to loop over all entries and store them in three
        temporary tables (for genes, transcripts, exons). The entries of type
        UTR and CDS are merged to correct exon entries by keeping a backlog
        of these entries that can still be modified before storing them in the
        database.

        The values from the three temporary tables are aggregated into the
        'MappingTemp' table.

        @arg mapping_file: Path to NCBI mapping information.
        @type mapping_file: string
        @arg assembly: Use only entries from this assembly (this is the
            'group_name' column in the NCBI mapping file).
        @type assembly: string
        """
        self._create_temporary_tables()
        self._import_mapping(mapping_file, assembly)
        self._aggregate_mapping()
        self._drop_temporary_tables()
    #load

    def _import_mapping(self, mapping_file, assembly):
        """
        Import mapping information from {mapping_file} into three temporary
        tables.

        @note: We issue a separate INSERT statement to the database for every
            entry. An alternative is to write everything to tab-separated
            files and load those into the database with LOAD DATA LOCAL INFILE
            statements. This alternative seems to be about twice as fast, but
            for now we stick with the simpler solution.
        """
        self.exon_backlog = {}

        with open(mapping_file, 'r') as mapping:
            for line in mapping:
                if line.startswith('#'):
                    continue
                entry = dict(zip(self.COLUMNS, line.rstrip().split('\t')))

                # Only use entries from the given assembly.
                if entry['group_label'] != assembly:
                    continue

                # Only use entries on the normal chromosomes.
                try:
                    int(entry['chromosome'])
                except ValueError:
                    if entry['chromosome'] not in 'XY':
                        continue

                if entry['feature_type'] == 'GENE':
                    self._import_gene(entry)
                elif entry['feature_type'] == 'RNA':
                    self._import_transcript(entry)
                elif entry['feature_type'] in ('UTR', 'CDS'):
                    self._import_exon(entry)

        self._import_exon_backlog()
    #_import_mapping

    def _import_gene(self, entry):
        """
        Insert a gene in the database.
        """
        self.db.ncbi_import_gene(entry['feature_id'], entry['feature_name'])
    #_import_gene

    def _import_transcript(self, entry):
        """
        Insert a transcript in the database.
        """
        self.db.ncbi_import_transcript(
            entry['feature_name'], entry['feature_id'], entry['chromosome'],
            int(entry['start']), int(entry['stop']), entry['orientation'])
    #_import_transcript

    def _import_exon(self, entry):
        """
        Instead of directly inserting each exon in the database, we keep them
        in a backlog of at most one exon per transcript. Exons are taken from
        the backlog and inserted in the database only when we passed their
        genomic stop location by more than one position.

        This way, exons in the backlog can be merged when they are on a
        UTR/CDS boundary.
        """
        cds = entry['feature_type'] == 'CDS'
        entry['start'] = int(entry['start'])
        entry['stop'] = int(entry['stop'])
        entry['protein'] = entry['feature_name'] if cds else None
        entry['cds_start'] = entry['start'] if cds else None
        entry['cds_stop'] = entry['stop'] if cds else None
        entry['cds'] = cds

        self._import_exon_backlog(entry['start'] - 1)

        try:
            previous = self.exon_backlog[entry['transcript']]
            if previous['cds'] != entry['cds'] \
                   and previous['stop'] == entry['start'] - 1:
                if previous['cds']:
                    entry['cds_start'] = previous['cds_start']
                    entry['cds_stop'] = previous['cds_stop']
                    entry['protein'] = previous['protein']
                entry['start'] = previous['start']
        except KeyError:
            pass

        self.exon_backlog[entry['transcript']] = entry
    #_import_exon

    def _import_exon_backlog(self, up_to_position=None):
        """
        Import exons from the backlog in the database. If the optional
        argument {up_to_position} is set, only import exons with a stop
        position before this value.

        We explicitely remove imported exons from the backlog, because it
        might be suboptimal to keep more than 30,000 exons in there.
        """
        for transcript, exon in self.exon_backlog.items():
            if not up_to_position or exon['stop'] < up_to_position:
                del self.exon_backlog[transcript]
                del exon['cds']
                self.db.ncbi_import_exon(
                    exon['transcript'], exon['chromosome'], exon['start'], exon['stop'],
                    exon['cds_start'], exon['cds_stop'], exon['protein'] or None)
    #_import_exon_backlog

    def _aggregate_mapping(self):
        """
        Aggregate the genes, transcripts and exons from their temporary
        tables into the 'MappingTemp' table.
        """
        self.db.ncbi_aggregate_mapping()
    #_aggregate_mapping

    def _create_temporary_tables(self):
        """
        Create temporary tables needed for loading the NCBI mapping data.
        """
        self.db.ncbi_create_temporary_tables()
    #_create_temporary_tables

    def _drop_temporary_tables(self):
        """
        Drop temporary tables needed for loading the NCBI mapping data.
        """
        self.db.ncbi_drop_temporary_tables()
    #_drop_temporary_tables
#NCBIUpdater
