
from datetime import datetime

from sqlalchemy import (Column, ForeignKey, Integer, String, UniqueConstraint)
from sqlalchemy.orm import relationship

from mutalyzer import dbgb


class Reference(dbgb.Base):
    """
    Cached information about a reference sequence.
    """
    __tablename__ = 'references'
    # __table_args__ = {'mysql_engine': 'InnoDB', 'mysql_charset': 'utf8'}
    __table_args__ = (UniqueConstraint('accession', 'version',
                                       'checksum_reference',
                                       'checksum_sequence',
                                       name='_accession_version_checksum'),
                      )

    id = Column(Integer, primary_key=True)

    transcripts = relationship("Transcript")
    transcripts_fuzzy = relationship("TranscriptFuzzy")

    #: Accession number for this reference.
    # If multiple accessions in the .gb file only the first one is considered.
    accession = Column(String(100), nullable=False, index=True)

    #: Version number.
    version = Column(String(5), nullable=False, index=True)

    #: Reference file MD5 checksum.
    checksum_reference = Column(String(32), nullable=False, index=True)

    #: Stored sequence file MD5 checksum.
    checksum_sequence = Column(String(32), nullable=False, index=True)

    #: From where was the reference gathered. At the moment only 'NCBI eutils'.
    source = Column(String)

    #: Entry creation date and time.
    date_added = Column(String)

    #: Date presend in the annotations section of the genbank file.
    date_annotation = Column(String(20))

    #: Sequence length.
    length = Column(Integer, nullable=False)

    #: In vivo sequence molecule type, found within the 'source' feature
    # qualifiers. Value format: "genomic DNA", "genomic RNA", "mRNA", "tRNA",
    # "rRNA", "other RNA", "other DNA", "transcribed RNA", "viral cRNA",
    # "unassigned DNA", "unassigned RNA".
    mol_type = Column(String(20), nullable=False)

    #: For a CDS feature the 'transl_table' qualifier defines the genetic code
    # table used if other than the universal genetic code table; genetic code
    # exceptions outside the range of the specified tables is reported in
    # /transl_except qualifier;
    # It should be extracted from a CDS feature of the record and should be
    # consistent among all the CDS features.
    transl_table = Column(String(20), nullable=False)

    def __init__(self, accession, version,
                 checksum_reference, checksum_sequence,
                 source, date_annotation, length, mol_type,
                 transl_table):
        self.accession = accession
        self.version = version
        self.checksum_reference = checksum_reference
        self.checksum_sequence = checksum_sequence
        self.source = source
        self.date_added = datetime.now()
        self.date_annotation = date_annotation
        self.length = length
        self.mol_type = mol_type
        self.transl_table = transl_table

    def __repr__(self):
        return '<Reference %r>' % self.accession


class Transcript(dbgb.Base):
    """
    Mapping of a gene transcript on a chromosome.
    """
    __tablename__ = 'transcripts'
    __table_args__ = {'mysql_engine': 'InnoDB', 'mysql_charset': 'utf8'}

    id = Column(Integer, primary_key=True)

    reference_id = Column(Integer,
                          ForeignKey('references.id'),
                          nullable=False)

    exons = relationship("Exon")

    #: Accession number for the mRNA transcript
    transcript_accession = Column(String(20), nullable=False, index=True)

    #: Version number for the mRNA transcript
    transcript_version = Column(String(5), index=True)

    #: Accession number for the CDS
    protein_accession = Column(String(20), nullable=False, index=True)

    #: Version number for the CDS
    protein_version = Column(String(5), index=True)

    #: Gene symbol (e.g., ``DMD``, ``PSEN1``, ``TRNS1``).
    gene = Column(String(20), nullable=False)

    #: Gene symbol (e.g., ``DMD``, ``PSEN1``, ``TRNS1``).
    gene_synonym = Column(String)

    #: The orientation of the transcript on the chromosome.
    #: For the moment: '+' for 5' to 3' and '-' for 3' to 5'.
    strand = Column(String(2), nullable=False)

    #: The start position of the transcript on the chromosome (one-based,
    #: inclusive, in chromosomal orientation).
    transcript_start = Column(Integer, nullable=False, index=True)

    #: The stop position of the transcript on the chromosome (one-based,
    #: inclusive, in chromosomal orientation).
    transcript_stop = Column(Integer, nullable=False, index=True)

    #: The CDS start position of the transcript on the chromosome (one-based,
    #: inclusive, in chromosomal orientation).
    cds_start = Column(Integer, nullable=False, index=True)

    #: The CDS stop position of the transcript on the chromosome (one-based,
    #: inclusive).
    cds_stop = Column(Integer, nullable=False, index=True)

    #: The mRNA product.
    transcript_product = Column(String)

    #: The CDS product.
    protein_product = Column(String)

    #: The exon start positions of the transcript on the chromosome
    #: (one-based, inclusive, in chromosomal orientation).
    exons_start = Column(String)

    #: The exon start positions of the transcript on the chromosome
    #: (one-based, inclusive, in chromosomal orientation).
    exons_stop = Column(String)

    #: The locus_tag qualifier, i.e., "a submitter-supplied, systematic,
    #: stable identifier for a gene and its associated features, used for
    #: tracking purposes
    locus_tag = Column(String(20))

    #: For CDS features the 'codon_start' qualifier has valid value of 1 or 2
    #: or 3, indicating the offset at which the first complete codon of a coding
    #: feature can be found, relative to the first base of that feature.
    #: It should be extracted from a CDS during transcript to protein linking.
    codon_start = Column(String(20))

    #: Entrez Gene Database (replaces NCBI Locus Link).
    #: Example: /db_xref="GeneID:3054987".
    #: From: http://www.insdc.org/db_xref.html
    cds_db_xref_geneid = Column(String(20))

    #: A unique ID provided by the HGNC (HUGO Gene Nomenclature Committee) for
    #: each gene with an approved symbol. IDs are of the format HGNC:n, where
    #: n is a unique number.
    #: Example: /db_xref="HGNC:11818" for approved gene symbol TIMM8B with the
    #: name: 'translocase of inner mitochondrial membrane 8 homolog B'.
    cds_db_xref_hgnc = Column(String(20))

    def __init__(self, transcript_accession, transcript_version,
                 protein_accession, protein_version, gene, gene_synonym,
                 strand, transcript_start, transcript_stop,
                 cds_start, cds_stop,
                 transcript_product, protein_product,
                 exons_start, exons_stop,
                 locus_tag, codon_start,
                 cds_db_xref_geneid, cds_db_xref_hgnc):
        self.transcript_accession = transcript_accession
        self.transcript_version = transcript_version
        self.protein_accession = protein_accession
        self.protein_version = protein_version
        self.gene = gene
        self.gene_synonym = gene_synonym
        self.strand = strand
        self.transcript_start = transcript_start
        self.transcript_stop = transcript_stop
        self.cds_start = cds_start
        self.cds_stop = cds_stop
        self.transcript_product = transcript_product
        self.protein_product = protein_product
        self.exons_start = exons_start
        self.exons_stop = exons_stop
        self.locus_tag = locus_tag
        self.codon_start = codon_start
        self.cds_db_xref_geneid = cds_db_xref_geneid
        self.cds_db_xref_hgnc = cds_db_xref_hgnc

    def __repr__(self):
        return ('<Transcript %s.%s gene=%s '
                'protein %s.%s>'
                % (self.transcript_accession, self.transcript_version,
                   self.gene, self.protein_accession, self.protein_version))

    @property
    def coding(self):
        """
        Set to `True` iff the transcript is coding.
        """
        return self.cds_start is not None and self.cds_stop is not None

    @property
    def cds(self):
        """
        Tuple of CDS start and stop positions on the chromosome, or `None` if
        the transcript is non-coding.
        """
        if self.coding:
            return self.cds_start, self.cds_stop
        return None

    @cds.setter
    def cds(self, cds):
        self.cds_start, self.cds_stop = cds or (None, None)

    def get_reference(self, include_version=True):
        """
        Get fully qualified reference for this transcript.

        You would usually want to use the simpler :attr:`reference` property
        instead, except if the accession number may not include version number
        (which we consider bad practice).
        """
        if include_version and self.version:
            accession = '%s.%i' % (self.accession, self.version)
        else:
            accession = self.accession

        if self.select_transcript:
            if self.reference_type == 'lrg':
                selector = 't%d' % self.transcript_accession
            elif self.transcript_accession:
                selector = '(%s_v%.3i)' % (self.gene, self.transcript_accession)
            else:
                selector = '(%s)' % self.gene
        else:
            selector = ''

        return '%s%s' % (accession, selector)

    @property
    def reference(self):
        """
        Fully qualified reference for this transcript.
        """
        return self.get_reference()


# Index('transcript_mapping_transcript',
#       TranscriptMapping.accession, TranscriptMapping.version,
#       TranscriptMapping.gene, TranscriptMapping.transcript,
#       TranscriptMapping.chromosome_id,
#       unique=True)


class TranscriptFuzzy(dbgb.Base):
    """
    Mapping of a gene transcript on a chromosome.
    """
    __tablename__ = 'transcripts_fuzzy'
    __table_args__ = {'mysql_engine': 'InnoDB', 'mysql_charset': 'utf8'}

    id = Column(Integer, primary_key=True)

    reference_id = Column(Integer,
                          ForeignKey('references.id'),
                          nullable=False)

    #: Accession number for the mRNA transcript
    transcript_accession = Column(String(20), nullable=False)

    #: Version number for the mRNA transcript
    transcript_version = Column(String(5))

    #: Accession number for the CDS
    protein_accession = Column(String(20), nullable=False)

    #: Version number for the CDS
    protein_version = Column(String(5))

    #: Gene symbol (e.g., ``DMD``, ``PSEN1``, ``TRNS1``).
    gene = Column(String(20), nullable=False)

    #: Gene symbol (e.g., ``DMD``, ``PSEN1``, ``TRNS1``).
    gene_synonym = Column(String)

    #: The orientation of the transcript on the chromosome.
    #: For the moment: '+' for 5' to 3' and '-' for 3' to 5'.
    strand = Column(String(2), nullable=False)

    #: The start position of the transcript on the chromosome (one-based,
    #: inclusive, in chromosomal orientation).
    transcript_start = Column(String(64), nullable=False)

    #: The stop position of the transcript on the chromosome (one-based,
    #: inclusive, in chromosomal orientation).
    transcript_stop = Column(String(64), nullable=False)

    #: The CDS start position of the transcript on the chromosome (one-based,
    #: inclusive, in chromosomal orientation).
    cds_start = Column(String(64))

    #: The CDS stop position of the transcript on the chromosome (one-based,
    #: inclusive).
    cds_stop = Column(String(64))

    #: The exons start position of the transcript on the chromosome (one-based,
    #: inclusive, in chromosomal orientation).
    exons_start = Column(String)

    #: The exons stop position of the transcript on the chromosome (one-based,
    #: inclusive).
    exons_stop = Column(String)

    #: The locus_tag qualifier, i.e., "a submitter-supplied, systematic,
    #: stable identifier for a gene and its associated features, used for
    #: tracking purposes
    locus_tag = Column(String(20))

    #: For CDS features the 'codon_start' qualifier has valid value of 1 or 2
    #: or 3, indicating the offset at which the first complete codon of a coding
    #: feature can be found, relative to the first base of that feature.
    #: It should be extracted from a CDS during transcript to protein linking.
    codon_start = Column(String(20))

    #: Entrez Gene Database (replaces NCBI Locus Link).
    #: Example: /db_xref="GeneID:3054987".
    #: From: http://www.insdc.org/db_xref.html
    cds_db_xref_geneid = Column(String(20))

    #: A unique ID provided by the HGNC (HUGO Gene Nomenclature Committee) for
    #: each gene with an approved symbol. IDs are of the format HGNC:n, where
    #: n is a unique number.
    #: Example: /db_xref="HGNC:11818" for approved gene symbol TIMM8B with the
    #: name: 'translocase of inner mitochondrial membrane 8 homolog B'.
    cds_db_xref_hgnc = Column(String(20))

    def __init__(self, transcript_accession, transcript_version,
                 protein_accession, protein_version, gene, gene_synonym,
                 strand, transcript_start, transcript_stop,
                 cds_start, cds_stop, exons_start, exons_stop,
                 locus_tag, codon_start,
                 cds_db_xref_geneid, cds_db_xref_hgnc
                 ):
        self.transcript_accession = transcript_accession
        self.transcript_version = transcript_version
        self.protein_accession = protein_accession
        self.protein_version = protein_version
        self.gene = gene
        self.gene_synonym = gene_synonym
        self.strand = strand
        self.transcript_start = transcript_start
        self.transcript_stop = transcript_stop
        self.cds_start = cds_start
        self.cds_stop = cds_stop
        self.exons_start = exons_start
        self.exons_stop = exons_stop
        self.locus_tag = locus_tag
        self.codon_start = codon_start
        self.cds_db_xref_geneid = cds_db_xref_geneid
        self.cds_db_xref_hgnc = cds_db_xref_hgnc


class Exon(dbgb.Base):
    """
    Mapping of a gene transcript on a chromosome.

    .. note:: All positions are ordered according to the chromosome,
       irrespective of transcript orientation. For example, `start` is always
       smaller than `stop`.
    """
    __tablename__ = 'exons'

    id = Column(Integer, primary_key=True)

    transcript_id = Column(Integer,
                           ForeignKey('transcripts.id'),
                           nullable=False,
                           index=True)

    #: Start position.
    start = Column(Integer, nullable=False, index=True)

    #: Stop/end position.
    stop = Column(Integer, nullable=False, index=True)
