"""
Mutalyzer interface for variant back translation (p. to c.).
"""


from __future__ import unicode_literals

from backtranslate.backtranslate import BackTranslate
from backtranslate.util import protein_letters_3to1
from Bio.Data import CodonTable
from extractor.variant import Allele, DNAVar, ISeq, ISeqList

from . import ncbi, Retriever
from .grammar import Grammar
from .output import Output


def backtranslate(output, description, table_id=1):
    """
    Back translation of an amino acid substitution to all possible
    causal one-nucleotide substitutions.

    :arg object output: Output object for logging.
    :arg unicode description: Amino acid substitution in HGVS format.
    :arg int table_id: Translation table id.

    :returns: List of DNA substitutions in HGVS format.
    :rtype: list(unicode)
    """
    # TODO: We currently only support NM and NP references, but in principle
    #   we could also support other references.
    # TODO: For NP references where we don't find a link to the corresponding
    #   NM, we don't check if the specified reference amino acid corresponds
    #   to the NP reference sequence.
    parse_tree = Grammar(output).parse(description)

    if not parse_tree:
        return []
    if parse_tree.RefType != 'p':
        output.addMessage(
            __file__, 3, 'ENOPROT', 'Reference type is not p. (protein).')
        return []
    if parse_tree.RawVar.MutationType != 'subst':
        output.addMessage(
            __file__, 3, 'ENOSUBST', 'Variant is not a substitution.')
        return []

    if parse_tree.Version:
        accession_number = '{}.{}'.format(parse_tree.RefSeqAcc, parse_tree.Version)
    else:
        accession_number = parse_tree.RefSeqAcc
    position = int(parse_tree.RawVar.Main)
    reference_amino_acid, amino_acid = parse_tree.RawVar.Args

    if len(reference_amino_acid) == 3:
        reference_amino_acid = protein_letters_3to1[reference_amino_acid]
    if len(amino_acid) == 3:
        amino_acid = protein_letters_3to1[amino_acid]

    bt = BackTranslate(table_id)

    # FIXME: Rancid workaround to silence fatal error raised by `loadrecord`.
    output_ = Output('')
    retriever = Retriever.GenBankRetriever(output_)

    # The genbank retriever does not (yet) support protein references, but we
    # cannot reliably distinguish between different reference types from the
    # variant description before downloading the reference.
    # Therefore, we just try to download the reference. This will succeed for
    # transcript references, but fail for protein references.
    # As a quick and dirty optimization, we shortcut this for accessions
    # starting with 'NP_', of which we know that they are protein references.
    # In the future we hope to support protein references directly.
    if accession_number.startswith('NP_'):
        genbank_record = None
    else:
        genbank_record = retriever.loadrecord(accession_number)

    if not genbank_record:
        # Assuming RefSeqAcc is an NP, try to get the corresponding NM.
        version = int(parse_tree.Version) if parse_tree.Version else None

        try:
            transcript = ncbi.protein_to_transcript(parse_tree.RefSeqAcc,
                                                    version)
        except ncbi.NoLinkError:
            pass
        else:
            if transcript[1] is not None:
                accession_number = '{}.{}'.format(*transcript)
            else:
                output.addMessage(
                    __file__, 2, 'WNOVERSION',
                    'Found corresponding nucleotide sequence, but note that '
                    'the version number is missing.')
                accession_number = transcript[0]
            genbank_record = retriever.loadrecord(accession_number)

    offset = (position - 1) * 3
    if genbank_record and genbank_record.molType == 'n':
        # Only NM for now.
        cds_loc = genbank_record.geneList[0].transcriptList[0].CDS.location
        codon = genbank_record \
            .seq[cds_loc[0] - 1:cds_loc[1]][offset:offset + 3]

        forward_table = CodonTable.unambiguous_dna_by_id[table_id]
        found_ref = forward_table.forward_table[unicode(codon)]
        if reference_amino_acid != found_ref:
            output.addMessage(
                __file__, 3, 'EREF',
                '{} not found at position {}, found {} instead.'.format(
                    reference_amino_acid, position, found_ref))

        substitutions = bt.with_dna(unicode(codon), amino_acid)
    else:
        # Assume NP.
        output.addMessage(
            __file__, 2, 'WNODNA',
            'Nucleotide reference sequence could not be found, using '
            'protein fallback method.')
        accession_number = 'UNKNOWN'

        substitutions = bt.without_dna(reference_amino_acid, amino_acid)
        if (reference_amino_acid, amino_acid) in bt.improvable():
            output.addMessage(
                __file__, 2, 'WIMPROVE',
                'The back translation for this variant can be improved by '
                'supplying a nucleotide reference sequence.')

    return ['{}:c.{}'.format(accession_number, v)
            for v in subst_to_hgvs(substitutions, offset)]


def subst_to_hgvs(substitutions, offset=0):
    """
    Translate a set of substitutions to HGVS.

    :arg dict substitutions: Set of single nucleotide substitutions indexed by
      position.
    :arg int offset: Codon position in the CDS.

    :returns: Substitutions in HGVS format.
    :rtype: set(Allele)
    """
    variants = set()

    for position in substitutions:
        for substitution in substitutions[position]:
            variants.add(Allele([DNAVar(
                start=position + offset + 1,
                end=position + offset + 1,
                sample_start=position + offset + 1,
                sample_end=position + offset + 1,
                type='subst',
                deleted=ISeqList([ISeq(sequence=substitution[0])]),
                inserted=ISeqList([ISeq(sequence=substitution[1])])
            )]))

    return variants
