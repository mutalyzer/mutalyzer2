"""
Update the database with mapping information for the given gene or genomic
reference.
"""


import argparse
from operator import attrgetter
import sys

from sqlalchemy import or_

from ..db import session
from ..db.models import Assembly, TranscriptMapping
from .. import output
from .. import Retriever


def import_gene(assembly, gene):
    """
    Update the database with information from the UCSC.

    .. todo: This is not really tested.
    """
    o = output.Output(__file__)
    o.addMessage(__file__, -1, 'INFO',
                 'Starting UCSC mapping data update (gene: %s)' % gene)

    connection = MySQLdb.connect(user='genome',
                                 host='genome-mysql.cse.ucsc.edu',
                                 db=assembly.alias)

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
    parameters = gene

    cursor = connection.cursor()
    cursor.execute(query, parameters)
    result = cursor.fetchall()
    cursor.close()

    for (acc, version, txStart, txEnd, cdsStart, cdsEnd, exonStarts, exonEnds,
         geneName, chrom, strand, protAcc) in result:
        orientation = 'reverse' if strand == '-' else 'forward'
        exon_starts = [int(i) + 1 for i in exonStarts.split(',') if i]
        exon_stops = [int(i) for i in exonEnds.split(',') if i]
        if cdsStart and cdsEnd:
            cds = cdsStart + 1, cdsEnd
        else:
            cds = None
        mapping = TranscriptMapping.create_or_update(
            chrom, 'refseq', acc, gene.name, orientation, txStart + 1, txEnd,
            exon_starts, exon_stops, 'ucsc', cds=cds, version=int(version))
        session.add(mapping)

    session.commit()

    o.addMessage(__file__, -1, 'INFO',
                 'UCSC mapping data update end (gene: %s)' % gene)


def import_reference(assembly, reference):
    """
    Update the database with information from the given reference.

    .. todo: Also report how much was added/updated.

    .. note: Currently no exon locations are supported, this has only been
       tested on mtDNA.
    """
    o = output.Output(__file__)
    o.addMessage(__file__, -1, 'INFO',
                 'Starting reference mapping data update (reference: %s)' % reference)

    chromosome = assembly.chromosomes.filter_by(name='chrM').one()

    retriever = Retriever.GenBankRetriever(o)
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

    o.addMessage(__file__, -1, 'INFO',
                 'Reference mapping data update end (reference: %s)' % reference)


def main():
    """
    Command-line interface to the mapping importer.
    """
    database_parser = argparse.ArgumentParser(add_help=False)
    database_parser.add_argument(
        '-a', '--assembly', metavar='ASSEMBLY', dest='assembly_name_or_alias',
        default='hg19', help='assembly to import to (default: hg19)')

    parser = argparse.ArgumentParser(
        description='Mutalyzer mapping importer.',
        epilog='This program is intended to be run manually whenever '
        'transcript mappings for specific genes are required that are not '
        'yet in our database (i.e., they are not yet available from the '
        'NCBI, or they are mtDNA genes). It will not overwrite '
        'transcript/version entries that are already in our database.',
        parents=[database_parser])

    subparsers = parser.add_subparsers(
        title='subcommands', dest='subcommand', help='subcommand help')

    p = subparsers.add_parser(
        'gene', help='import gene', parents=[database_parser],
        description='Import gene mapping from the UCSC.')
    p.add_argument(
        'gene', metavar='GENE_SYMBOL',
        help='gene to import all transcript mappings for from the UCSC '
        'database (example: TTN)')

    p = subparsers.add_parser(
        'reference', help='import reference', parents=[database_parser],
        description='Import genomic reference file')
    p.add_argument(
        'reference', metavar='FILE',
        help='genomic reference to import all genes from (example: '
        'NC_012920.1)')

    args = parser.parse_args()

    assembly = Assembly.query \
        .filter(or_(Assembly.name == args.assembly_name_or_alias,
                    Assembly.alias == args.assembly_name_or_alias)).first()
    if not assembly:
        sys.stderr.write('Invalid assembly: %s\n' % (assembly_name_or_alias,))
        sys.exit(1)

    if args.subcommand == 'gene':
        import_gene(assembly, args.gene)
    if args.subcommand == 'reference':
        import_reference(assembly, args.reference)


if __name__ == '__main__':
    main()
