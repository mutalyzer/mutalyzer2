"""
Update the database with mapping information from the NCBI.

This program is intended to be run daily from cron. Example:

    25 6 * * *  mutalyzer-mapping-update hg19 /tmp/seq_gene.md reference
"""


# Todo: Merge this script with mapping_import, the difference between the two
#   makes no sense.


import argparse
from collections import defaultdict
from itertools import groupby
from operator import itemgetter
import sys

from sqlalchemy import or_

from ..db import session
from ..db.models import Assembly, Chromosome, TranscriptMapping


COLUMNS = ['taxonomy', 'chromosome', 'start', 'stop', 'orientation', 'contig',
           'ctg_start', 'ctg_stop', 'ctg_orientation', 'feature_name',
           'feature_id', 'feature_type', 'group_label', 'transcript',
           'evidence_code']


def update_mapping(mappings_path, group_label, assembly_name_or_alias):
    """
    Update the database with information from a NCBI `seq_gene.md` file.

    We require that this file is first sorted on the `feature_id` column
    (#11), which always contains the gene identifier, and then on the
    `chromosome` column (#2).

        sort -k 11,11 -k 2,2 seq_gene.md > seq_gene.by_gene.md

    Load NCBI mapping information from {mapping_file} into the database.

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
    """
    assembly = Assembly.query \
        .filter(or_(Assembly.name == assembly_name_or_alias,
                    Assembly.alias == assembly_name_or_alias)).first()
    if not assembly:
        sys.stderr.write('Invalid assembly: %s\n' % (assembly_name_or_alias,))
        sys.exit(1)

    chromosomes = assembly.chromosomes.all()

    def read_records(mappings_file):
        for line in mappings_file:
            if line.startswith('#'):
                continue
            record = dict(zip(COLUMNS, line.rstrip().split('\t')))

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

            yield TranscriptMapping.create_or_update(
                chromosome, 'refseq', accession, gene, orientation, start,
                stop, exon_starts, exon_stops, 'ncbi', cds=cds,
                version=version)

    with open(mappings_path) as mappings_file:
        for _, records in groupby(read_records(mappings_file),
                                  itemgetter('feature_id', 'chromosome')):
            for mapping in build_mappings(records):
                session.add(mapping)

    session.commit()


def main():
    """
    Command-line interface to the mapping updater.
    """
    parser = argparse.ArgumentParser(
        description='Mutalyzer mapping updater.')

    parser.add_argument(
        'mappings_path', metavar='FILE',
        help='Path to the NCBI mapping information (example: seq_gene.md)')
    parser.add_argument(
        'group_label', metavar='GROUP_LABEL',
        help='use only entries with this group label')
    parser.add_argument(
        'assembly', metavar='ASSEMBLY',
        help='import mappings into this assembly (example: hg19)')

    args = parser.parse_args()
    update_mapping(args.mappings_path, args.group_label, args.assembly)


if __name__ == '__main__':
    main()
