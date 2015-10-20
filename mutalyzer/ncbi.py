"""
Communication with the NCBI.
"""


import functools

from Bio import Entrez

from .config import settings
from .redisclient import client as redis


def _get_link(source_accession, source_db, target_db, match_link_name):
    """
    Retrieve a linked accession number from the NCBI.

    :arg str source_accession: Accession number for which we want to find a
      link (without version number).
    :arg str source_db: NCBI source database.
    :arg str target_db: NCBI target database.
    :arg function match_link_name: For each link found, this function is
      called with the link name (`str`) and it should return `True` iff the
      link is to be used.

    :returns: Linked accession number (without version number) or `None` if no
      link can be found.
    :rtype: str
    """
    Entrez.email = settings.EMAIL
    handle = Entrez.esearch(db=source_db, term=source_accession)
    try:
        result = Entrez.read(handle)
    except Entrez.Parser.ValidationError:
        return None
    finally:
        handle.close()

    try:
        source_gi = unicode(result['IdList'][0])
    except IndexError:
        return None

    handle = Entrez.elink(dbfrom=source_db, db=target_db, id=source_gi)
    try:
        result = Entrez.read(handle)
    except Entrez.Parser.ValidationError:
        return None
    finally:
        handle.close()

    if not result[0]['LinkSetDb']:
        return None

    for link in result[0]['LinkSetDb']:
        if match_link_name(unicode(link['LinkName'])):
            target_gi = unicode(link['Link'][0]['Id'])
            break
    else:
        return None

    handle = Entrez.efetch(
        db=target_db, id=target_gi, rettype='acc', retmode='text')
    target_accession = unicode(handle.read()).split('.')[0]
    handle.close()
    return target_accession


def cache_link(source, target):
    """
    Decorator to add caching to link retrieval.

    :arg str source: Source database (used to construct cache key).
    :arg str target: Target database (used to construct cache key).
    """
    forward_key = 'ncbi:%s-to-%s:%%s' % (source, target)
    reverse_key = 'ncbi:%s-to-%s:%%s' % (target, source)

    def cache_source_to_target(f):
        @functools.wraps(f)
        def cached_f(accession):
            result = redis.get(forward_key % accession)
            if result is not None:
                # The empty string is a cached negative result, which we return as
                # `None`.
                return result or None

            result = f(accession)

            if result is None:
                redis.setex(forward_key % accession,
                            settings.NEGATIVE_LINK_CACHE_EXPIRATION, '')
                return None

            # We store the resulting link in both directions.
            redis.set(forward_key % accession, result)
            redis.set(reverse_key % result, accession)
            return result

        return cached_f

    return cache_source_to_target


@cache_link('transcript', 'protein')
def transcript_to_protein(transcript_accession):
    """
    Try to find the protein linked to a transcript.

    Links are retrieved from the NCBI using their Entrez API and cached in
    Redis. Negative results (accession or link could not be found) are also
    cached, but expire after `NEGATIVE_LINK_CACHE_EXPIRATION` seconds.

    :arg str transcript_accession: Accession number of the transcript for
      which we want to find the protein (without version number).

    :returns: Accession number of a protein (without version number) or `None`
      if no link can be found.
    :rtype: str
    """
    return _get_link(
        transcript_accession, 'nucleotide', 'protein',
        lambda link: link in ('nuccore_protein', 'nuccore_protein_cds'))


@cache_link('protein', 'transcript')
def protein_to_transcript(protein_accession):
    """
    Try to find the transcript linked to a protein.

    Links are retrieved from the NCBI using their Entrez API and cached in
    Redis. Negative results (accession or link could not be found) are also
    cached, but expire after `NEGATIVE_LINK_CACHE_EXPIRATION` seconds.

    :arg str protein_accession: Accession number of the protein for which we
      want to find the transcript (without version number).

    :returns: Accession number of a transcript (without version number) or
      `None` if no link can be found.
    :rtype: str
    """
    return _get_link(
        protein_accession, 'protein', 'nucleotide',
        lambda link: link == 'protein_nuccore_mrna')
