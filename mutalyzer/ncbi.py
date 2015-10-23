"""
Communication with the NCBI.
"""


from Bio import Entrez

from .config import settings
from .redisclient import client as redis


def _get_link(source_accession, source_db, target_db, match_link_name,
              source_version=None, match_version=True):
    """
    Retrieve a linked accession number from the NCBI.

    :arg str source_accession: Accession number for which we want to find a
      link (without version number).
    :arg str source_db: NCBI source database.
    :arg str target_db: NCBI target database.
    :arg function match_link_name: For each link found, this function is
      called with the link name (`str`) and it should return `True` iff the
      link is to be used.
    :arg int source_version: Optional version number for `source_accession`.
    :arg bool match_version: If `False`, the link does not have to match
      `source_version`.

    :returns: Tuple of `(target_accession, target_version)` representing the
      link target, or `None` if no link can be found. If `source_version` is
      not specified or `match_version` is `False`, `target_version` can be
      `None`.
    :rtype: tuple(str, int)
    """
    Entrez.email = settings.EMAIL

    # If we are currently strictly matching on version, we can try again if
    # no result is found. Otherwise, we just report failure.
    def fail_or_retry():
        if source_version is None or match_version:
            return None
        return _get_link(source_accession, source_db, target_db,
                         match_link_name, source_version=None,
                         match_version=False)

    if source_version is None:
        source = source_accession
    else:
        source = '%s.%d' % (source_accession, source_version)

    # Find source record.
    handle = Entrez.esearch(db=source_db, term=source)
    try:
        result = Entrez.read(handle)
    except Entrez.Parser.ValidationError:
        return fail_or_retry()
    finally:
        handle.close()

    try:
        source_gi = unicode(result['IdList'][0])
    except IndexError:
        return fail_or_retry()

    # Find link from source record to target record.
    handle = Entrez.elink(dbfrom=source_db, db=target_db, id=source_gi)
    try:
        result = Entrez.read(handle)
    except Entrez.Parser.ValidationError:
        return fail_or_retry()
    finally:
        handle.close()

    if not result[0]['LinkSetDb']:
        return fail_or_retry()

    for link in result[0]['LinkSetDb']:
        if match_link_name(unicode(link['LinkName'])):
            target_gi = unicode(link['Link'][0]['Id'])
            break
    else:
        return fail_or_retry()

    # Get target record.
    handle = Entrez.efetch(
        db=target_db, id=target_gi, rettype='acc', retmode='text')
    target = unicode(handle.read()).strip().split('.')
    handle.close()

    target_accession = target[0]
    target_version = int(target[1]) if source_version is not None else None
    return target_accession, target_version


def _get_link_cached(forward_key, reverse_key, source_accession, source_db,
                     target_db, match_link_name, source_version=None,
                     match_version=True):
    """
    Version of :func:`_get_link` with caching.

    :arg str forward_key: Cache key format string for the forward direction.
      The source term will be substituted in this template.
    :arg str reverse_key: Cache key format string for the reverse direction.
      The target term will be substituted in this template.

    The cache value for a negative result (no link found) is the empty string
    and expires in `NEGATIVE_LINK_CACHE_EXPIRATION` seconds.
    """
    if source_version is not None:
        # Query cache for link with version.
        target = redis.get(forward_key %
                           ('%s.%d' % (source_accession, source_version)))
        if target == '':
            return None
        if target:
            target_accession, target_version = target.split('.')
            return target_accession, int(target_version)

    if source_version is None or not match_version:
        # Query cache for link without version.
        target = redis.get(forward_key % source_accession)
        if target == '':
            return None
        if target is not None:
            return target, None

    # Query NCBI service.
    try:
        target_accession, target_version = _get_link(
            source_accession, source_db, target_db, match_link_name,
            source_version=source_version, match_version=match_version)
    except TypeError:
        # No link was found.
        if source_version is not None:
            # Store a negative forward link with version.
            redis.setex(forward_key %
                        ('%s.%d' % (source_accession, source_version)),
                        settings.NEGATIVE_LINK_CACHE_EXPIRATION, '')
        if source_version is None or not match_version:
            # Store a negative forward link without version.
            redis.setex(forward_key % source_accession,
                        settings.NEGATIVE_LINK_CACHE_EXPIRATION, '')
        return None

    # Store the link without version in both directions.
    redis.set(forward_key % source_accession, target_accession)
    redis.set(reverse_key % target_accession, source_accession)

    if source_version is not None and target_version is not None:
        # Store the link with version in both directions.
        redis.set(forward_key % ('%s.%d' % (source_accession, source_version)),
                  '%s.%d' % (target_accession, target_version))
        redis.set(reverse_key % ('%s.%d' % (target_accession, target_version)),
                  '%s.%d' % (source_accession, source_version))

    return target_accession, target_version


def transcript_to_protein(transcript_accession, transcript_version=None,
                          match_version=True):
    """
    Try to find the protein linked to a transcript.

    Links are retrieved from the NCBI using their Entrez API and cached in
    Redis. Negative results (accession or link could not be found) are also
    cached, but expire after `NEGATIVE_LINK_CACHE_EXPIRATION` seconds.

    :arg str transcript_accession: Accession number of the transcript for
      which we want to find the protein (without version number).
    :arg int transcript_version: Transcript version number. Please provide
      this if available, also if it does not need to match. This will enrich
      the cache.
    :arg bool match_version: If `False`, the link does not have to match
      `transcript_version`.

    :returns: Tuple of `(protein_accession, protein_version)` representing the
      linked protein, or `None` if no link can be found. If `match_version` is
      `False`, `protein_version` can be `None`.  TODO: can or will?
    :rtype: tuple(str, int)
    """
    return _get_link_cached(
        'ncbi:transcript-to-protein:%s', 'ncbi:protein-to-transcript:%s',
        transcript_accession, 'nucleotide', 'protein',
        lambda link: link in ('nuccore_protein', 'nuccore_protein_cds'),
        source_version=transcript_version, match_version=match_version)


def protein_to_transcript(protein_accession, protein_version=None,
                          match_version=True):
    """
    Try to find the transcript linked to a protein.

    Links are retrieved from the NCBI using their Entrez API and cached in
    Redis. Negative results (accession or link could not be found) are also
    cached, but expire after `NEGATIVE_LINK_CACHE_EXPIRATION` seconds.

    :arg str protein_accession: Accession number of the protein for which we
      want to find the transcript (without version number).
    TODO

    :returns: Accession number of a transcript (without version number) or
      `None` if no link can be found.
    :rtype: str
    """
    return _get_link_cached(
        'ncbi:protein-to-transcript:%s', 'ncbi:transcript-to-protein:%s',
        protein_accession, 'protein', 'nucleotide',
        lambda link: link == 'protein_nuccore_mrna',
        source_version=protein_version, match_version=match_version)
