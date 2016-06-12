"""
Communication with the NCBI.
"""


from urllib2 import HTTPError

from Bio import Entrez

from .config import settings
from .redisclient import client as redis


class _NegativeLinkError(Exception):
    """
    Raised when no transcript-protein link exists (used for cached negative
    links).
    """
    pass


class NoLinkError(Exception):
    """
    Raised when no transcript-protein link can be found.
    """
    pass


def _get_link_from_ncbi(source_db, target_db, match_link_name,
                        source_accession, source_version=None,
                        match_version=True):
    """
    Retrieve a linked accession number from the NCBI.

    :arg str source_db: NCBI source database.
    :arg str target_db: NCBI target database.
    :arg function match_link_name: For each link found, this function is
      called with the link name (`str`) and it should return `True` iff the
      link is to be used.
    :arg str source_accession: Accession number for which we want to find a
      link (without version number).
    :arg int source_version: Optional version number for `source_accession`.
    :arg bool match_version: If `False`, the link does not have to match
      `source_version`.

    :raises NoLinkError: If no link could be retrieved from the NCBI.

    :returns: Tuple of `(target_accession, target_version)` representing the
      link target. If `source_version` is not specified or `match_version` is
      `False`, `target_version` can be `None`.
    :rtype: tuple(str, int)
    """
    # This procedure uses sequence GI numbers as returned by the Entrez
    # webservice. The NCBI has announced phasing out the use of GIs in favor
    # of `accession.version` combinations.
    #
    # https://www.ncbi.nlm.nih.gov/news/03-02-2016-phase-out-of-GI-numbers/
    #
    # From that page:
    #
    # > NCBI services that accept GI's as input will continue to be supported,
    # > and NCBI will be adding support for accession.version identifiers to
    # > all services that currently do not support them.
    #
    # At the moment (2016-06-01) only GIs are returned by the calls we use
    # below, so we cannot move to `accession.version` here. This is fine for
    # now, but should be reconsidered at some point.
    Entrez.email = settings.EMAIL

    # If we are currently strictly matching on version, we can try again if
    # no result is found. Otherwise, we just report failure.
    def fail_or_retry():
        if source_version is None or match_version:
            raise NoLinkError()
        return _get_link_from_ncbi(source_db, target_db, match_link_name,
                                   source_accession, source_version=None,
                                   match_version=False)

    if source_version is None:
        source = source_accession
    else:
        source = '%s.%d' % (source_accession, source_version)

    # Find source record.
    try:
        handle = Entrez.esearch(db=source_db, term=source)
    except HTTPError:
        return fail_or_retry()

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
    try:
        handle = Entrez.elink(dbfrom=source_db, db=target_db, id=source_gi)
    except HTTPError:
        return fail_or_retry()

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
    try:
        handle = Entrez.efetch(
            db=target_db, id=target_gi, rettype='acc', retmode='text')
    except HTTPError:
        return fail_or_retry()

    target = unicode(handle.read()).strip().split('.')
    handle.close()

    target_accession = target[0]
    target_version = int(target[1]) if source_version is not None else None
    return target_accession, target_version


def _get_link_from_cache(forward_key, reverse_key, source_accession,
                         source_version=None, match_version=True):
    """
    Retrieve a linked accession number from the cache.

    :arg str forward_key: Cache key format string for the forward direction.
      The source term will be substituted in this template.
    :arg str reverse_key: Cache key format string for the reverse direction.
      The target term will be substituted in this template.
    :arg str source_accession: Accession number for which we want to find a
      link (without version number).
    :arg int source_version: Optional version number for `source_accession`.
    :arg bool match_version: If `False`, the link does not have to match
      `source_version`.

    :raises _NegativeLinkError: If a negative link was found.
    :raises NoLinkError: If no link could be found.

    :returns: Tuple of `(target_accession, target_version)` representing the
      link target. If `source_version` is not specified or `match_version` is
      `False`, `target_version` can be `None`.
    :rtype: tuple(str, int)
    """
    if source_version is not None:
        # Query cache for link with version.
        target = redis.get(forward_key %
                           ('%s.%d' % (source_accession, source_version)))
        if target == '':
            raise _NegativeLinkError()
        if target:
            target_accession, target_version = target.split('.')
            return target_accession, int(target_version)

    if source_version is None or not match_version:
        # Query cache for link without version.
        target = redis.get(forward_key % source_accession)
        if target == '':
            raise _NegativeLinkError()
        if target is not None:
            return target, None

    raise NoLinkError()


def _cache_negative_link(forward_key, source_accession, source_version=None,
                         match_version=True):
    """
    Store a negative transcript-protein link (a "no link found" result) in the
    cache.

    The cache value for a negative link is the empty string and expires in
    `NEGATIVE_LINK_CACHE_EXPIRATION` seconds.
    """
    if source_version is not None:
        # Store a negative forward link with version.
        redis.setex(forward_key %
                    ('%s.%d' % (source_accession, source_version)),
                    settings.NEGATIVE_LINK_CACHE_EXPIRATION, '')
    if source_version is None or not match_version:
        # Store a negative forward link without version.
        redis.setex(forward_key % source_accession,
                    settings.NEGATIVE_LINK_CACHE_EXPIRATION, '')


def _cache_link(forward_key, reverse_key, source_accession, target_accession,
                source_version=None, target_version=None):
    """
    Store a transcript-protein link in the cache.
    """
    # Store the link without version in both directions.
    redis.set(forward_key % source_accession, target_accession)
    redis.set(reverse_key % target_accession, source_accession)

    if source_version is not None and target_version is not None:
        # Store the link with version in both directions.
        redis.set(forward_key % ('%s.%d' % (source_accession, source_version)),
                  '%s.%d' % (target_accession, target_version))
        redis.set(reverse_key % ('%s.%d' % (target_accession, target_version)),
                  '%s.%d' % (source_accession, source_version))


def _get_link(forward_key, reverse_key, source_db, target_db, match_link_name,
              source_accession, source_version=None, match_version=True):
    """
    Combines :func:`_get_link_from_ncbi` with :func:`_get_link_from_cache` to
    add caching to transcript-protein-link retrieval.
    """
    try:
        return _get_link_from_cache(
            forward_key, reverse_key, source_accession,
            source_version=source_version, match_version=match_version)
    except _NegativeLinkError:
        # If a negative link was in the cache, we report no link found.
        raise NoLinkError()
    except NoLinkError:
        # If no link was in the cache, we continue by querying the NCBI.
        pass

    # Query NCBI service.
    try:
        target_accession, target_version = _get_link_from_ncbi(
            source_db, target_db, match_link_name, source_accession,
            source_version=source_version, match_version=match_version)
    except NoLinkError:
        # No link found, store this negative result in the cache and re-raise
        # the exception.
        _cache_negative_link(
            forward_key, source_accession, source_version=source_version,
            match_version=match_version)
        raise

    # Store the link in the cache and return the target value.
    _cache_link(
        forward_key, reverse_key, source_accession, target_accession,
        source_version=source_version, target_version=target_version)
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

    :raises NoLinkError: If no link could be found.

    :returns: Tuple of `(protein_accession, protein_version)` representing the
      linked protein. If `transcript_version` is not specified or
      `match_version` is `False`, `protein_version` can be `None`.
    :rtype: tuple(str, int)
    """
    return _get_link(
        'ncbi:transcript-to-protein:%s', 'ncbi:protein-to-transcript:%s',
        'nucleotide', 'protein',
        lambda link: link in ('nuccore_protein', 'nuccore_protein_cds'),
        transcript_accession, source_version=transcript_version,
        match_version=match_version)


def protein_to_transcript(protein_accession, protein_version=None,
                          match_version=True):
    """
    Try to find the transcript linked to a protein.

    Links are retrieved from the NCBI using their Entrez API and cached in
    Redis. Negative results (accession or link could not be found) are also
    cached, but expire after `NEGATIVE_LINK_CACHE_EXPIRATION` seconds.

    :arg str protein_accession: Accession number of the protein for which we
      want to find the transcript (without version number).
    :arg int protein_version: Protein version number. Please provide this if
      available, also if it does not need to match. This will enrich the
      cache.
    :arg bool match_version: If `False`, the link does not have to match
      `protein_version`.

    :raises NoLinkError: If no link could be found.

    :returns: Tuple of `(transcript_accession, transcript_version)`
      representing the linked transcript. If `protein_version` is not
      specified or `match_version` is `False`, `transcript_version` can be
      `None`.
    :rtype: tuple(str, int)
    """
    return _get_link(
        'ncbi:protein-to-transcript:%s', 'ncbi:transcript-to-protein:%s',
        'protein', 'nucleotide', lambda link: link == 'protein_nuccore_mrna',
        protein_accession, source_version=protein_version,
        match_version=match_version)
