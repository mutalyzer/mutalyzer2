from Bio import Entrez

from .config import settings
from .db import queries


def transcript_to_protein(transcript_accession):
    """
    Try to find the protein linked to a transcript.

    First look in our database. If a link cannot be found, try to retrieve it
    from the NCBI. Add the result to our database.

    :arg str transcript_accession: Accession number of the transcript for
      which we want to find the protein (without version number).

    :returns: Accession number of a protein (without version number) or `None`
      if no link can be found.
    :rtype: str
    """
    Entrez.email = settings.EMAIL

    link = queries.get_transcript_protein_link(transcript_accession)
    if link is not None:
        return link.protein_accession

    handle = Entrez.esearch(db='nucleotide', term=transcript_accession)
    try:
        result = Entrez.read(handle)
    except Entrez.Parser.ValidationError:
        # Todo: Log this error.
        return None
    finally:
        handle.close()

    transcript_gi = unicode(result['IdList'][0])

    handle = Entrez.elink(dbfrom='nucleotide', db='protein', id=transcript_gi)
    try:
        result = Entrez.read(handle)
    except Entrez.Parser.ValidationError:
        # Todo: Log this error.
        return None
    finally:
        handle.close()

    if not result[0]['LinkSetDb']:
        # We also cache negative results.
        queries.update_transcript_protein_link(
            transcript_accession=transcript_accession)
        return None

    protein_gi = unicode(result[0]['LinkSetDb'][0]['Link'][0]['Id'])

    handle = Entrez.efetch(
        db='protein', id=protein_gi, rettype='acc', retmode='text')
    protein_accession = unicode(handle.read()).split('.')[0]
    handle.close()

    queries.update_transcript_protein_link(
        transcript_accession=transcript_accession,
        protein_accession=protein_accession)
    return protein_accession


def protein_to_transcript(protein_accession):
    """
    Try to find the transcript linked to a protein.

    First look in our database. If a link cannot be found, try to retrieve it
    from the NCBI. Add the result to our database.

    :arg str protein_accession: Accession number of the protein for which we
      want to find the transcript (without version number).

    :returns: Accession number of a transcript (without version number) or
      `None` if no link can be found.
    :rtype: str
    """
    Entrez.email = settings.EMAIL

    link = queries.get_transcript_protein_link(protein_accession, reverse=True)
    if link is not None:
        return link.transcript_accession

    handle = Entrez.esearch(db='protein', term=protein_accession)
    try:
        result = Entrez.read(handle)
    except Entrez.Parser.ValidationError:
        # Todo: Log this error.
        return None
    finally:
        handle.close()

    if not result['IdList']:
        return None
    protein_gi = unicode(result['IdList'][0])

    handle = Entrez.elink(dbfrom='protein', db='nucleotide', id=protein_gi)
    try:
        result = Entrez.read(handle)
    except Entrez.Parser.ValidationError:
        # Todo: Log this error.
        return None
    finally:
        handle.close()

    if not result[0]['LinkSetDb']:
        # We also cache negative results.
        queries.update_transcript_protein_link(
            protein_accession=protein_accession)
        return None

    transcript_gi = ''
    for link in result[0]['LinkSetDb']:
        if unicode(link['LinkName']) == 'protein_nuccore_mrna':
            transcript_gi = unicode(link['Link'][0]['Id'])
            break

    handle = Entrez.efetch(
        db='nucleotide', id=transcript_gi, rettype='acc', retmode='text')
    transcript_accession = unicode(handle.read()).split('.')[0]
    handle.close()

    queries.update_transcript_protein_link(
        transcript_accession=transcript_accession,
        protein_accession=protein_accession)
    return transcript_accession
