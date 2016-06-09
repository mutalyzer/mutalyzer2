"""
Module for retrieving files from either the cache or the NCBI.

A hash of every retrieved file is stored in the internal database. If a
requested file is not found, but its hash is, we use additional information
to re-download the file.
"""


from __future__ import unicode_literals

import bz2
import chardet
import hashlib
import io
import os
import time
import urllib2

from Bio import Entrez
from Bio import SeqIO
from Bio.Alphabet import ProteinAlphabet
from Bio.Seq import UnknownSeq
from httplib import HTTPException, IncompleteRead
from sqlalchemy.orm.exc import NoResultFound
from xml.dom import DOMException, minidom
from xml.parsers import expat

from mutalyzer import util
from mutalyzer.config import settings
from mutalyzer.db import session
from mutalyzer.db.models import Reference
from mutalyzer.parsers import genbank
from mutalyzer.parsers import lrg


ENTREZ_MAX_TRIES = 4
ENTREZ_SLEEP = 1  # In seconds.


class Retriever(object):
    """
    Retrieve a record from either the cache or the NCBI.
    """
    def __init__(self, output):
        """
        Use variables from the configuration file for some simple
        settings. Make the cache directory if it does not exist yet.

        :arg object output: The output object.
        """
        self._output = output
        if not os.path.isdir(settings.CACHE_DIR):
            os.mkdir(settings.CACHE_DIR)
        Entrez.email = settings.EMAIL
        self.file_type = None

    def _name_to_file(self, name):
        """
        Convert an accession number to a filename.

        :arg unicode name: The accession number.

        :returns: A filename.
        :rtype: unicode
        """
        return os.path.join(
            settings.CACHE_DIR, '{}.{}.bz2'.format(name, self.file_type))

    def _write(self, raw_data, filename):
        """
        Write raw data to a compressed file.

        :arg str raw_data: The raw_data to be compressed and written.
        :arg unicode filename: The intended name of the output filename.

        :returns: The full path and name of the file written.
        :rtype: unicode
        """
        result = chardet.detect(raw_data)
        if result['confidence'] > 0.5:
            encoding = unicode(result['encoding'])
        else:
            encoding = 'utf-8'

        if not util.is_utf8_alias(encoding):
            try:
                raw_data = raw_data.decode(encoding).encode('utf-8')
            except UnicodeDecodeError:
                self._output.addMessage(
                    __file__, 4, 'ENOPARSE',
                    'Could not decode file (using {} encoding).'.format(
                        encoding))
                return None

        # Compress the data to save disk space.
        comp = bz2.BZ2Compressor()
        data = comp.compress(raw_data)
        data += comp.flush()
        out_handle = open(self._name_to_file(filename), 'wb')
        out_handle.write(data)
        out_handle.close()

        # Return the full path to the file.
        return out_handle.name

    def _calculate_hash(self, content):
        """
        Calculate the md5sum of a piece of text.

        :arg unicode content: Arbitrary text.

        :returns: The md5sum of 'content'.
        :rtype: unicode
        """
        hash_func = hashlib.md5()
        hash_func.update(content)
        md5sum = hash_func.hexdigest()

        return unicode(md5sum)

    def _new_ud(self):
        """
        Make a new UD number based on the current time (seconds since 1970).

        :returns: A new UD number.
        :rtype: unicode
        """
        ud = util.generate_id()
        return 'UD_' + unicode(ud)

    def _update_db_md5(self, raw_data, name, gi):
        """
        :arg str raw_data:
        :arg unicode name:
        :arg unicode gi:

        :returns: filename
        :rtype: unicode
        """
        # TODO: Documentation.
        try:
            reference = Reference.query.filter_by(accession=name).one()
            current_md5sum = reference.checksum
        except NoResultFound:
            current_md5sum = None

        if current_md5sum:
            md5sum = self._calculate_hash(raw_data)
            if md5sum != current_md5sum:
                self._output.addMessage(
                    __file__, -1, 'WHASH',
                    'Warning: Hash of {} changed from {} to {}.'.format(
                        name, current_md5sum, md5sum))
                Reference.query.filter_by(accession=name).update(
                    {'checksum': md5sum})
                session.commit()
        else:
            reference = Reference(
                name, self._calculate_hash(raw_data), geninfo_identifier=gi)
            session.add(reference)
            session.commit()
        return self._name_to_file(name)

    def snpConvert(self, rs_id):
        """
        Search for an rsId in dbSNP and return all annotated HGVS notations of
        it.

        :arg unicode rsId: The rsId of the SNP (example: 'rs9919552').

        :returns: A list of HGVS notations.
        :rtype: list(unicode)
        """
        # A simple input check.
        id = rs_id[2:]
        if rs_id[:2] != 'rs' or not id.isdigit():
            self._output.addMessage(
                __file__, 4, 'ESNPID', 'This is not a valid dbSNP id.')
            return []

        # Query dbSNP for the SNP. The following weird construct is to catch
        # any glitches in our Entrez connections. We try up to ENTREZ_MAX_TRIES
        # and only then give up.
        # Todo: maybe also implement this for other Entrez queries?
        for i in range(ENTREZ_MAX_TRIES - 1):
            try:
                response = Entrez.efetch(
                    db='snp', id=id, rettype='flt', retmode='xml')
                break
            except (IOError, HTTPException):
                time.sleep(ENTREZ_SLEEP)
        else:
            try:
                response = Entrez.efetch(
                    db='snp', id=id, rettype='flt', retmode='xml')
            except (IOError, HTTPException) as e:
                # Could not parse XML.
                self._output.addMessage(
                    __file__, 4, 'EENTREZ', 'Error connecting to dbSNP.')
                self._output.addMessage(
                    __file__, -1, 'INFO', 'IOError: {}'.format(unicode(e)))
                return []

        try:
            response_text = response.read()
        except IncompleteRead as e:
            self._output.addMessage(
                __file__, 4, 'EENTREZ', 'Error reading from dbSNP.')
            self._output.addMessage(
                __file__, -1, 'INFO', 'IncompleteRead: {}'.format(unicode(e)))
            return []

        if response_text.strip() == b'\n':
            # This is apparently what dbSNP returns for non-existing dbSNP id
            self._output.addMessage(
                __file__, 4, 'EENTREZ',
                'ID rs{} could not be found in dbSNP.'.format(id))
            return []

        try:
            # Parse the output.
            doc = minidom.parseString(response_text)
            rs = doc.getElementsByTagName('Rs')[0]
        except expat.ExpatError as e:
            # Could not parse XML.
            self._output.addMessage(
                __file__, 4, 'EENTREZ',
                'Unknown dbSNP error. Error parsing result XML.')
            self._output.addMessage(
                __file__, -1, 'INFO', 'ExpatError: {}'.format(unicode(e)))
            self._output.addMessage(
                __file__, -1, 'INFO', 'Result from dbSNP: {}'.format(
                    unicode(response_text, 'utf-8')))
            return []
        except IndexError:
            # The expected root element is not present.
            self._output.addMessage(
                __file__, 4, 'EENTREZ',
                'Unknown dbSNP error. Result XML was not as expected.')
            self._output.addMessage(
                __file__, -1, 'INFO', 'Result from dbSNP: {}'.format(
                    unicode(response_text, 'utf-8')))
            return []

        snps = []
        for i in rs.getElementsByTagName('hgvs'):
            snps.append(i.lastChild.data)

        return snps


class GenBankRetriever(Retriever):
    """
    """
    def __init__(self, output):
        """
        Initialise the class.

        :arg object output: The output object.
        """
        Retriever.__init__(self, output)
        self.file_type = 'gb'
        # TODO documentation

    def write(self, raw_data, filename, extract):
        """
        Write raw data to a file. The data is parsed before writing, if a
        parse error occurs an error is returned and the function exits.
        If 'filename' is set and 'extract' is set to 0, then 'filename' is
        used for output.
        If 'extract' is set to 1, then the filename is constructed from the
        id of the GenBank record. Additionally the id and GI number are
        returned for further processing (putting them in the internal
        database).

        :arg str raw_data: The data.
        :arg unicode filename: The intended name of the file.
        :arg int extract: Flag that indicates whether to extract the record ID
            and GI number:
            - 0 ; Do not extract, use 'filename'
            - 1 ; Extract

        :returns: Depending on the value of 'extract':
            - 0 ; ('filename', None)
            - 1 ; (id, gi)
        :rtype: tuple(unicode, unicode)
        """
        if raw_data.strip() == b'Nothing has been found':
            self._output.addMessage(
                __file__, 4, 'ENORECORD', 'The record could not be retrieved.')
            return None

        # BioPython needs a file handle.
        fake_handle = io.BytesIO()
        fake_handle.write(raw_data)
        fake_handle.seek(0)
        try:
            record = SeqIO.read(fake_handle, 'genbank')
        except (ValueError, AttributeError):
            self._output.addMessage(
                __file__, 4, 'ENOPARSE', 'The file could not be parsed.')
            return None

        if type(record.seq) == UnknownSeq:
            self._output.addMessage(
                __file__, 4, 'ENOSEQ',
                'This record contains no sequence. Chromosomal or contig '
                'records should be uploaded with the GenBank uploader.')
            return None

        out_filename = filename
        gi = None
        if extract:
            out_filename = unicode(record.id)
            gi = unicode(record.annotations['gi'])
            if out_filename != filename:
                # Add the reference (incl version) to the reference output
                # This differs if the original reference lacks a version
                self._output.addOutput('reference', unicode(record.id))
                self._output.addOutput(
                    'BatchFlags',
                    ('A1', (filename, out_filename, filename+'.')))
                self._output.addMessage(
                    __file__, 2, 'WNOVER',
                    'No version number is given, using {}. Please use this '
                    'number to reduce downloading overhead.'.format(
                        unicode(record.id)))

        if not self._write(raw_data, out_filename):
            return None

        return out_filename, gi

    def fetch(self, name):
        """
        Todo: Documentation.

        Todo: A better implementation would probably use an esummary query
            first to get the length of the sequence. If this is within limits,
            use efetch with rettype=gbwithparts to download the GenBank file.
        """
        try:
            net_handle = Entrez.efetch(
                db='nuccore', id=name, rettype='gb', retmode='text')
            raw_data = net_handle.read()
            net_handle.close()
        except (IOError, urllib2.HTTPError, HTTPException) as e:
            self._output.addMessage(
                __file__, -1, 'INFO',
                'Error connecting to Entrez nuccore database: {}'.format(
                    unicode(e)))
            self._output.addMessage(
                __file__, 4, 'ERETR', 'Could not retrieve {}.'.format(name))
            return None

        # Check if the file is empty or not.
        if raw_data.strip() == b'':
            self._output.addMessage(
                __file__, 4, 'ERETR', 'Could not retrieve {}.'.format(name))
            return None

        if b'Resource temporarily unavailable' in raw_data:
            self._output.addMessage(
                __file__, 4, 'ERETR',
                'Resource temporarily unavailable from NCBI servers: '
                '{}.'.format(name))
            return None

        # This is a hack to detect constructed references, the proper way to
        # do this would be to check the data_file_division attribute of the
        # parsed GenBank file (it would be 'CON').
        if b'\nCONTIG' in raw_data:
            try:
                # Get the length in base pairs
                length = int(
                    raw_data[:raw_data.index(b' bp', 0, 500)].split()[-1])
            except (ValueError, IndexError):
                self._output.addMessage(
                    __file__, 4, 'ERETR', 'Could not retrieve {}.'.format(
                        name))
                return None
            if length > settings.MAX_FILE_SIZE:
                self._output.addMessage(
                    __file__, 4, 'ERETR',
                    'Could not retrieve {} (exceeds maximum file size of {} '
                    'megabytes).'.format(
                        name, settings.MAX_FILE_SIZE // 1048576))
                return None
            try:
                net_handle = Entrez.efetch(
                    db='nuccore', id=name, rettype='gbwithparts',
                    retmode='text')
                raw_data = net_handle.read()
                net_handle.close()
            except (IOError, urllib2.HTTPError, HTTPException) as e:
                self._output.addMessage(
                    __file__, -1, 'INFO',
                    'Error connecting to Entrez nuccore database: {}'.format(
                        unicode(e)))
                self._output.addMessage(
                    __file__, 4, 'ERETR', 'Could not retrieve {}.'.format(
                        name))
                return None

        result = self.write(raw_data, name, 1)
        if not result:
            return None
        name, gi = result

        if name:
            # Processing went okay.
            return self._update_db_md5(raw_data, name, gi)
        else:
            # Parse error in the GenBank file.
            return None

    def retrieveslice(self, accno, start, stop, orientation):
        """
        Retrieve a slice of a chromosome.
        If the arguments are recognised (found in the internal database),
        we look if the associated file is still present and if so: return
        its UD number.
        If the arguments are recognised but no file was found, we download
        the new slice and update the hash (and log if the hash changes).
        If the arguments are not recognised, we download the new slice and
        make a new UD number.
        The content of the slice is placed in the cache with the UD number
        as filename.

        :arg unicode accno: The accession number of the chromosome.
        :arg int start: Start position of the slice (one-based, inclusive, in
          reference orientation).
        :arg int stop: End position of the slice (one-based, inclusive, in
          reference orientation).
        :arg int orientation: Orientation of the slice:
            - 1 ; Forward.
            - 2 ; Reverse complement.

        :returns: An UD number.
        :rtype: unicode
        """
        # Not a valid slice.
        if start > stop:
            self._output.addMessage(__file__, 4, 'ERETR',
                                    'Could not retrieve slice for start '
                                    'position greater than stop position.')
            return None

        # The slice can not be too big.
        if stop - start + 1 > settings.MAX_FILE_SIZE:
            self._output.addMessage(__file__, 4, 'ERETR',
                                    'Could not retrieve slice (request '
                                    'exceeds maximum of %d bases)' %
                                    settings.MAX_FILE_SIZE)
            return None

        slice_orientation = ['forward', 'reverse'][orientation - 1]

        # Check whether we have seen this slice before.
        try:
            reference = Reference.query.filter_by(
                slice_accession=accno, slice_start=start, slice_stop=stop,
                slice_orientation=slice_orientation).one()
        except NoResultFound:
            reference = None
        else:
            if os.path.isfile(self._name_to_file(reference.accession)):
                # It's still present.
                return reference.accession

        # It's not present, so download it.
        try:
            # EFetch `seq_start` and `seq_stop` are one-based, inclusive, and
            # in reference orientation.
            handle = Entrez.efetch(
                db='nuccore', rettype='gb', retmode='text', id=accno,
                seq_start=start, seq_stop=stop, strand=orientation)
            raw_data = handle.read()
            handle.close()
        except (IOError, urllib2.HTTPError, HTTPException) as e:
            self._output.addMessage(
                __file__, -1, 'INFO',
                'Error connecting to Entrez nuccore database: {}'.format(
                    unicode(e)))
            self._output.addMessage(
                __file__, 4, 'ERETR', 'Could not retrieve slice.')
            return None

        # Calculate the hash of the downloaded file.
        md5sum = self._calculate_hash(raw_data)

        if reference is not None:
            # We have seen this one before.
            current_md5sum = reference.checksum

            if md5sum != current_md5sum:
                self._output.addMessage(
                    __file__, -1, 'WHASH',
                    'Warning: Hash of {} changed from {} to {}.'.format(
                        reference.accession, current_md5sum, md5sum))
                Reference.query.filter_by(
                    accession=reference.accession).update({'checksum': md5sum})
                session.commit()
        else:
            # We haven't seen it before, so give it a name.
            ud = self._new_ud()
            slice_orientation = ['forward', 'reverse'][orientation - 1]
            reference = Reference(
                ud, md5sum, slice_accession=accno, slice_start=start,
                slice_stop=stop, slice_orientation=slice_orientation)
            session.add(reference)
            session.commit()

        if self.write(raw_data, reference.accession, 0):
            return reference.accession

    def retrievegene(self, gene, organism, upstream=0, downstream=0):
        """
        Query the NCBI for the chromosomal location of a gene and make a
        slice if the gene can be found.

        :arg unicode gene: Name of the gene.
        :arg unicode organism: The organism in which we search.
        :arg int upstream: Number of upstream nucleotides for the slice.
        :arg int downstream: Number of downstream nucleotides for the slice.

        :returns: GenBank record.
        :rtype: object
        """
        # Search the NCBI for a specific gene in an organism.
        query = '{}[Gene] AND {}[Orgn]'.format(gene, organism)
        try:
            handle = Entrez.esearch(db='gene', term=query)
            try:
                search_result = Entrez.read(handle)
            except Entrez.Parser.ValidationError:
                self._output.addMessage(
                    __file__, -1, 'INFO',
                    'Error reading Entrez esearch result.')
                self._output.addMessage(
                    __file__, 4, 'ERETR',
                    'Could not search for gene {}.'.format(gene))
                return None
            finally:
                handle.close()
        except (IOError, urllib2.HTTPError, HTTPException) as e:
            self._output.addMessage(
                __file__, -1, 'INFO',
                'Error connecting to Entrez esearch: {}'.format(unicode(e)))
            self._output.addMessage(
                __file__, 4, 'ERETR',
                'Could not search for gene {}: {}'.format(gene, unicode(e)))
            return None

        chr_acc_ver = None  # We did not find anything yet.
        aliases = []  # A list of aliases in case we find them.

        for i in search_result['IdList']:
            # Inspect all results.
            try:
                handle = Entrez.esummary(db='gene', id=i)
                try:
                    summary = Entrez.read(handle)
                except Entrez.Parser.ValidationError:
                    self._output.addMessage(
                        __file__, -1, 'INFO',
                        'Error reading Entrez esummary result.')
                    self._output.addMessage(
                        __file__, 4, 'ERETR',
                        'Could not get mapping information for gene '
                        '{}.'.format(gene))
                    return None
                finally:
                    handle.close()
            except (IOError, urllib2.HTTPError, HTTPException) as e:
                self._output.addMessage(
                    __file__, -1, 'INFO',
                    'Error connecting to Entrez esummary: {}'.format(
                        unicode(e)))
                self._output.addMessage(
                    __file__, 4, 'ERETR',
                    'Could not get mapping information for gene {}: {}'.format(
                        gene, unicode(e)))
                return None

            try:
                document = summary['DocumentSummarySet']['DocumentSummary'][0]
            except (KeyError, IndexError):
                self._output.addMessage(
                    __file__, -1, 'INFO',
                    'Error parsing Entrez esummary result.')
                self._output.addMessage(
                    __file__, 4, 'ERETR',
                    'Could not get mapping information for gene {}.'.format(
                        gene))
                return None

            # For the available fields and their meaning, see Table 9 in:
            # http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.359.4838&rep=rep1&type=pdf
            if gene.lower() in [unicode(document[key]).lower()
                                for key in ('NomenclatureSymbol', 'Name')]:
                # Found it.
                if not document['GenomicInfo']:
                    self._output.addMessage(
                        __file__, 4, 'ENOMAPPING',
                        'No mapping information found for gene {}.'.format(
                            gene))
                    return None

                # Positions are zero-based, inclusive and in gene orientation.
                chr_acc_ver = unicode(document['GenomicInfo'][0]['ChrAccVer'])
                chr_start = int(document['GenomicInfo'][0]['ChrStart'])
                chr_stop = int(document['GenomicInfo'][0]['ChrStop'])

                # Convert to one-based, inclusive, reference orientation. We
                # also add the flanking regions.
                chr_start += 1
                chr_stop += 1
                if chr_start <= chr_stop:
                    # Gene on forward strand.
                    orientation = 1
                    chr_start -= upstream
                    chr_stop += downstream
                else:
                    # Gene on reverse strand.
                    orientation = 2
                    chr_start, chr_stop = chr_stop, chr_start
                    chr_start -= downstream
                    chr_stop += upstream

                break

            # Collect official symbols that has this gene as alias in case we
            # can not find anything.
            if (gene in unicode(document['OtherAliases']).split(',') and
                (document['NomenclatureSymbol'] or document['Name'])):
                aliases.append(unicode(document['NomenclatureSymbol'] or document['Name']))

        if not chr_acc_ver:
            # We did not find any genes.
            if aliases:
                self._output.addMessage(
                    __file__, 4, 'ENOGENE',
                    'Gene {} not found, found aliases: {}'.format(
                        gene, ', '.join(aliases)))
                return None
            self._output.addMessage(
                __file__, 4, 'ENOGENE', 'Gene {} not found.'.format(gene))
            return None

        # And retrieve the slice.
        return self.retrieveslice(
            chr_acc_ver, chr_start, chr_stop, orientation)

    def downloadrecord(self, url):
        """
        Download a GenBank record from a URL.
        If the downloaded file is recognised by its hash, the old UD number
        is used.

        :arg unicode url: Location of a GenBank record.

        :returns: UD or None.
        :rtype: unicode
        """
        if not (url.startswith('http://') or url.startswith('https://') or
                url.startswith('ftp://')):
            self._output.addMessage(
                __file__, 4, 'ERECPARSE',
                'Only HTTP(S) or FTP locations are allowed.')
            return None

        handle = urllib2.urlopen(url)
        info = handle.info()
        if info.gettype() == 'text/plain':
            length = int(info['Content-Length'])
            if 512 < length < settings.MAX_FILE_SIZE:
                raw_data = handle.read()
                md5sum = self._calculate_hash(raw_data)

                ud = None
                try:
                    reference = Reference.query.filter_by(
                        checksum=md5sum).one()
                except NoResultFound:
                    ud = self._new_ud()
                    if not os.path.isfile(self._name_to_file(ud)):
                        ud = self.write(raw_data, ud, 0) and ud
                    if ud:
                        # Parsing went OK, add to DB.
                        reference = Reference(ud, md5sum, download_url=url)
                        session.add(reference)
                        session.commit()
                else:
                    if (os.path.isfile(self._name_to_file(reference.accession)) or
                            self.write(raw_data, reference.accession, 0)):
                        ud = reference.accession

                # Returns the UD or None.
                return ud
            else:
                self._output.addMessage(
                    __file__, 4, 'EFILESIZE',
                    'Filesize is not within the allowed boundaries.')
                return None
        else:
            self._output.addMessage(
                __file__, 4, 'ERECPARSE', 'This is not a GenBank record.')
            return None

    def uploadrecord(self, raw_data):
        """
        Write an uploaded record to a file.
        If the downloaded file is recognised by its hash, the old UD number
        is used.

        :arg str raw_data: A GenBank record.

        :returns: Accession number for the uploaded file.
        :rtype: unicode
        """
        md5sum = self._calculate_hash(raw_data)

        try:
            reference = Reference.query.filter_by(checksum=md5sum).one()
        except NoResultFound:
            ud = self._new_ud()
            if self.write(raw_data, ud, 0):
                reference = Reference(ud, md5sum)
                session.add(reference)
                session.commit()
                return ud
        else:
            if os.path.isfile(self._name_to_file(reference.accession)):
                return reference.accession
            else:
                return (self.write(raw_data, reference.accession, 0) and
                        reference.accession)

    def loadrecord(self, identifier):
        """
        Load a RefSeq record and return it.

        The record is found by trying the following options in order:

        1. Returned from the cache if it is there.
        2. Re-created (if it was created by slicing) or re-downloaded (if it
           was created by URL) if we have information on its source in the
           database.
        3. Fetched from the NCBI.

        :arg unicode identifier: A RefSeq accession number or geninfo
            identifier (GI).

        :returns: A parsed RefSeq record or `None` if no record could be found
          for the given identifier.
        :rtype: object
        """
        if identifier[0].isdigit():
            # This is a GI number (geninfo identifier).
            reference = Reference.query \
                .filter_by(geninfo_identifier=identifier) \
                .first()
        else:
            # This is a RefSeq accession number.
            reference = Reference.query \
                .filter_by(accession=identifier) \
                .first()

        if reference is None:
            # We don't know it, fetch it from NCBI.
            filename = self.fetch(identifier)

        else:
            # We have seen it before.
            filename = self._name_to_file(reference.accession)

            if os.path.isfile(filename):
                # It is still in the cache, so filename is valid.
                pass

            elif reference.slice_accession:
                # It was previously created by slicing.
                cast_orientation = {
                    None: None, 'forward': 1, 'reverse': 2}
                if not self.retrieveslice(
                        reference.slice_accession, reference.slice_start,
                        reference.slice_stop,
                        cast_orientation[reference.slice_orientation]):
                    filename = None

            elif reference.download_url:
                # It was previously created by URL.
                if not self.downloadrecord(reference.download_url):
                    filename = None

            elif reference.geninfo_identifier:
                # It was previously fetched from NCBI.
                filename = self.fetch(reference.accession)

            else:
                # It was previously created by uploading.
                self._output.addMessage(
                    __file__, 4, 'ERETR', 'Please upload this sequence again.')
                filename = None

        # If filename is None, we could not retrieve the record.
        if filename is None:
            # Notify batch job to skip all instance of identifier.
            self._output.addOutput('BatchFlags', ('S1', identifier))
            return None

        # Now we have the file, so we can parse it.
        genbank_parser = genbank.GBparser()
        record = genbank_parser.create_record(filename)

        if reference:
            record.id = reference.accession
        else:
            record.id = record.source_id

        # Todo: This will change once we support protein references.
        if isinstance(record.seq.alphabet, ProteinAlphabet):
            self._output.addMessage(
                __file__, 4, 'ENOTIMPLEMENTED',
                'Protein reference sequences are not supported.')
            return None

        return record


class LRGRetriever(Retriever):
    """
    Retrieve a LRG record from either the cache or the web.
    """
    def __init__(self, output):
        """
        Initialize the class.

        :arg object output: The output object.
        """
        Retriever.__init__(self, output)
        self.file_type = 'xml'

    def loadrecord(self, identifier):
        """
        Load and parse a LRG file based on the identifier.

        :arg unicode identifier: The name of the LRG file to read.

        :returns: GenRecord.Record of LRG file or None in case of failure.
        :rtype: object
        """
        # Make a filename based upon the identifier.
        filename = self._name_to_file(identifier)

        if not os.path.isfile(filename):
            # We can't find the file.
            filename = self.fetch(identifier)

        if filename is None:
            # Notify batch to skip all instance of identifier.
            self._output.addOutput('BatchFlags', ('S1', identifier))
            # return None in case of error.
            return None

        # Now we have the file, so we can parse it.
        file_handle = bz2.BZ2File(filename, 'r')

        # Create GenRecord.Record from LRG file.
        record = lrg.create_record(file_handle.read())
        file_handle.close()

        # We don't create LRGs from other sources, so id is always the same
        # as source_id.
        record.id = identifier
        record.source_id = identifier

        return record

    def fetch(self, name):
        """
        Fetch the LRG file and store in the cache directory. First try to
        grab the file from the confirmed section, if this fails, get it
        from the pending section.

        :arg unicode name: The name of the LRG file to fetch.

        :returns: the full path to the file; None in case of an error.
        :rtype: unicode
        """
        prefix = settings.LRG_PREFIX_URL
        url = prefix + '{}.xml'.format(name)
        pending_url = prefix + 'pending/{}.xml'.format(name)

        try:
            return self.downloadrecord(url, name)
        except urllib2.URLError:
            # Catch error: file not found.
            pass

        try:
            # Try to get the file from the pending section.
            filename = self.downloadrecord(pending_url, name)
            self._output.addMessage(
                __file__, 2, 'WPEND',
                'Warning: LRG file {} is a pending entry.'.format(name))
            return filename
        except urllib2.URLError:
            self._output.addMessage(
                __file__, 4, 'ERETR', 'Could not retrieve {}.'.format(name))
            # Explicit return in case of an Error.
            return None

    def downloadrecord(self, url, name=None):
        """
        Download an LRG record from an URL.

        :arg unicode url: Location of the LRG record.

        :returns: The full path to the file or Nonein case of failure.
        :rtype: unicode
        """
        lrg_id = name or os.path.splitext(os.path.split(url)[1])[0]
        # if not lrg_id.startswith('LRG'):
        #     return None
        filename = self._name_to_file(lrg_id)

        # TODO: Properly read the file contents to a unicode string and write
        # it utf-8 encoded.
        handle = urllib2.urlopen(url)
        info = handle.info()

        if (info['Content-Type'] == 'application/xml' and
                'Content-length' in info):
            # Looks like a valid LRG file.

            length = int(info['Content-Length'])
            if 512 < length < settings.MAX_FILE_SIZE:
                raw_data = handle.read()
                handle.close()

                # Do an md5 check.
                md5sum = self._calculate_hash(raw_data)
                try:
                    reference = Reference.query.filter_by(
                        accession=lrg_id).one()
                    md5_db = reference.checksum
                except NoResultFound:
                    md5_db = None

                if md5_db is None:
                    reference = Reference(lrg_id, md5sum, download_url=url)
                    session.add(reference)
                    session.commit()
                elif md5sum != md5_db:
                    # Hash has changed for the LRG ID.
                    self._output.addMessage(
                        __file__, -1, 'WHASH',
                        'Warning: Hash of {} changed from {} to {}.'.format(
                            lrg_id, md5_db, md5sum))
                    Reference.query.filter_by(accession=lrg_id).update(
                        {'checksum': md5sum})
                    session.commit()
                else:
                    # Hash the same as in db.
                    pass

                if not os.path.isfile(filename):
                    return self.write(raw_data, lrg_id)
                else:
                    # This can only occur if synchronus calls to mutalyzer are
                    # made to recover a file that did not exist. Still leaves
                    # a window in between the check and the write.
                    return filename
            else:
                self._output.addMessage(
                    __file__, 4, 'EFILESIZE',
                    'Filesize is not within the allowed boundaries.')
        else:
            self._output.addMessage(
                __file__, 4, 'ERECPARSE', 'This is not an LRG record.')
        handle.close()

    def write(self, raw_data, filename):
        """
        Write raw LRG data to a file. The data is parsed before writing,
        if a parse error occurs None is returned.

        :arg str raw_data: The data.
        :arg unicode filename: The intended name of the file.

        :returns: The full path and name of the file written, None in case of
          an error.
        :rtype: unicode
        """
        # Dirty way to test if a file is valid,
        # Parse the file to see if it's a real LRG file.
        try:
            lrg.create_record(raw_data)
        except DOMException:
            self._output.addMessage(
                __file__, 4, 'ERECPARSE', 'Could not parse file.')
            # Explicit return on Error.
            return None

        # Returns full path.
        return self._write(raw_data, filename)
