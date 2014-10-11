"""
Module for retrieving files from either the cache or the NCBI.

A hash of every retrieved file is stored in the internal database. If a
requested file is not found, but its hash is, we use additional information
to re-download the file.

Public classes:
- Retriever ; Retrieve a record from either the cache or the NCBI.
"""


from __future__ import unicode_literals

import codecs
import os              # path.isfile(), link() path.isdir(), path.mkdir(),
                       # walk(), path.getsize(), path.join(), stat(), remove()
import time
import bz2             # BZ2Compressor(), BZ2File()
import hashlib         # md5(), update(), hexdigest()
import urllib2         # urlopen()
import StringIO        # StringIO()
from Bio import SeqIO  # read()
from Bio import Entrez # efetch(), read(), esearch(), esummary()
from Bio.Seq import UnknownSeq
from Bio.Alphabet import ProteinAlphabet
from xml.dom import DOMException, minidom
from xml.parsers import expat
from httplib import HTTPException, IncompleteRead
from sqlalchemy.orm.exc import NoResultFound

from mutalyzer import util
from mutalyzer.config import settings
from mutalyzer.db import session
from mutalyzer.db.models import Reference
from mutalyzer.parsers import lrg
from mutalyzer.parsers import genbank


ENTREZ_MAX_TRIES = 4
ENTREZ_SLEEP = 1      # in seconds


class Retriever(object) :
    """
    Retrieve a record from either the cache or the NCBI.

    Special methods:
        - __init__(output, database) ; Use variables from the
        configuration file to initialise the class private variables.

    Private methods:
        - _nametofile(name)   ; Convert a name to a filename.
        - _write(raw_data, filename, extract) ; Write a record to a file.
        - _calcHash(content)  ; Calculate the md5sum of 'content'.
        - _newUD()            ; Generate a new UD number.

    Public methods:
        - retrieveslice(accno, start, stop, orientation) ; Retrieve a chromosome
        slice from the NCBI.
        - retrievegene(gene, organism, upstream, downstream) ; Retrieve a gene
        from the NCBI.
        - downloadrecord(url)    ; Download a GenBank file.
        - uploadrecord(raw_data) ; Let someone upload a GenBank file.
        - loadrecord(identifier) ; Load a record, store it in the cache, manage
        the cache and return the record.
    """
    def __init__(self, output) :
        """
        Use variables from the configuration file for some simple
        settings. Make the cache directory if it does not exist yet.

        @arg output:
        @type output:
        @arg database:
        @type database:
        """
        self._output = output
        if not os.path.isdir(settings.CACHE_DIR) :
            os.mkdir(settings.CACHE_DIR)
        Entrez.email = settings.EMAIL
        self.fileType = None
    #__init__

    def _nametofile(self, name) :
        """
        Convert an accession number to a filename.

        @arg name: The accession number
        @type name: unicode

        @return: A filename
        @rtype: unicode
        """
        return os.path.join(settings.CACHE_DIR, name + "." + self.fileType + ".bz2")
    #_nametofile

    def _write(self, raw_data, filename) :
        """
        Write raw data to a compressed file.

        @arg raw_data: The raw_data to be compressed and written
        @type raw_data: string
        @arg filename: The intended name of the outfile
        @type filename: unicode

        @return: outfile ; The full path and name of the file written
        @rtype: unicode
        """
        # Todo: Should we write a utf-8 encoded genbank file? Not even sure
        #   what type `raw_data` is...
        # Compress the data to save disk space.
        comp = bz2.BZ2Compressor()
        data = comp.compress(raw_data)
        data += comp.flush()
        out_handle = open(self._nametofile(filename), "w")
        out_handle.write(data)
        out_handle.close()

        return out_handle.name      # return the full path to the file
    #_write

    # Todo: check callers; argument should be a byte string
    def _calcHash(self, content) :
        """
        Calculate the md5sum of a piece of text.

        @arg content: Arbitrary text
        @type content: byte string

        @return: The md5sum of 'content'
        @rtype: unicode
        """

        hashfunc = hashlib.md5()
        hashfunc.update(content)
        md5sum = hashfunc.hexdigest()
        del hashfunc

        return unicode(md5sum)
    #_calcHash

    def _newUD(self) :
        """
        Make a new UD number based on the current time (seconds since 1970).

        @return: A new UD number
        @rtype: unicode
        """

        UD = util.generate_id()
        return "UD_" + unicode(UD)
    #_newUD

    def _updateDBmd5(self, raw_data, name, GI):
        #TODO documentation
        """
        @todo: documentation

        @arg raw_data:
        @type raw_data:
        @arg name:
        @type name:
        @arg GI:
        @type GI:

        @return: filename
        @rtype: unicode
        """
        try:
            reference = Reference.query.filter_by(accession=name).one()
            currentmd5sum = reference.checksum
        except NoResultFound:
            currentmd5sum = None

        if currentmd5sum :
            md5sum = self._calcHash(raw_data)
            if md5sum != currentmd5sum :
                self._output.addMessage(__file__, -1, "WHASH",
                    "Warning: Hash of %s changed from %s to %s." % (
                    name, currentmd5sum, md5sum))
                Reference.query.filter_by(accession=name).update({'checksum': md5sum})
                session.commit()
            #if
        else :
            reference = Reference(name, self._calcHash(raw_data),
                                  geninfo_identifier=GI)
            session.add(reference)
            session.commit()
        return self._nametofile(name)
    #_updateDBmd5


    def snpConvert(self, rs_id) :
        """
        Search for an rsId in dbSNP and return all annotated HGVS notations of
        it.

        @arg rsId: The rsId of the SNP (example: 'rs9919552').
        @type rsId: unicode

        @return: A list of HGVS notations.
        @rtype: list(unicode)
        """
        # A simple input check.
        id = rs_id[2:]
        if rs_id[:2] != 'rs' or not id.isdigit():
            self._output.addMessage(__file__, 4, 'ESNPID',
                                    'This is not a valid dbSNP id.')
            return []

        # Query dbSNP for the SNP. The following weird construct is to catch
        # any glitches in our Entrez connections. We try up to ENTREZ_MAX_TRIES
        # and only then give up.
        # Todo: maybe also implement this for other Entrez queries?
        for i in range(ENTREZ_MAX_TRIES - 1):
            try:
                response = Entrez.efetch(db='snp', id=id, rettype='flt',
                                         retmode='xml')
                break
            except (IOError, HTTPException):
                time.sleep(ENTREZ_SLEEP)
        else:
            try:
                response = Entrez.efetch(db='snp', id=id, rettype='flt',
                                         retmode='xml')
            except (IOError, HTTPException) as e:
                # Could not parse XML.
                self._output.addMessage(__file__, 4, 'EENTREZ',
                                        'Error connecting to dbSNP.')
                self._output.addMessage(__file__, -1, 'INFO',
                                        'IOError: %s' % unicode(e))
                return []

        try:
            response_text = response.read()
        except IncompleteRead as e:
            self._output.addMessage(__file__, 4, 'EENTREZ',
                                    'Error reading from dbSNP.')
            self._output.addMessage(__file__, -1, 'INFO',
                                    'IncompleteRead: %s' % unicode(e))
            return []

        if response_text == '\n':
            # This is apparently what dbSNP returns for non-existing dbSNP id
            self._output.addMessage(__file__, 4, 'EENTREZ',
                                    'ID rs%s could not be found in dbSNP.' \
                                    % id)
            return []

        try:
            # Parse the output.
            doc = minidom.parseString(response_text)
            rs = doc.getElementsByTagName('Rs')[0]
        except expat.ExpatError as e:
            # Could not parse XML.
            self._output.addMessage(__file__, 4, 'EENTREZ', 'Unknown dbSNP ' \
                                    'error. Error parsing result XML.')
            self._output.addMessage(__file__, -1, 'INFO',
                                    'ExpatError: %s' % unicode(e))
            self._output.addMessage(__file__, -1, 'INFO',
                                    'Result from dbSNP: %s' % response_text)
            return []
        except IndexError:
            # The expected root element is not present.
            self._output.addMessage(__file__, 4, 'EENTREZ', 'Unknown dbSNP ' \
                                    'error. Result XML was not as expected.')
            self._output.addMessage(__file__, -1, 'INFO',
                                    'Result from dbSNP: %s' % response_text)
            return []

        snps = []
        for i in rs.getElementsByTagName('hgvs'):
            snps.append(i.lastChild.data)

        return snps
    #snpConvert
#Retriever

class GenBankRetriever(Retriever):
    # TODO documentation
    """
    """

    def __init__(self, output):
        """
        @todo: Documentation.
        """
        # Recall init of parent
        Retriever.__init__(self, output)
        self.fileType = "gb"
        # Child specific init
    #__init__

    # todo: raw_data must always be a byte string
    def write(self, raw_data, filename, extract) :
        """
        Write raw data to a file. The data is parsed before writing, if a
        parse error occurs an error is returned and the function exits.
        If 'filename' is set and 'extract' is set to 0, then 'filename' is
        used for output.
        If 'extract' is set to 1, then the filename is constructed from the
        id of the GenBank record. Additionally the id and GI number are
        returned for further processing (putting them in the internal
        database).

        @arg raw_data: The data
        @type raw_data: string
        @arg filename: The intended name of the file.
        @type filename: unicode
        @arg extract: Flag that indicates whether to extract the record ID and
        GI number:
            - 0 ; Do not extract, use 'filename'
            - 1 ; Extract
        @type extract: integer

        @return: tuple ; Depending on the value of 'extract':
            - 0 ; ('filename', None)
            - 1 ; (id, GI)
        @rtype: tuple (unicode, unicode)
        """

        if raw_data == "\nNothing has been found\n" :
            self._output.addMessage(__file__, 4, "ENORECORD",
                "The record could not be retrieved.")
            return None
        #if

        fakehandle = StringIO.StringIO() # Unfortunately, BioPython needs a
        fakehandle.write(raw_data)       # file handle.
        fakehandle.seek(0)
        try :
            record = SeqIO.read(fakehandle, "genbank")
        except (ValueError, AttributeError):  # An error occured while parsing.
            self._output.addMessage(__file__, 4, "ENOPARSE",
                "The file could not be parsed.")
            fakehandle.close()
            return None
        #except

        if type(record.seq) == UnknownSeq :
            fakehandle.close()
            self._output.addMessage(__file__, 4, "ENOSEQ",
                "This record contains no sequence. Chromosomal or contig " \
                "records should be uploaded with the GenBank uploader.")
            return None
        #if

        outfile = filename
        GI = None
        if extract :
            outfile = record.id
            GI = record.annotations["gi"]
            if outfile != filename :
                # Add the reference (incl version) to the reference output
                # This differs if the original reference lacks a version
                self._output.addOutput("reference", record.id)
                self._output.addOutput(
                        "BatchFlags", ("A1",(
                            filename,
                            outfile,
                            filename+"." )))
                self._output.addMessage(__file__, 2, "WNOVER",
                    "No version number is given, using %s. Please use this " \
                    "number to reduce downloading overhead." % record.id)
        #if
        fakehandle.close()

        self._write(raw_data, outfile)

        return outfile, GI
    #write

    def fetch(self, name) :
        """
        Todo: Documentation.

        Todo: A better implementation would probably use an esummary query
            first to get the length of the sequence. If this is within limits,
            use efetch with rettype=gbwithparts to download the GenBank file.
        """
        try:
            net_handle = Entrez.efetch(db='nuccore', id=name, rettype='gb', retmode='text')
            raw_data = net_handle.read()
            net_handle.close()
        except (IOError, urllib2.HTTPError, HTTPException) as e:
            self._output.addMessage(__file__, -1, 'INFO',
                                    'Error connecting to Entrez nuccore database: %s' % unicode(e))
            self._output.addMessage(__file__, 4, 'ERETR',
                                    'Could not retrieve %s.' % name)
            return None

        if raw_data == '\n' :       # Check if the file is empty or not.
            self._output.addMessage(__file__, 4, 'ERETR',
                                    'Could not retrieve %s.' % name)
            return None

        # This is a hack to detect constructed references, the proper way to
        # do this would be to check the data_file_division attribute of the
        # parsed GenBank file (it would be 'CON').
        if '\nCONTIG' in raw_data:
            try:
                # Get the length in base pairs
                length = int(raw_data[:raw_data.index(' bp', 0, 500)].split()[-1])
            except ValueError, IndexError:
                self._output.addMessage(__file__, 4, 'ERETR',
                                        'Could not retrieve %s.' % name)
                return None
            if length > settings.MAX_FILE_SIZE:
                self._output.addMessage(__file__, 4, 'ERETR',
                                        'Could not retrieve %s.' % name)
                return None
            try:
                net_handle = Entrez.efetch(db='nuccore', id=name, rettype='gbwithparts', retmode='text')
                raw_data = net_handle.read()
                net_handle.close()
            except (IOError, urllib2.HTTPError, HTTPException) as e:
                self._output.addMessage(__file__, -1, 'INFO',
                                        'Error connecting to Entrez nuccore database: %s' % unicode(e))
                self._output.addMessage(__file__, 4, 'ERETR',
                                        'Could not retrieve %s.' % name)
                return None

        result = self.write(raw_data, name, 1)
        if not result:
            return None
        name, GI = result
        if name:               # Processing went okay.
            return self._updateDBmd5(raw_data, name, GI)
        else:                  # Parse error in the GenBank file.
            return None
    #fetch

    def retrieveslice(self, accno, start, stop, orientation) :
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

        @arg accno: The accession number of the chromosome
        @type accno: unicode
        @arg start: Start position of the slice
        @type start: integer
        @arg stop: End position of the slice.
        @type stop: integer
        @arg orientation:
        Orientation of the slice:
            - 1 ; Forward
            - 2 ; Reverse complement
        @type orientation: integer

        @return: An UD number
        @rtype: unicode
        """

        # Not a valid slice.
        if start >= stop :
            return None

        # The slice can not be too big.
        if stop - start > settings.MAX_FILE_SIZE:
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
            if os.path.isfile(self._nametofile(reference.accession)) : # It's still present.
                return reference.accession

        # It's not present, so download it.
        try:
            handle = Entrez.efetch(db='nuccore', rettype='gb', retmode='text',
                                   id=accno, seq_start=start, seq_stop=stop,
                                   strand=orientation)
            raw_data = handle.read()
            handle.close()
        except (IOError, urllib2.HTTPError, HTTPException) as e:
            self._output.addMessage(__file__, -1, 'INFO',
                                    'Error connecting to Entrez nuccore database: %s' % unicode(e))
            self._output.addMessage(__file__, 4, 'ERETR',
                                    'Could not retrieve slice.')
            return None

        # Calculate the hash of the downloaded file.
        md5sum = self._calcHash(raw_data)

        if reference is not None: # We have seen this one before.
            currentmd5sum = reference.checksum

            if md5sum != currentmd5sum :
                self._output.addMessage(__file__, -1, "WHASH",
                    "Warning: Hash of %s changed from %s to %s." % (
                    reference.accession, currentmd5sum, md5sum))
                Reference.query.filter_by(accession=reference.accession).update({'checksum': md5sum})
                session.commit()
            #if
        else : # We haven't seen it before, so give it a name.
            UD = self._newUD()
            slice_orientation = ['forward', 'reverse'][orientation - 1]
            reference = Reference(UD, md5sum, slice_accession=accno,
                                  slice_start=start, slice_stop=stop,
                                  slice_orientation=slice_orientation)
            session.add(reference)
            session.commit()
        #else

        if self.write(raw_data, reference.accession, 0):
            return reference.accession
    #retrieveslice

    def retrievegene(self, gene, organism, upstream, downstream) :
        """
        Query the NCBI for the chromosomal location of a gene and make a
        slice if the gene can be found.

        @arg gene: Name of the gene
        @type gene: unicode
        @arg organism: The organism in which we search.
        @type organism: unicode
        @arg upstream: Number of upstream nucleotides for the slice.
        @type upstream: integer
        @arg downstream: Number of downstream nucleotides for the slice.
        @type downstream: integer

        @return: slice
        @rtype:
        """

        # Search the NCBI for a specific gene in an organism.
        query = "%s[Gene] AND %s[Orgn]" % (gene, organism)
        try:
            handle = Entrez.esearch(db = "gene", term = query)
            try:
                searchresult = Entrez.read(handle)
            except Entrez.Parser.ValidationError:
                self._output.addMessage(__file__, -1, 'INFO',
                                        'Error reading Entrez esearch result.')
                self._output.addMessage(__file__, 4, 'ERETR',
                                        'Could not search for gene %s.' % gene)
                return None
            finally:
                handle.close()
        except (IOError, urllib2.HTTPError, HTTPException) as e:
            self._output.addMessage(__file__, -1, 'INFO',
                                    'Error connecting to Entrez esearch: %s' % unicode(e))
            self._output.addMessage(__file__, 4, 'ERETR',
                                    'Could not search for gene %s.' % gene)
            return None

        ChrAccVer = None        # We did not find anything yet.
        aliases = []            # A list of aliases in case we find them.
        for i in searchresult["IdList"] :                 # Inspect all results.
            try:
                handle = Entrez.esummary(db = "gene", id = i)
                try:
                    summary = Entrez.read(handle)
                except Entrez.Parser.ValidationError:
                    self._output.addMessage(__file__, -1, 'INFO',
                                            'Error reading Entrez esummary result.')
                    self._output.addMessage(__file__, 4, 'ERETR',
                                            'Could not get mapping information for gene %s.' % gene)
                    return None
                finally:
                    handle.close()
            except (IOError, urllib2.HTTPError, HTTPException) as e:
                self._output.addMessage(__file__, -1, 'INFO',
                                        'Error connecting to Entrez esummary: %s' % unicode(e))
                self._output.addMessage(__file__, 4, 'ERETR',
                                        'Could not get mapping information for gene %s.' % gene)
                return None

            if summary[0]["NomenclatureSymbol"].lower() == gene.lower() : # Found it.
                if not summary[0]["GenomicInfo"] :
                    self._output.addMessage(__file__, 4, "ENOMAPPING",
                        "No mapping information found for gene %s." % gene)
                    return None
                #if
                ChrAccVer = summary[0]["GenomicInfo"][0]["ChrAccVer"]
                ChrLoc = summary[0]["GenomicInfo"][0]["ChrLoc"]
                ChrStart = summary[0]["GenomicInfo"][0]["ChrStart"]
                ChrStop = summary[0]["GenomicInfo"][0]["ChrStop"]
                break;
            #if

            # Collect official symbols that has this gene as alias in case we
            # can not find anything.
            if gene in summary[0]["OtherAliases"] and \
                summary[0]["NomenclatureSymbol"] :
                aliases.append(summary[0]["NomenclatureSymbol"]);
        #for

        if not ChrAccVer : # We did not find any genes.
            if aliases :
                self._output.addMessage(__file__, 4, "ENOGENE",
                    "Gene %s not found, found aliases: %s" % (gene, aliases))
                return None
            #if
            self._output.addMessage(__file__, 4, "ENOGENE",
                "Gene %s not found." % gene)
            return None
        #if

        # Figure out the orientation of the gene.
        orientation = 1
        if ChrStart > ChrStop :             # Swap start and stop.
            orientation = 2
            temp = ChrStart
            ChrStart = ChrStop - downstream # Also take care of the flanking
            ChrStop = temp + upstream + 1   # sequences.
        #if
        else :
            ChrStart -= upstream - 1
            ChrStop += downstream + 2
        #else

        # And retrieve the slice.
        return self.retrieveslice(ChrAccVer, ChrStart, ChrStop, orientation)
    #retrievegene

    def downloadrecord(self, url) :
        """
        Download a GenBank record from a URL.
        If the downloaded file is recognised by its hash, the old UD number
        is used.

        @arg url: Location of a GenBank record
        @type url: unicode

        @return: UD or None
        @rtype: unicode
        """
        handle = urllib2.urlopen(url)
        info = handle.info()
        if info["Content-Type"] == "text/plain" :
            length = int(info["Content-Length"])
            if 512 < length < settings.MAX_FILE_SIZE:
                raw_data = handle.read()
                md5sum = self._calcHash(raw_data)

                UD = None

                try:
                    reference = Reference.query.filter_by(checksum=md5sum).one()
                except NoResultFound:
                    UD = self._newUD()
                    if not os.path.isfile(self._nametofile(UD)):
                        UD = self.write(raw_data, UD, 0) and UD
                    if UD:      #Parsing went OK, add to DB
                        reference = Reference(UD, md5sum, download_url=url)
                        session.add(reference)
                        session.commit()
                else:
                    if not os.path.isfile(self._nametofile(reference.accession)):
                        UD = self.write(raw_data, reference.accession, 0) and reference.accession

                return UD #Returns the UD or None
            #if
            else :
                self._output.addMessage(__file__, 4, "EFILESIZE",
                    "Filesize is not within the allowed boundaries.")
                return None
            #else
        #if
        else :
            self._output.addMessage(__file__, 4, "ERECPARSE",
                                     "This is not a GenBank record.")
            return None
        #else
    #downloadrecord

    def uploadrecord(self, raw_data) :
        """
        Write an uploaded record to a file.
        If the downloaded file is recognised by its hash, the old UD number
        is used.

        @arg raw_data: A GenBank record
        @type raw_data: byte string

        @return: Accession number for the uploaded file.
        @rtype: unicode
        """
        md5sum = self._calcHash(raw_data)

        try:
            reference = Reference.query.filter_by(checksum=md5sum).one()
        except NoResultFound:
            UD = self._newUD()
            if self.write(raw_data, UD, 0):
                reference = Reference(UD, md5sum)
                session.add(reference)
                session.commit()
                return UD
        else:
            if os.path.isfile(self._nametofile(reference.accession)):
                return reference.accession
            else:
                return self.write(raw_data, reference.accession, 0) and reference.accession
    #uploadrecord

    def loadrecord(self, identifier):
        """
        Load a RefSeq record and return it.

        The record is found by trying the following options in order:

        1. Returned from the cache if it is there.
        2. Re-created (if it was created by slicing) or re-downloaded (if it
           was created by URL) if we have information on its source in the
           database.
        3. Fetched from the NCBI.

        :arg identifier: A RefSeq accession number or geninfo identifier (GI).
        :type identifier: unicode

        :return: A parsed RefSeq record or `None` if no record could be found
            for the given identifier.
        :rtype: mutalyzer.GenRecord.Record
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
            filename = self._nametofile(reference.accession)

            if os.path.isfile(filename):
                # It is still in the cache, so filename is valid.
                pass

            elif reference.slice_accession:
                # It was previously created by slicing.
                cast_orientation = {None: None,
                                    'forward': 1,
                                    'reverse': 2}
                if not self.retrieveslice(reference.slice_accession,
                                          reference.slice_start,
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
                self._output.addMessage(__file__, 4, 'ERETR',
                                        'Please upload this sequence again.')
                filename = None

        # If filename is None, we could not retrieve the record.
        if filename is None:
            # Notify batch job to skip all instance of identifier.
            self._output.addOutput('BatchFlags', ('S1', identifier))
            return None

        # Now we have the file, so we can parse it.
        GenBankParser = genbank.GBparser()
        record = GenBankParser.create_record(filename)

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
    #loadrecord
#GenBankRetriever

class LRGRetriever(Retriever):
    """
    Retrieve a LRG record from either the cache or the web.

    Public methods:
        - loadrecord(identifier) ; Load a record, store it in the cache, manage
                                   the cache and return the record.
    """

    def __init__(self, output):
        #TODO documentation
        """
        Initialize the class.

        @todo: documentation
        @arg  output:
        @type  output:
        @arg  database:
        @type  database:
        """
        # Recall init of parent
        Retriever.__init__(self, output)
        self.fileType = "xml"
        # Child specific init
    #__init__

    def loadrecord(self, identifier):
        """
        Load and parse a LRG file based on the identifier

        @arg identifier: The name of the LRG file to read
        @type identifier: unicode

        @return: record ; GenRecord.Record of LRG file
                   None ; in case of failure
        @rtype:
        """

        # Make a filename based upon the identifier.
        filename = self._nametofile(identifier)

        if not os.path.isfile(filename) :   # We can't find the file.
            filename = self.fetch(identifier)

        if filename is None:                # return None in case of error
            #Notify batch to skip all instance of identifier
            self._output.addOutput("BatchFlags", ("S1", identifier))
            return None

        # Now we have the file, so we can parse it.
        file_handle = bz2.BZ2File(filename, "r")
        file_handle = codecs.getreader('utf-8')(file_handle)

        #create GenRecord.Record from LRG file
        record = lrg.create_record(file_handle.read())
        file_handle.close()

        # We don't create LRGs from other sources, so id is always the same
        # as source_id.
        record.id = identifier
        record.source_id = identifier

        return record
    #loadrecord

    def fetch(self, name):
        """
        Fetch the LRG file and store in the cache directory. First try to
        grab the file from the confirmed section, if this fails, get it
        from the pending section.

        @arg name: The name of the LRG file to fetch
        @type name: unicode

        @return: the full path to the file; None in case of an error
        @rtype: unicode
        """

        prefix = settings.LRG_PREFIX_URL
        url        = prefix + "%s.xml"          % name
        pendingurl = prefix + "pending/%s.xml"  % name

        try:
            return self.downloadrecord(url, name)
        except urllib2.URLError: #Catch error: file not found
            pass

        try:                # Try to get the file from the pending section
            filename = self.downloadrecord(pendingurl, name)
            self._output.addMessage(__file__, 2, "WPEND",
                "Warning: LRG file %s is a pending entry." % name)
            return filename
        except urllib2.URLError:
            self._output.addMessage(__file__, 4, "ERETR",
                                 "Could not retrieve %s." % name)
            return None             #Explicit return in case of an Error
    #fetch

    def downloadrecord(self, url, name = None) :
        """
        Download an LRG record from an URL.

        @arg url: Location of the LRG record
        @type url: unicode

        @return:
            - filename    ; The full path to the file
            - None        ; in case of failure
        @rtype: unicode
        """

        lrgID = name or os.path.splitext(os.path.split(url)[1])[0]
        #if not lrgID.startswith("LRG"):
        #    return None
        filename = self._nametofile(lrgID)

        # Todo: Properly read the file contents to a unicode string and write
        #   it utf-8 encoded.
        handle = urllib2.urlopen(url)
        info = handle.info()
        if info["Content-Type"] == "application/xml" and info.has_key("Content-length"):

            length = int(info["Content-Length"])
            if 512 < length < settings.MAX_FILE_SIZE:
                raw_data = handle.read()
                handle.close()

                #Do an md5 check
                md5sum = self._calcHash(raw_data)
                try:
                    reference = Reference.query.filter_by(accession=lrgID).one()
                    md5db = reference.checksum
                except NoResultFound:
                    md5db = None

                if md5db is None:
                    reference = Reference(lrgID, md5sum, download_url=url)
                    session.add(reference)
                    session.commit()
                elif md5sum != md5db:       #hash has changed for the LRG ID
                    self._output.addMessage(__file__, -1, "WHASH",
                        "Warning: Hash of %s changed from %s to %s." % (
                        lrgID, md5db, md5sum))
                    Reference.query.filter_by(accession=lrgID).update({'checksum': md5sum})
                    session.commit()
                else:                       #hash the same as in db
                    pass

                if not os.path.isfile(filename) :
                    return self.write(raw_data, lrgID)
                else:
                    # This can only occur if synchronus calls to mutalyzer are
                    # made to recover a file that did not exist. Still leaves
                    # a window in between the check and the write.
                    return filename
            #if
            else :
                self._output.addMessage(__file__, 4, "EFILESIZE",
                    "Filesize is not within the allowed boundaries.")
        #if
        else :
            self._output.addMessage(__file__, 4, "ERECPARSE",
                                     "This is not an LRG record.")
        handle.close()
    #downloadrecord

    def write(self, raw_data, filename) :
        """
        Write raw LRG data to a file. The data is parsed before writing,
        if a parse error occurs None is returned.

        @arg raw_data: The data
        @type raw_data: string
        @arg filename: The intended name of the file
        @type filename: unicode

        @return:
            - filename ; The full path and name of the file written
            - None     ; In case of an error
        @rtype: unicode
        """
        # Dirty way to test if a file is valid,
        # Parse the file to see if it's a real LRG file.
        try:
            lrg.create_record(raw_data)
        except DOMException:
            self._output.addMessage(__file__, 4, "ERECPARSE",
                                      "Could not parse file.")
            return None             # Explicit return on Error

        return self._write(raw_data, filename) #returns full path
    #write
#LargeRetriever

if __name__ == "__main__" :
    pass
#if
