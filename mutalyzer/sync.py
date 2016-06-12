"""
Synchronizing the reference file cache with other Mutalyzer instances.
"""


from __future__ import unicode_literals

from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

from datetime import datetime, timedelta
import os
import re
import urllib2

from sqlalchemy.orm.exc import NoResultFound
from suds.client import Client

from mutalyzer.config import settings
from mutalyzer.db import session
from mutalyzer.db.models import Reference
from mutalyzer import Retriever


#: Default number of days to go back in synchronization.
DEFAULT_CREATED_SINCE_DAYS = 7


class CacheSync(object):
    """
    Synchronize the database cache with other Mutalyzer instances.
    """
    def __init__(self, output):
        """
        Instantiate the object.

        :arg output: An output object.
        :type output: :class:`mutalyzer.output.Output`
        """
        self._output = output

    def local_cache(self, created_since=None):
        """
        Get all entries in the local cache with creation date `created_since`
        or later.

        :arg created_since: Only entries with this creation date or later
          are returned.
        :type created_since: datatime.datetime

        :return: List of cache entries.
        :rtype: list(dict)
        """
        if not created_since:
            created_since = datetime.today() - \
                            timedelta(days=DEFAULT_CREATED_SINCE_DAYS)

        references = Reference.query.filter(Reference.added >= created_since)
        cache = []

        # Translate each entry to a dictionary and check if it is cached on
        # our filesystem.
        for reference in references:
            # Note that this way we only include Genbank files, not LRG files.
            # But we are only really interested in manually uploaded files
            # anyway, which can currently only be Genbank files.
            cached = None
            if os.path.isfile(os.path.join(settings.CACHE_DIR,
                                           '%s.gb.bz2' % reference.accession)):
                cached = '%s.gb' % reference.accession
            cache.append({'name':                  reference.accession,
                          'source':                reference.source,
                          'source_data':           reference.source_data,
                          'hash':                  reference.checksum,
                          'created':               reference.added,
                          'cached':                cached})

        return cache

    def remote_cache(self, remote_wsdl, created_since=None):
        """
        Get all entries in the remote cache with creation date `created_since`
        or later.

        :arg remote_wsdl: The url of the remote SOAP WSDL description.
        :type remote_wsdl: unicode
        :arg created_since: Only entries with this creation date or later
          are returned.
        :type created_since: datatime.datetime

        :return: List of cache entries.
        :rtype: list(dict)
        """
        self._output.addMessage(__file__, -1, 'INFO', 'Getting remote cache'
                                ' from %s' % remote_wsdl)

        if not created_since:
            created_since = datetime.today() - \
                            timedelta(days=DEFAULT_CREATED_SINCE_DAYS)
        client = Client(remote_wsdl, cache=None)
        cache = client.service.getCache(created_since)

        def cache_entry_from_soap(entry):
            """
            Create a nice dictionary out of the CacheEntry object.
            """
            entry_dict = {'name':    entry.name,
                          'hash':    entry.hash,
                          'created': entry.created}
            for attribute in ('source', 'source_data', 'cached'):
                entry_dict[attribute] = entry[attribute] \
                                        if attribute in entry else None
            return entry_dict

        return map(cache_entry_from_soap, cache.CacheEntry)

    def store_remote_file(self, name, url):
        """
        Download a remote file located at `url` and store it as `name`.

        :arg name: Name to store the file under.
        :type name: unicode
        :arg url: Url to the remote file.
        :type url: unicode
        """
        if not re.match('^[\da-zA-Z\._-]+$', name):
            return

        # Download remote data
        handle = urllib2.urlopen(url)
        data = handle.read()
        handle.close()

        # Store remote data
        retriever = Retriever.GenBankRetriever(self._output)
        retriever.write(data, name, 0)

    def sync_with_remote(self, remote_wsdl, url_template,
                         days=DEFAULT_CREATED_SINCE_DAYS):
        """
        Synchronize the local cache with the remote cache.

        ::

            >>> wsdl = 'https://mutalyzer.nl/mutalyzer/services/?wsdl'
            >>> template = 'https://mutalyzer.nl/mutalyzer/Reference/{file}'
            >>> self.sync_with_remote(wsdl, template)
            (14, 3)

        :arg remote_wsdl: The url of the remote SOAP WSDL description.
        :type remote_wsdl: unicode
        :arg url_template: Formatting string containing a ``{file}``
          occurence, see example usage above.
        :string url_template: unicode
        :arg days: Only remote entries added this number of days ago or
          later are considered.
        :type days: int

        :return: The number of entries added to the local cache and the number
          cache files downloaded from the remote site.
        :rtype: tuple(int, int)
        """
        self._output.addMessage(__file__, -1, 'INFO', 'Starting cache sync')

        created_since = datetime.today() - timedelta(days=days)
        remote_cache = self.remote_cache(remote_wsdl, created_since)

        inserted = downloaded = 0

        for entry in remote_cache:
            try:
                reference = Reference.query.filter_by(accession=entry['name']).one()
                if reference.checksum is not None:
                    continue
            except NoResultFound:
                pass

            if Reference.query.filter_by(checksum=entry['hash']).count() > 0:
                continue

            reference = Reference(entry['name'], entry['hash'], entry['source'],
                                  source_data=entry['source_data'])
            session.add(reference)
            session.commit()
            inserted += 1
            if entry['source'] == 'upload' and entry['cached']:
                url = url_template.format(file=entry['cached'])
                self.store_remote_file(entry['name'], url)
                downloaded += 1

        self._output.addMessage(__file__, -1, 'INFO',
                                'Inserted %d entries in the cache,'
                                ' downloaded %d files.' \
                                % (inserted, downloaded))
        self._output.addMessage(__file__, -1, 'INFO', 'Finished cache sync')

        return inserted, downloaded
