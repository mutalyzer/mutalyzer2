"""
Module for synchronizing the database with other Mutalyzer instances.

Todo: add some logging to the output object.
"""


from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import os
import re
from datetime import datetime, timedelta
import urllib2
from suds.client import Client

from mutalyzer import Retriever


DEFAULT_CREATED_SINCE_DAYS = 7


class CacheSync(object):
    """
    Synchronize the database cache with other Mutalyzer instances.
    """
    def __init__(self, config, output, database):
        """
        Instantiate the object.

        @arg config: A configuration object.
        @type config: mutalyzer.config.Config
        @arg output: An output object.
        @type output: mutalyzer.output.Output
        @arg database: A database object.
        @type database: mutalyzer.Db.Cache
        """
        self._config = config
        self._output = output
        self._database = database

    def local_cache(self, created_since=None):
        """
        Get all entries in the local cache with creation date {created_since}
        or later.

        @kwarg created_since: Only entries with this creation date or later
            are returned.
        @type created_since: datatime.datetime

        @return: List of cache entries.
        @rtype: list(dictionary)
        """
        if not created_since:
            created_since = datetime.today() - \
                            timedelta(days=DEFAULT_CREATED_SINCE_DAYS)

        entries = self._database.getGBSince(created_since)
        cache = []

        # Translate each entry to a dictionary and check if it is cached on
        # our filesystem.
        for entry in entries:
            # Note that this way we only include Genbank files, not LRG files.
            cached = None
            if os.path.isfile(os.path.join(self._config.Retriever.cache,
                                           '%s.gb.bz2' % entry[0])):
                cached = '%s.gb' % entry[0]
            cache.append({'name':                  entry[0],
                          'gi':                    entry[1],
                          'hash':                  entry[2],
                          'chromosomeName':        entry[3],
                          'chromosomeStart':       entry[4],
                          'chromosomeStop':        entry[5],
                          'chromosomeOrientation': entry[6],
                          'url':                   entry[7],
                          'created':               entry[8],
                          'cached':                cached}

        return cache

    def remote_cache(self, remote_wsdl, created_since=None):
        """
        Get all entries in the remote cache with creation date {created_since}
        or later.

        @arg remote_wsdl: The url of the remote SOAP WSDL description.
        @type remote_wsdl: string
        @kwarg created_since: Only entries with this creation date or later
            are returned.
        @type created_since: datatime.datetime

        @return: List of cache entries.
        @rtype: list(dictionary)
        """
        if not created_since:
            created_since = datetime.today() - \
                            timedelta(days=DEFAULT_CREATED_SINCE_DAYS)
        client = Client(remote_wsdl, cache=None)
        cache = client.service.getCache(created_since)

        def cache_entry_from_soap(entry):
            """
            Create a nice dictionary out of the CacheEntry object.
            """
            entry_dict =  {'name':    str(entry.name),
                           'hash':    str(entry.hash),
                           'created': entry.created}
            for attribute in ('gi', 'chromosomeName', 'url', 'cached'):
                entry_dict[attribute] = str(entry[attribute]) \
                                        if attribute in entry else None
            for attribute in ('chromosomeStart', 'chromosomeStop',
                              'chromosomeOrientation'):
                entry_dict[attribute] = int(entry[attribute]) \
                                        if attribute in entry else None
            return entry_dict

        return map(cache_entry_from_soap, cache.CacheEntry)

    def sync_with_remote(self, remote_wsdl, url_template, created_since=None):
        """
        Synchronize the local cache with the remote cache.

            >>> wsdl = 'http://mutalyzer.nl/mutalyzer/services/?wsdl'
            >>> template = 'http://mutalyzer.nl/mutalyzer/Reference/{file}'
            >>> self.sync_with_remote(wsdl, template)
            (14, 3)

        @arg remote_wsdl: The url of the remote SOAP WSDL description.
        @type remote_wsdl: string
        @arg url_template: Formatting string containing a {file} occurence,
            see examle usage above.
        @string url_template: string
        @kwarg created_since: Only remote entries with this creation date or
            later are considered.
        @type created_since: datatime.datetime

        @return: The number of entries added to the local cache and the number
            cache files downloaded from the remote site.
        @rtype: tuple(int, int)
        """
        remote_cache = self.remote_cache(remote_wsdl, created_since)

        inserted = downloaded = 0

        for entry in remote_cache:
            if self._database.getHash(entry['name']):
                continue
            if self._database.getGBFromHash(entry['hash']):
                continue
            if entry['gi'] and self._database.getGBFromGI(entry['gi']):
                continue
            self._database.insertGB(entry['name'],
                                    entry['gi'],
                                    entry['hash'],
                                    entry['chromosomeName'],
                                    entry['chromosomeStart'],
                                    entry['chromosomeStop'],
                                    entry['chromosomeOrientation'],
                                    entry['url'])
            inserted += 1
            if not entry['chromosomeName'] and not entry['url'] \
                   and entry['cached']:
                url = url_template.format(file=entry['cached'])
                self.store_remote_file(entry['name'], url)
                downloaded += 1

        return inserted, downloaded

    def store_remote_file(self, name, url):
        """
        Download a remote file located at {url} and store it as {name}.

        @arg name: Name to store the file under.
        @type name: string
        @arg url: Url to the remote file.
        @type url: string
        """
        if not re.match('^[\da-zA-Z\._-]+$', name):
            return

        # Download remote data
        handle = urllib2.urlopen(url)
        data = handle.read()
        handle.close()

        # Store remote data
        retriever = Retriever.GenBankRetriever(self._config.Retriever,
                                               self._output,
                                               self._database)
        retriever.write(data, name, 0)
