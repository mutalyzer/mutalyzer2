"""
Module for synchronizing the database with other Mutalyzer instances.
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
    Todo.
    """
    def __init__(self, config, output, database):
        """
        Todo.
        """
        self._config = config
        self._output = output
        self._database = database

    def local_cache(self, created_since=None):
        """
        Todo.
        """
        if not created_since:
            created_since = datetime.today() - \
                            timedelta(days=DEFAULT_CREATED_SINCE_DAYS)
        cache = self._database.getGBSince(created_since)

        entries = []

        # For each entry, check if it is cached on our filesystem.
        # Todo: refactor
        for entry in cache:
            e = list(entry)
            # Note that this way we only include Genbank files, not LRG files.
            file_name = '%s.gb.bz2' % entry[0]
            file_path = os.path.join(self._config.Retriever.cache, file_name)
            if os.path.isfile(file_path):
                e.append('%s.gb' % entry[0])
            else:
                e.append(None)
            entries.append(e)

        return entries

    def remote_cache(self, remote_wsdl, created_since=None):
        """
        Todo.
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
        Todo.
        """
        remote_cache = self.remote_cache(remote_wsdl, created_since)

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
            print 'inserting %s' % entry['name']
            if not entry['chromosomeName'] and not entry['url']:
                if entry['cached']:
                    print 'downloading file from remote cache: %s' % (url_template % str(entry['cached']))
                    self.store_remote_file(entry['name'], url_template % entry['cached'])
                else:
                    print 'cannot download this file from remote cache'

    def store_remote_file(self, name, url):
        """
        Todo.
        """
        if not re.match('^[\da-zA-Z\._-]+$', name):
            return

        # Download remote data
        handle = urllib2.urlopen(url)
        data = handle.read()
        handle.close()

        # Store remote data
        retriever = Retriever.GenBankRetriever(self._config.Retriever, self._output, self._database)
        retriever.write(data, name, 0)
