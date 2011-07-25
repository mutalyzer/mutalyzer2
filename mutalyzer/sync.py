"""
Module for synchronizing the database with other Mutalyzer instances.
"""


from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()
from suds.client import Client
from datetime import datetime, timedelta


DEFAULT_CREATED_SINCE_DAYS = 7


class CacheSync(object):
    """
    Todo.
    """
    def __init__(self, config, database):
        """
        Todo.
        """
        self._config = config
        self._database = database

    def local_cache(self, created_since=None):
        """
        Todo.
        """
        if not created_since:
            created_since = datetime.today() - \
                            timedelta(days=DEFAULT_CREATED_SINCE_DAYS)
        return self._database.getGB(created_since)

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
            entry_dict =  {'name':    entry.name,
                           'hash':    entry.hash,
                           'created': entry.created}
            for attribute in ('gi', 'chromosomeName', 'chromosomeStart'
                              'chromosomeStop', 'chromosomeOrientation',
                              'url'):
                entry_dict[attribute] = entry[attribute] \
                                        if attribute in entry else None
            return entry_dict

        return map(cache_entry_from_soap, cache.CacheEntry)

    def sync_with_remote(self, remote_wsdl, created_since=None):
        """
        Todo.
        """
        remote_cache = self.remote_cache(remote_wsdl, created_since)

        for entry in remote_cache:
            if self._database.getHash(entry['name']):
                continue
            #self._database.insertGB(entry['name'],
            #                        entry['gi'],
            #                        entry['hash'],
            #                        entry['chromosomeName'],
            #                        entry['chromosomeStart'],
            #                        entry['chromosomeStop'],
            #                        entry['chromosomeOrientation'],
            #                        entry['url'])
            print 'inserting %s' % entry['name']
            print entry
            if not entry['chromosomeName'] and not entry['url']:
                print '(downloading file from remote cache...)'
