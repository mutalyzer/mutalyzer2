"""
Module for synchronizing the database with other Mutalyzer instances.
"""


from datetime import datetime, timedelta


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
            created_since = datetime.today() - timedelta(days=7)
        return self._database.getGB(created_since)
