"""
Module for synchronizing the database with other Mutalyzer instances.
"""


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

    def local_cache(self):
        """
        Todo.
        """
        return self._database.getGB()
