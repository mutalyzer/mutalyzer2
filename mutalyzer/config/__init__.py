"""
Module for reading the configuration values from a configuration file.

Default values will be read from the :mod:`mutalyzer.config.default_settings`
module and overridden by any values from the module specified by the
`MUTALYZER_SETTINGS`.

Alternatively, the default values can be overridden manually using the
:meth:`settings.configure` method. If this is done before the first use of a
configuration value, the `MUTALYZER_SETTINGS` environment variable will never
be used.
"""


from __future__ import unicode_literals

import collections
import os
import warnings


import flask.config

from mutalyzer import util


#: Environment variable used for locating the settings module.
ENVIRONMENT_VARIABLE = 'MUTALYZER_SETTINGS'


class Settings(flask.config.Config, util.AttributeDictMixin):
    """
    Dictionary with some extra ways to fill it from files or special
    dictionaries (see `flask.config.Config`) and attribute access.
    """
    def __init__(self):
        # We fix the root_path argument to the current working directory.
        super(Settings, self).__init__(os.getcwd())


class LazySettings(util.LazyObject):
    """
    A lazy proxy for a settings object.

    Taken from `Django <https://www.djangoproject.com/>`_
    (`django.conf.LazySettings`).

    Configuration settings can be updated with the :meth:`configure` method.

    The user can register callbacks to configuration keys that are called
    whenever the value for that key is updated with :meth:`on_update`.

    .. note:: Django also does some logging config magic here, we did not copy
        that.
    """
    def __init__(self, *args, **kwargs):
        # Assign to __dict__ to avoid __setattr__ call.
        self.__dict__['_callbacks'] = collections.defaultdict(list)
        super(LazySettings, self).__init__(*args, **kwargs)

    def _setup(self, from_environment=True):
        """
        Load the settings module pointed to by the environment variable. This
        is used the first time we need any settings at all, if the user has not
        previously configured the settings manually with :meth:`configure`.
        """
        self._wrapped = Settings()
        self._wrapped.from_object('mutalyzer.config.default_settings')
        if from_environment:
            if ENVIRONMENT_VARIABLE in os.environ:
                self._wrapped.from_envvar(ENVIRONMENT_VARIABLE)
            else:
                warnings.warn('The environment variable \'%s\' is not set '
                              'and as such default configuration settings '
                              'are used. Set this variable and make it point '
                              'to a configuration file to customize '
                              'configuration.' % ENVIRONMENT_VARIABLE,
                              RuntimeWarning)

    def configure(self, settings):
        """
        Called to manually configure the settings.
        """
        if self._wrapped is util.empty:
            self._setup(from_environment=False)
        self._wrapped.update(settings)

        # Callbacks for specific keys.
        for key, callbacks in self._callbacks.items():
            if key in settings:
                for callback in callbacks:
                    callback(settings[key])

        # General callbacks.
        for callback in self._callbacks[None]:
            callback(settings)

    @property
    def configured(self):
        """
        Returns True if the settings have already been configured.
        """
        return self._wrapped is not util.empty

    def on_update(self, callback, key=None):
        """
        Register a callback for the update of a key (or any update if `key` is
        `None`).

        The callback is called with as argument the new value for the updated
        key (or a dictionary with all updated key-value pairs if `key` is
        `None`).
        """
        self._callbacks[key].append(callback)


#: Global :class:`LazySettings` instance. Use this for querying configuration
#: settings.
settings = LazySettings()
