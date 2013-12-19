"""
Module for reading the configuration values from a configuration file.

Default values will be read from the :mod:`mutalyzer.config.default_settings`
module and overridden by any values from the module specified by the
`MUTALYZER_SETTINGS`.

Alternatively, the default values can be overridden manually using the
:meth:`settings.configure` method before the first use of a configuration
value, in which case the `MUTALYZER_SETTINGS` environment variable will not be
used.
"""


import flask.config
import os

from mutalyzer import util


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

    .. note:: Django also does some logging config magic here, we did not copy
        that.
    """
    def _setup(self, from_environment=True):
        """
        Load the settings module pointed to by the environment variable. This
        is used the first time we need any settings at all, if the user has not
        previously configured the settings manually with :meth:`configure`.
        """
        self._wrapped = Settings()
        self._wrapped.from_object('mutalyzer.config.default_settings')
        if from_environment:
            self._wrapped.from_envvar(ENVIRONMENT_VARIABLE)

    def configure(self, settings):
        """
        Called to manually configure the settings.
        """
        if self._wrapped is None:
            self._setup(from_environment=False)
        self._wrapped.update(settings)

    @property
    def configured(self):
        """
        Returns True if the settings have already been configured.
        """
        return self._wrapped is not None


settings = LazySettings()
