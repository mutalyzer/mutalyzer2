"""
Module for reading the configuration values from configuration files.

All communication with this module should be done by using the get function
which returns a configuration value, given a name.

Reading the configuration file is implemented lazily and as such done upon the
first call to the get function.

Configuration values are read from two locations, in this order:
1) /etc/mutalyzer/config
2) $XDG_CONFIG_HOME/mutalyzer/config

If both files exist, values defined in the second overwrite values defined in
the first.
"""


import os
from configobj import ConfigObj

from mutalyzer.util import singleton


SYSTEM_CONFIGURATION = '/etc/mutalyzer/config'
USER_CONFIGURATION = os.path.join(
    os.environ.get('XDG_CONFIG_HOME', None) or \
    os.path.join(os.path.expanduser('~'), '.config'),
    'mutalyzer', 'config')


class ConfigurationError(Exception):
    """
    Raised when a configuration file cannot be read.
    """
    pass


def get(name):
    """
    Get a configuration value by name.

    @arg name: Name for the configuration value.
    @type name: string

    @raise ConfigurationError: If configuration value could not be read.
        Reasons are:
        - Configuration file could not be parsed.
        - Not all variables are present in configuration file.
        - Given configuration value name does not exist.
    """
    return _Config().get(name)


@singleton
class _Config():
    """
    Read the configuration file and provide access to its values.

    Please note the limitations from the use of the @singleton decorator as
    described in its docstring.
    """
    def __init__(self, filename=None):
        """
        Initialise the class with variables read from the configuration
        file. In principle, this is the only place in the code where a
        hard coded constant is used (the name and path to the configuration
        file).

        Configuration values are read from two locations, in this order:
        1) /etc/mutalyzer/config
        2) $XDG_CONFIG_HOME/mutalyzer/config

        If both files exist, values defined in the second overwrite values
        defined in the first.

        An exception to this system is when the optional {filename} argument
        is set. In that case, the locations listed above are ignored and the
        configuration is read from {filename}.

        @kwarg filename: Optional filename to read configuration from. If
            present, this overrides automatic detection of configuration file
            location.
        @type filename: string

        @raise ConfigurationError: If configuration could not be read.
            Reasons are:
            - Supplied argument {filename} could not be opened.
            - Configuration file could not be parsed.
            - Not all variables are present in configuration file.
        """
        config = None

        if filename:
            config = self._load_config(filename)
        else:
            if os.path.isfile(SYSTEM_CONFIGURATION):
                config = self._load_config(SYSTEM_CONFIGURATION)
            if os.path.isfile(USER_CONFIGURATION):
                user_config = self._load_config(USER_CONFIGURATION)
                if config:
                    config.merge(user_config)
                else:
                    config = user_config

        if not config:
            raise ConfigurationError('Could not locate configuration.')

        try:

            # We explicitely read all configuration values ad store them in
            # our own dictionary. This makes sure we notice missing or
            # incorrect values upon instantiation.

            # A few 'special' values.
            self._values = {'autoReconnect': config.as_bool('autoReconnect'),
                            'debug':         config.as_bool('debug'),
                            'piwik':         config.as_bool('piwik'),
                            'threshold':     float(config['threshold'])}

            # Simple string values.
            for name in ('email', 'cache', 'lrgurl', 'internalDb', 'dbNames',
                         'LocalMySQLuser', 'LocalMySQLhost', 'log',
                         'datestring', 'mailFrom', 'mailSubject',
                         'resultsDir', 'nameCheckOutHeader', 'defaultDb',
                         'syntaxCheckOutHeader', 'positionConverterOutHeader',
                         'snpConverterOutHeader', 'PIDfile', 'header',
                         'piwikBase'):
                self._values[name] = config[name]

            # Simple integer values.
            for name in ('minDldSize', 'loglevel', 'outputlevel', 'flanksize',
                         'maxvissize', 'flankclipsize', 'bufSize',
                         'spliceAlarm', 'spliceWarn', 'piwikSite',
                         'proteinLinkLifetime', 'proteinLinkNoneLifetime'):
                self._values[name] = int(config[name])

            # File sizes (given in megabytes, stored in bytes).
            for name in ('cachesize', 'maxDldSize', 'batchInputMaxSize'):
                self._values[name] = int(config[name]) * 1048576

        except KeyError as e:
            raise ConfigurationError('Missing configuration value: %s' % e)
    #__init__

    def get(self, name):
        """
        Get a configuration value by name.

        @arg name: Name for the configuration value.
        @type name: string

        @raise ConfigurationError: If given configuration value name does not
            exist.
        """
        try:
            return self._values[name]
        except KeyError:
            raise ConfigurationError('No such configuration value: %s' % name)
    #get

    def _load_config(self, filename):
        """
        Create a ConfigObj from the configuration in {filename}.
        """
        try:
            return ConfigObj(filename)
        except IOError:
            raise ConfigurationError('Could not open configuration file: %s' \
                                     % filename)
        except SyntaxError:
            raise ConfigurationError('Could not parse configuration file: %s' \
                                     % filename)
    #_load_config
#_Config
