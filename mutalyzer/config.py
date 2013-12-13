"""
Module for reading the configuration values from a configuration file.

All communication with this module should be done by using the get function
which returns a configuration value, given a name.

Reading the configuration file is implemented lazily and as such done upon the
first call to the get function.

Configuration is read from the file specified by the `MUTALYZER_SETTINGS`
environment variable, or `mutalyzer.conf` in the current directory if it is
not set.
"""


import os

from configobj import ConfigObj

from mutalyzer.util import singleton


class ConfigurationError(Exception):
    """
    Raised when a configuration file cannot be read.
    """
    pass


def get(name):
    """
    Get a configuration value by name.

    :arg name: Name for the configuration value.
    :type name: string

    :raise ConfigurationError: If configuration value could not be read.
    """
    return _Config().get(name)


@singleton
class _Config():
    """
    Read the configuration file and provide access to its values.

    Please note the limitations from the use of the @singleton decorator as
    described in its docstring.
    """
    def __init__(self):
        """
        Initialise the class with variables read from the configuration
        file.

        Configuration values are read from the file specified by the
        `MUTALYZER_SETTINGS` environment variable, or `mutalyzer.conf` in the
        current directory if it is not set.

        :raises ConfigurationError: If configuration could not be read.
        """
        filename = os.environ.get('MUTALYZER_SETTINGS', 'mutalyzer.conf')
        config = self._load_config(filename)

        # We define default values for many configuration settings (except for
        # some that are mandatory for the user to define, i.e. those in the
        # extras/config.user.example file).
        # Todo: Do not duplicate default values here and in the example config
        #     file template.
        config.setdefault('cachesize', 50)
        config.setdefault('maxDldSize', 10)
        config.setdefault('minDldSize', 512)
        config.setdefault('lrgurl', 'ftp://ftp.ebi.ac.uk/pub/databases/lrgex/')
        config.setdefault('internalDb', 'mutalyzer')
        config.setdefault('dbNames', ['hg18', 'hg19', 'mm10'])
        config.setdefault('defaultDb', 'hg19')
        config.setdefault('LocalMySQLuser', 'mutalyzer')
        config.setdefault('LocalMySQLhost', 'localhost')
        config.setdefault('autoReconnect', False)
        config.setdefault('proteinLinkLifetime', 30)
        config.setdefault('proteinLinkNoneLifetime', 5)
        config.setdefault('datestring', '%Y-%m-%d %H:%M:%S')
        config.setdefault('loglevel', 3)
        config.setdefault('outputlevel', 1)
        config.setdefault('debug', True)
        config.setdefault('flanksize', 25)
        config.setdefault('maxvissize', 25)
        config.setdefault('flankclipsize', 6)
        config.setdefault('mailFrom', 'noreply@humgen.nl')
        config.setdefault('mailSubject', 'Result of Mutalyzer batch check.')
        config.setdefault('resultsDir', config['cache'])
        config.setdefault('PIDfile', '/var/run/mutalyzer/mutalyzer-batchd.pid')
        config.setdefault('batchInputMaxSize', 5)
        config.setdefault('nameCheckOutHeader',
            ['Input', 'Errors | Messages', 'AccNo', 'Genesymbol', 'Variant',
             'Reference Sequence Start Descr.', 'Coding DNA Descr.',
             'Protein Descr.', 'GeneSymbol Coding DNA Descr.',
             'GeneSymbol Protein Descr.', 'Genomic Reference',
             'Coding Reference', 'Protein Reference', 'Affected Transcripts',
             'Affected Proteins', 'Restriction Sites Created',
             'Restriction Sites Deleted'])
        config.setdefault('syntaxCheckOutHeader', ['Input', 'Status'])
        config.setdefault('positionConverterOutHeader',
            ['Input Variant', 'Errors', 'Chromosomal Variant', 'Coding Variant(s)'])
        config.setdefault('snpConverterOutHeader',
            ['Input Variant', 'HGVS description(s)', 'Errors | Messages'])
        config.setdefault('bufSize', 32768)
        config.setdefault('header', ['AccNo', 'Genesymbol', 'Mutation'])
        config.setdefault('threshold', 0.05)
        config.setdefault('spliceAlarm', 2)
        config.setdefault('spliceWarn', 5)
        config.setdefault('piwik', False)
        config.setdefault('piwikBase', 'https://piwik.example.com')
        config.setdefault('piwikSite', 1)

        try:
            # We explicitely read all configuration values ad store them in
            # our own dictionary. This makes sure we notice missing or
            # incorrect values upon instantiation.

            # A few 'special' values.
            self._values = {'autoReconnect': config.as_bool('autoReconnect'),
                            'debug':         config.as_bool('debug'),
                            'piwik':         config.as_bool('piwik'),
                            'threshold':     config.as_float('threshold')}

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
                self._values[name] = config.as_int(name)

            # File sizes (given in megabytes, stored in bytes).
            for name in ('cachesize', 'maxDldSize', 'batchInputMaxSize'):
                self._values[name] = config.as_int(name) * 1048576

        except KeyError as e:
            raise ConfigurationError('Missing configuration value: %s' % e)

    def get(self, name):
        """
        Get a configuration value by name.

        :arg name: Name for the configuration value.
        :type name: string

        :raises ConfigurationError: If given configuration value name does not
          exist.
        """
        try:
            return self._values[name]
        except KeyError:
            raise ConfigurationError('No such configuration value: %s' % name)

    def _load_config(self, filename):
        """
        Create a `ConfigObj` from the configuration in `filename`.
        """
        try:
            return ConfigObj(filename)
        except IOError:
            raise ConfigurationError('Could not open configuration file: %s'
                                     % filename)
        except SyntaxError:
            raise ConfigurationError('Could not parse configuration file: %s'
                                     % filename)
