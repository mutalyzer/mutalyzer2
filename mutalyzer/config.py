"""
Module for reading the config file and splitting up the variables into
subclasses. Each of these subclasses are used to configure a specific
module.
"""


import os
from configobj import ConfigObj


SYSTEM_CONFIGURATION = '/etc/mutalyzer/config'
USER_CONFIGURATION = os.path.join(
    os.environ.get('XDG_CONFIG_HOME', None) or \
    os.path.join(os.path.expanduser('~'), '.config'),
    'mutalyzer', 'config')


class ConfigurationError(Exception):
    pass


class Config():
    """
    Read the configuration file and store the data in subclasses.
    """
    class Retriever(): pass
    class Db(): pass
    class Output(): pass
    class Mutator(): pass
    class Scheduler(): pass
    class Batch(): pass
    class File(): pass
    class GenRecord(): pass

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

        By the DRY-principle, we don't enumerate the configuration variables
        for each class in documentation. Instead, what variables are used by
        each class is easy to see from the code below.

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

            # Set the variables needed by the Retriever module.
            self.Retriever.email = config["email"]
            self.Retriever.cache = config["cache"]
            self.Retriever.cachesize = int(config["cachesize"]) * 1048576
            self.Retriever.maxDldSize = int(config["maxDldSize"]) * 1048576
            self.Retriever.minDldSize = int(config["minDldSize"])
            self.Retriever.lrgURL = config["lrgurl"]

            # Set the variables needed by the Db module.
            self.Db.internalDb = config["internalDb"]
            self.Db.dbNames = config["dbNames"]
            self.Db.LocalMySQLuser = config["LocalMySQLuser"]
            self.Db.LocalMySQLhost = config["LocalMySQLhost"]

            # Set the variables needed by the Output module.
            self.Output.log = config["log"]
            self.Output.datestring = config["datestring"]
            self.Output.loglevel = int(config["loglevel"])
            self.Output.outputlevel = int(config["outputlevel"])
            self.Output.debug = config.as_bool('debug')

            # Set the variables needed by the Mutator module.
            self.Mutator.flanksize = int(config["flanksize"])
            self.Mutator.maxvissize = int(config["maxvissize"])
            self.Mutator.flankclipsize = int(config["flankclipsize"])

            # Set the variables needed by the Scheduler module.
            self.Scheduler.mailFrom = config["mailFrom"]
            self.Scheduler.mailSubject = config["mailSubject"]
            self.Scheduler.resultsDir = config["resultsDir"]
            self.Scheduler.nameCheckOutHeader = config["nameCheckOutHeader"]
            self.Scheduler.syntaxCheckOutHeader = config["syntaxCheckOutHeader"]
            self.Scheduler.positionConverterOutHeader = config["positionConverterOutHeader"]
            self.Scheduler.snpConverterOutHeader = config["snpConverterOutHeader"]

            # Set thte variables neede for the Batch module.
            self.Batch.PIDfile = config["PIDfile"]
            self.Batch.batchInputMaxSize = int(config["batchInputMaxSize"]) * 1048576

            # Set the variables needed by the File module.
            self.File.bufSize = int(config["bufSize"])
            self.File.header = config["header"]
            self.File.threshold = float(config["threshold"])

            # Set the variables needed by the File module.
            self.GenRecord.spliceAlarm = int(config["spliceAlarm"])
            self.GenRecord.spliceWarn = int(config["spliceWarn"])

        except KeyError as e:
            raise ConfigurationError('Missing configuration value: %s' % e)
    #__init__

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
#Config
