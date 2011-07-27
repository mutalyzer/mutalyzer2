"""
Module for reading the config file and splitting up the variables into
subclasses. Each of these subclasses are used to configure a specific
module.
"""


import os
import tempfile
from configobj import ConfigObj

import mutalyzer
from mutalyzer import util


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

        The configuration file location is automatically detected, in the
        following order:
        1) $XDG_CONFIG_HOME/mutalyzer/config
        2) /etc/mutalyzer/config

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

        @todo: Be able to have /etc/mutalyzer/config as 'base' configuration
            and some additional configuration in ~/.mutalyzer/config to
            overwrite the base configuration.
            For example, a normal user could use a different cache directory
            (writable by the user) than the system wide Mutalyzer config.
        """
        if filename is None:
            base = os.environ.get('XDG_CONFIG_HOME', None)
            if base is None:
                base = os.path.join(os.path.expanduser('~'), '.config')
            filename = os.path.join(base, 'mutalyzer', 'config')

            if not os.path.isfile(filename):
                filename = '/etc/mutalyzer/config'
                if not os.path.isfile(filename):
                    raise ConfigurationError('Could not locate configuration.')

        try:
            config = ConfigObj(filename)
        except IOError:
            raise ConfigurationError('Could not open configuration file: %s' \
                                     % filename)
        except SyntaxError:
            raise ConfigurationError('Could not parse configuration file.')

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
            self.Db.RemoteMySQLuser = config["RemoteMySQLuser"]
            self.Db.RemoteMySQLhost = config["RemoteMySQLhost"]
            self.Db.UpdateInterval = int(config["UpdateInterval"])

            # Set the variables needed by the Output module.
            self.Output.log = config["log"]
            self.Output.datestring = config["datestring"]
            self.Output.loglevel = int(config["loglevel"])
            self.Output.outputlevel = int(config["outputlevel"])

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


            # If we are in a testing environment, use a temporary file for
            # logging and a temporary directory for the cache.
            # We don't remove these after the tests, since they might be
            # useful for debugging.
            if mutalyzer.is_test():
                # Todo:
                #
                # This needs some refactoring. The problem with the temporary
                # file and dir names is that they will not be used by the
                # (running) batch daemon, which will thus save its results to
                # to 'normal' directory.
                # Furthermore, subsequent web requests from a unit test will
                # use different configuration instantiations, so might not
                # see results from previous requests.
                #
                # We need a more robust solution for different configurations,
                # depending of the running user/setting (e.g. unit tests).
                #
                # Idea: Don't create a local instance of the website in the
                # unit tests, but only use running instances of all servers
                # (website, webservice, batch daemon). They will use their
                # own 'normal' configuration.
                # All other parts of the unit tests will use temporary test
                # configuration values. We might even decorate the tests
                # needing server access as such and provide the option of
                # skipping these.

                #handle, filename = tempfile.mkstemp(suffix='.log',
                #                                    prefix='mutalyzer-tests-')
                #os.close(handle)
                #self.Output.log = filename
                #dirname = tempfile.mkdtemp(suffix='.cache',
                #                           prefix='mutalyzer-tests-')
                #self.Retriever.cache = dirname
                #self.Scheduler.resultsDir = dirname
                pass

        except KeyError as e:
            raise ConfigurationError('Missing configuration value: %s' % e)
    #__init__
#Config
