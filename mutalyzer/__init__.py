"""
HGVS variant nomenclature checker.
"""


import os


# On the event of a new release, we update the __version_info__ and __date__
# package globals and set RELEASE to True.
# Before a release, a development version is denoted by a __version_info__
# ending with a 'dev' item. Also, RELEASE is set to False (indicating that
# the __date__ value is to be ignored).
#
# We follow a versioning scheme compatible with setuptools [1] where the
# __version_info__ variable always contains the version of the upcomming
# release (and not that of the previous release), post-fixed with a 'dev'
# item. Only in a release commit, this 'dev' item is removed (and added
# again in the next commit).
#
# [1] http://peak.telecommunity.com/DevCenter/setuptools#specifying-your-project-s-version

RELEASE = False

__version_info__ = ('2', '0', 'beta-17', 'dev')
__date__ = '1 Mar 2012'


__version__ = '.'.join(__version_info__)
__author__ = 'Leiden University Medical Center'
__contact__ = 'humgen@lumc.nl'
__homepage__ = 'http://mutalyzer.nl'


NOMENCLATURE_VERSION_INFO = ('2', '0')
NOMENCLATURE_VERSION = '.'.join(NOMENCLATURE_VERSION_INFO)

COPYRIGHT_YEARS = (2007, int(__date__[-4:]))

SOAP_NAMESPACE = 'http://mutalyzer.nl/2.0/services'


def package_root():
    """
    Get the absolute path to the mutalyzer package. This is usefull for
    things like locating HTML templates (which are in a subdirectory of the
    package).

    @return: Absolute path to the mutalyzer package.
    @rtype:  string
    """
    return os.path.realpath(os.path.split(__file__)[0])


def is_test():
    """
    Check if we are in a test environment. This is determined by the
    MUTALYZER_ENV environment variable, which should then be set to 'test'.

    @return: True if we are in a test environment, False otherwise.
    @rtype:  bool
    """
    return 'MUTALYZER_ENV' in os.environ \
           and os.environ['MUTALYZER_ENV'] == 'test'
