"""HGVS variant nomenclature checker."""


# On the event of a new release, we update the __version_info__ and __date__
# package globals and set RELEASE to True.
# After a release, a development version is denoted by a __version_info__
# ending with a 'dev' item. Also, RELEASE is set to False (indicating that
# the __date__ value is to be ignored).

RELEASE = False

__version_info__ = ('2', '0', 'beta-9', 'dev')
__date__ = '31 Jan 2011'


__version__ = '.'.join(__version_info__)
__author__ = 'Leiden University Medical Center'
__contact__ = 'humgen@lumc.nl'
__homepage__ = 'http://mutalyzer.nl'


NOMENCLATURE_VERSION_INFO = ('2', '0')
NOMENCLATURE_VERSION = '.'.join(NOMENCLATURE_VERSION_INFO)
