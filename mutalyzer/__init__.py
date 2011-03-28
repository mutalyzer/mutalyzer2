"""
Mutalyzer, a HGVS variant nomenclature checker.

Developed at the Human Genetics department of the Leiden University
Medical Center.
"""


# On the event of a new release, we update the __version_info__ and __date__
# package globals and set RELEASE to True.
# After a release, a development version is denoted by a __version_info__
# ending with a 'dev' item. Also, RELEASE is set to False (indicating that
# the __date__ value is to be ignored).

__author__ = 'Leiden University Medical Center'
__version_info__ = ('2', '0', 'beta-9', 'dev')
__version__ = '.'.join(__version_info__)
__date__ = '31 Jan 2011'

RELEASE = False

NOMENCLATURE_VERSION_INFO = ('2', '0')
NOMENCLATURE_VERSION = '.'.join(NOMENCLATURE_VERSION_INFO)
