"""
HGVS variant nomenclature checker.

Mutalyzer is licensed under the `GNU Affero General Public License
<http://www.gnu.org/licenses/agpl-3.0.html>`_. Please contact the authors if
you want to discuss custom licensing.
"""


from __future__ import unicode_literals


# We follow a versioning scheme compatible with setuptools [1] where the
# package version is always that of the upcoming release (and not that of the
# previous release), post-fixed with ``.dev``. Only in a release commit, the
# ``.dev`` is removed (and added again in the next commit).
#
# Note that this scheme is not 100% compatible with SemVer [2] which would
# require ``-dev`` instead of ``.dev``.
#
# [1] http://peak.telecommunity.com/DevCenter/setuptools#specifying-your-project-s-version
# [2] http://semver.org/

__version_info__ = ('2', '0', '13')
__date__ = '1 Oct 2015'


__version__ = '.'.join(__version_info__)
__author__ = 'Leiden University Medical Center'
__contact__ = 'humgen@lumc.nl'
__homepage__ = 'http://mutalyzer.nl'


NOMENCLATURE_VERSION_INFO = ('2', '0')
NOMENCLATURE_VERSION = '.'.join(NOMENCLATURE_VERSION_INFO)

COPYRIGHT_YEARS = (2009, int(__date__[-4:]))

SOAP_NAMESPACE = 'http://mutalyzer.nl/2.0/services'
