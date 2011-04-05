import sys
from setuptools import setup, find_packages

if sys.version_info < (2, 6):
    raise Exception('Mutalyzer requires Python 2.6 or higher.')

import mutalyzer as distmeta

setup(
    name='mutalyzer',
    version=distmeta.__version__,
    description=distmeta.__doc__,
    author=distmeta.__author__,
    author_email=distmeta.__contact__,
    url=distmeta.__homepage__,
    packages=find_packages(exclude=['doc', 'extras', 'tests']),
    include_package_data=True,
    scripts=['bin/batch_daemon'],
    zip_safe=False
    #data_files=[('/etc/init.d', ['examples/init.d/mutalyzer-batchd'])]
    # Todo: templates, doc, bin
)
