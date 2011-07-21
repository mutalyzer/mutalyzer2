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
    license='Not distributable',
    platforms=['any'],
    packages=find_packages(exclude=['doc', 'extras', 'tests']),
    include_package_data=True,
    scripts=['bin/mutalyzer',
             'bin/mutalyzer-cache-sync',
             'bin/mutalyzer-batchd',
             'bin/mutalyzer-ucsc-update',
             'bin/mutalyzer-website.wsgi',
             'bin/mutalyzer-webservice.wsgi'],
    zip_safe=False
)

# Things not handled by this setup.py:
# - Copy extras/config.example to /etc/mutalyzer/config
# - Database setup
# - Chown /var/log/mutalyzer.log and /var/cache/mutalyzer
# - Copy extras/init.d/mutalyzer-batchd to /etc/init.d/mutalyzer-batchd
# - Copy doc to /usr/share/doc
# Check extras/post-install.sh for these.
