import os
import sys
from setuptools import setup

if sys.version_info < (2, 6):
    raise Exception('Mutalyzer requires Python 2.6 or higher.')

install_requires = []

try:
    with open('README.md') as readme:
        long_description = readme.read()
except IOError:
    long_description = 'See https://mutalyzer.nl'

# This is quite the hack, but we don't want to import our package from here
# since that's recipe for disaster (it might have some uninstalled
# dependencies, or we might import another already installed version).
distmeta = {}
for line in open(os.path.join('mutalyzer', '__init__.py')):
    try:
        field, value = (x.strip() for x in line.split('='))
    except ValueError:
        continue
    if field == '__version_info__':
        value = value.strip('[]()')
        value = '.'.join(x.strip(' \'"') for x in value.split(','))
    else:
        value = value.strip('\'"')
    distmeta[field] = value

setup(
    name='mutalyzer',
    version=distmeta['__version_info__'],
    description='HGVS variant nomenclature checker',
    long_description=long_description,
    author=distmeta['__author__'],
    author_email=distmeta['__contact__'],
    url=distmeta['__homepage__'],
    license='Not distributable',
    platforms=['any'],
    install_requires=install_requires,
    packages=['mutalyzer', 'mutalyzer.parsers', 'mutalyzer.services'],
    include_package_data=True,
    scripts=['bin/mutalyzer',
             'bin/mutalyzer-batchd',
             'bin/mutalyzer-cache-sync',
             'bin/mutalyzer-mapping-update',
             'bin/mutalyzer-mapping-import',
             'bin/mutalyzer-soap-service.wsgi',
             'bin/mutalyzer-json-service.wsgi',
             'bin/mutalyzer-website.wsgi'],
    zip_safe=False
)

# Things not handled by this setup.py:
# - Copy extras/config.example to /etc/mutalyzer/config
# - Database setup
# - Chown /var/log/mutalyzer.log and /var/cache/mutalyzer
# - Copy extras/init.d/mutalyzer-batchd to /etc/init.d/mutalyzer-batchd
# - Copy doc to /usr/share/doc
# Check extras/post-install.sh for these.
