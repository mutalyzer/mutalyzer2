import os
import sys
from setuptools import setup

if sys.version_info < (2, 7):
    raise Exception('Mutalyzer requires Python 2.7 or higher.')

install_requires = []

try:
    with open('README.rst') as readme:
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
    license='AGPL',
    platforms=['any'],
    install_requires=install_requires,
    packages=['mutalyzer',
              'mutalyzer.config',
              'mutalyzer.db',
              'mutalyzer.entrypoints',
              'mutalyzer.parsers',
              'mutalyzer.services',
              'mutalyzer.website'],
    include_package_data=True,
    entry_points={'console_scripts': [
        'mutalyzer = mutalyzer.entrypoints.mutalyzer:main',
        'mutalyzer-admin = mutalyzer.entrypoints.admin:main',
        'mutalyzer-batch-processor = mutalyzer.entrypoints.batch_processor:main',
        'mutalyzer-service-json = mutalyzer.entrypoints.service_json:main',
        'mutalyzer-service-soap = mutalyzer.entrypoints.service_soap:main',
        'mutalyzer-website = mutalyzer.entrypoints.website:main']},
    zip_safe=False,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    keywords='bioinformatics'
)
