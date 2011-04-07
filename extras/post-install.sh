#!/bin/bash

# Post-install script for Mutalyzer. Run after the setuptools installation
# (python setup.py install).
#
# Usage (from the source root directory):
#   sudo bash extras/post-install.sh

# Todo:
# - Database setup
# - Apache setup
# - Copy doc to /usr/share/doc

set -e

# The 'cd /' is a hack to prevent the mutalyzer package under the current
# directory to be used.
PACKAGE_ROOT=$(cd / && python -c 'import mutalyzer; print mutalyzer.package_root()')
BIN_BATCHD=$(which mutalyzer-batchd)

if [ ! -e /etc/mutalyzer/config ]; then
    mkdir -p /etc/mutalyzer
    cp extras/config.example /etc/mutalyzer/config
    chmod -R u=rwX,go=rX /etc/mutalyzer
fi

touch /var/log/mutalyzer.log
chown www-data:www-data /var/log/mutalyzer.log
chmod u=rw,go=r /var/log/mutalyzer.log

mkdir -p /var/cache/mutalyzer
chown -R www-data:www-data /var/cache/mutalyzer
chmod -R u=rwX,go=rX /var/cache/mutalyzer

if [ ! -e /etc/init.d/mutalyzer-batchd ]; then
    cp extras/init.d/mutalyzer-batchd /etc/init.d/mutalyzer-batchd
    sed -i -e "s@<MUTALYZER_PACKAGE_ROOT>@${PACKAGE_ROOT}@g" -e "s@<MUTALYZER_BIN_BATCHD>@${BIN_BATCHD}@g" /etc/init.d/mutalyzer-batchd
    chmod u=rwx,go=rx /etc/init.d/mutalyzer-batchd
fi

if [ ! -e /etc/apache2/conf.d/mutalyzer.conf ]; then
    cp extras/apache/mutalyzer.conf /etc/apache2/conf.d/mutalyzer.conf
    sed -i -e "s@<MUTALYZER_PACKAGE_ROOT>@${PACKAGE_ROOT}@g" -e "s@<MUTALYZER_BIN_BATCHD>@${BIN_BATCHD}@g" /etc/apache2/conf.d/mutalyzer.conf
    chmod u=rw,go=r /etc/apache2/conf.d/mutalyzer.conf
fi
