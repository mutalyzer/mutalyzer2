#!/bin/bash

# Post-install script for Mutalyzer. Run after the setuptools installation
# (python setup.py install).
#
# Usage (from the source root directory):
#   sudo bash extras/post-install.sh

# Todo:
# - Database setup
# - Copy doc to /usr/share/doc

if [ ! -f /etc/mutalyzer/config ]; do
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

if [ ! -f /etc/init.d/mutalyzer-batchd ]; do
    cp extras/init.d/mutalyzer-batchd /etc/init.d/mutalyzer-batchd
    chmod u=rwx,go=rx /etc/init.d/mutalyzer-batchd
fi
