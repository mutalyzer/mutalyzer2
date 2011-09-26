#!/bin/bash

# Post-upgrade script for Mutalyzer. Run after the setuptools installation
# (python setup.py install).
#
# Notice: The definitions in this file are quite specific to the standard
# Mutalyzer environment. This consists of a Debian stable (Squeeze) system
# with Apache and Mutalyzer using its mod_wsgi module. Debian conventions are
# used throughout. See the README file for more information.
#
# Usage (from the source root directory):
#   sudo bash extras/post-upgrade.sh

set -e

COLOR_INFO='\033[32m'
COLOR_WARNING='\033[33m'
COLOR_ERROR='\033[31m'
COLOR_END='\033[0m'

# The 'cd /' is a hack to prevent the mutalyzer package under the current
# directory to be used.
PACKAGE_ROOT=$(cd / && python -c 'import mutalyzer; print mutalyzer.package_root()')
BIN_WEBSITE=$(which mutalyzer-website.wsgi)
BIN_WEBSERVICE=$(which mutalyzer-webservice.wsgi)

if [ ! -e /var/www/mutalyzer ]; then
    mkdir -p /var/www/mutalyzer
fi

if [ -e /var/www/mutalyzer/base ]; then
    echo "Removing /var/www/mutalyzer/base"
    rm /var/www/mutalyzer/base
fi

echo -e "${COLOR_INFO}Symlinking /var/www/mutalyzer/base to $PACKAGE_ROOT/templates/base${COLOR_END}"
ln -s $PACKAGE_ROOT/templates/base /var/www/mutalyzer/base

echo "Running any needed migrations"
for MIGRATION in extras/migrations/*.migration; do
    echo "Checking migration $(basename $MIGRATION)"
    $MIGRATION migrate
done

echo -e "${COLOR_INFO}Assuming mod_wsgi daemon mode, not restarting Apache${COLOR_END}"
#/etc/init.d/apache2 restart

echo "Touching WSGI entry to reload application"
touch $BIN_WEBSITE
touch $BIN_WEBSERVICE

echo "Restarting Mutalyzer batch daemon"
/etc/init.d/mutalyzer-batchd restart
