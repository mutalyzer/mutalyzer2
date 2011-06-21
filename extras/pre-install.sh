#!/bin/bash

# Pre-install script for Mutalyzer on Debian or Debian-like systems. Run
# before the setuptools installation (python setup.py install).
#
# Notice: The definitions in this file are quite specific to the standard
# Mutalyzer environment. This consists of a Debian stable (Squeeze) system
# with Apache and Mutalyzer using its mod_wsgi module. Debian conventions are
# used throughout. See the README file for more information.
#
# Usage (from the source root directory):
#   sudo bash extras/pre-install.sh

set -e

echo "Installing packages with apt"

apt-get install -y \
  mysql-server \
  python \
  python-mysqldb \
  python-biopython \
  python-pyparsing \
  python-configobj \
  python-simpletal \
  python-soappy \
  python-magic \
  python-psutil \
  python-xlrd \
  python-daemon \
  python-webpy \
  python-webtest \
  python-nose \
  apache2 \
  libapache2-mod-wsgi \
  python-setuptools \
  git-core

echo "Installing latest soaplib from git master"

mkdir -p /tmp/mutalyzer-install
pushd /tmp/mutalyzer-install

git clone https://github.com/soaplib/soaplib.git
cd soaplib
python setup.py install

popd
rm -Rf /tmp/mutalyzer-install

echo "Installing suds using easy_install"

easy_install suds

echo "kthxbye"
