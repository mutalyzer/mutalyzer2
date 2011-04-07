#!/bin/bash

# Pre-install script for Mutalyzer. Run before the setuptools installation
# (python setup.py install).
#
# Usage (from the source root directory):
#   sudo bash extras/pre-install.sh

set -e

echo "Installing packages with apt"

apt-get install \
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
  python-suds \
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
sudo python setup.py install

popd
rm -Rf /tmp/mutalyzer-install

echo "kthxbye"
