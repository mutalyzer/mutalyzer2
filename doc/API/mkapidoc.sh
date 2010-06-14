#!/bin/sh

epydoc --config api.conf
mv api/api.pdf .
rm -rf api/
