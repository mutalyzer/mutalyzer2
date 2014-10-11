#!/usr/bin/env python

"""
Search log for bugs.

Finds occurrences of 'Received' in the log file that are not followed by an
occurrence of 'Finished'. These are probably runs of the namechecker that
crashed.
"""


from __future__ import unicode_literals

import os
from mutalyzer import config


handle = open(config.get('log'), 'r')

scanning = False
line = handle.readline()

while line:
    if not scanning:
        if ' Received ' in line:
            message = line
            scanning = True
    else:
        if ' Received ' in line:
            print message,
            scanning = False
        if ' Finished ' in line:
            scanning = False
    line = handle.readline()

handle.close()
