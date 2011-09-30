#!/usr/bin/env python

"""
Search log for bugs.

Finds occurrences of 'Received' in the log file that are not followed by an
occurrence of 'Finished'. These are probably runs of the namechecker that
crashed.
"""


import os
from mutalyzer.config import Config


config = Config()
handle = open(config.Output.log, 'r')

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
