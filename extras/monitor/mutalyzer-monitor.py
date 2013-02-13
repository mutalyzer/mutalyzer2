#!/usr/bin/env python
"""
Monitor a Mutalyzer server by some basic uptime checks.

Run with no arguments for usage info. The script returns a 0 exit status
on success, 1 on failure. It writes the error to standard error.

Example usage:
  ./mutalyzer-monitor.py && echo ok || echo problem

Currently implemented checks:
- Website homepage exists.
- Name checker SOAP web service can be called.
- Batch name checker can be called.
"""


import argparse
import logging
import sys
import time
import urllib2

from suds import WebFault
from suds.client import Client


DEFAULT_MUTALYZER_URL = 'https://mutalyzer.nl'
BATCH_MAX_POLL = 25
BATCH_POLL_WAIT = 5


logging.getLogger('suds.client').setLevel(logging.CRITICAL)


def main(mutalyzer_url):
    """
    Run all checks.
    """
    checks = check_website, check_soap, check_batch
    try:
        [check(mutalyzer_url) for check in checks]
    except BaseException as e:
        sys.stderr.write(str(e) + '\n')
        sys.exit(1)


def check_website(mutalyzer_url):
    """
    Check if website homepage exists.
    """
    response = urllib2.urlopen(mutalyzer_url)
    html = response.read()
    assert 'Welcome to the Mutalyzer web site' in html


def check_soap(mutalyzer_url):
    """
    Check if the name checker SOAP web service can be called.
    """
    wsdl_url = mutalyzer_url + '/services/?wsdl'
    service = Client(wsdl_url, cache=None).service

    result = service.runMutalyzer('AB026906.1:c.274G>T')
    assert result.genomicDescription == 'AB026906.1:g.7872G>T'


def check_batch(mutalyzer_url):
    """
    Check if the batch name checker can be called.
    """
    wsdl_url = mutalyzer_url + '/services/?wsdl'
    service = Client(wsdl_url, cache=None).service

    variants = ['AB026906.1(SDHD):g.7872G>T',
                'NM_003002.1:c.3_4insG',
                'AL449423.14(CDKN2A_v002):c.5_400del']
    data = '\n'.join(variants).encode('base64')
    result = service.submitBatchJob(data, 'NameChecker')

    job_id = int(result)

    for _ in range(BATCH_MAX_POLL):
        try:
            result = service.getBatchJob(job_id)
            break
        except WebFault:
            time.sleep(BATCH_POLL_WAIT)
    else:
        sys.exit(1)

    response = result.decode('base64')
    assert len([_ for line in response.split('\n') if line.strip()]) - 1 == len(variants)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('mutalyzer_url', metavar='MUTALYZER_URL',
                        default=DEFAULT_MUTALYZER_URL, nargs='?',
                        help='Mutalyzer instance URL (default: %(default)s')
    args = parser.parse_args()
    main(args.mutalyzer_url)
