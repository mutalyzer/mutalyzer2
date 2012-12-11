#!/usr/bin/env python
"""
Submit a batch job to a Mutalyzer installation.

Usage:
  {command} file [type] [argument]

  file: Batch job input file to upload.
  type: Optional type of the batch job, choose from NameChecker (default),
        SyntaxChecker, PositionConverter, SnpConverter.
  argument: Additional argument. Currently only used if batch_type is
            PositionConverter, denoting the human genome build.

The file is uploaded to the Mutalyzer SOAP API, which is then polled for
completeness of the batch job after which the result is retrieved and printed
to standard output.
"""


from mutalyzer.util import monkey_patch_suds; monkey_patch_suds()

import sys
from suds import WebFault
from suds.client import Client
import time

from mutalyzer.util import format_usage


WSDL_LOCATION = 'http://localhost/mutalyzer/services/?wsdl'

MAX_TRIES = 10
TRY_WAIT = 5


def main(batch_input, batch_type='NameChecker', batch_argument=''):
    """
    Submit a batch job, poll for completeness, and retrieve the result.
    """
    service = Client(WSDL_LOCATION, cache=None).service

    data = open(batch_input, 'rb').read().encode('base64')
    result = service.submitBatchJob(data, batch_type, batch_argument)

    job_id = int(result)

    for _ in range(MAX_TRIES):
        try:
            result = service.getBatchJob(job_id)
            break
        except WebFault:
            # Note that calling getBatchJob *and* monitorBatchJob is a bit
            # superfluous, but we like to illustrate the use of both of them
            # here.
            result = service.monitorBatchJob(job_id)
            sys.stderr.write('Waiting... (%d entries left)\n' % int(result))
            time.sleep(TRY_WAIT)
    else:
        sys.stderr.write('No result after trying %d times.\n' % MAX_TRIES)
        sys.exit(1)

    sys.stdout.write(result.decode('base64'))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print format_usage()
        sys.exit(1)
    main(*sys.argv[1:])
