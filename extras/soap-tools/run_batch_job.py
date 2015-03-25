#!/usr/bin/env python
"""
Mutalyzer batch job automation.


Copyright (c) 2015 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2015 Jeroen F.J. Laros <j.f.j.laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""


import argparse
import time

from suds.client import Client


URL = 'http://127.0.0.1:8081/?wsdl'
job_types = ["NameChecker", "SyntaxChecker", "PositionConverter",
             "SnpConverter"]
RETRY_WAIT = 1


def run_batchjob(input_handle, output_handle, job_type, build="hg19"):
    """
    Run a Mutalyzer batch job.

    :arg stream input_handle: Open readable handle to a batch job file.
    :arg stream output_handle: Opren writable handle for the results.
    :arg str job_type: Type of the job.
    :arg str build: Optional build for the Position Converter.
    """
    client = Client(URL, cache=None)

    batchfile = input_handle.read().encode("base64")
    job_id = client.service.submitBatchJob(batchfile, job_type, build)
    while client.service.monitorBatchJob(job_id):
        time.sleep(RETRY_WAIT)
    result = client.service.getBatchJob(job_id).decode("base64")
    output_handle.write(result)


def main():
    """
    Main entry point.
    """
    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])

    parser.add_argument("-v", action='version', version="0.1")
    parser.add_argument("input_handle", metavar="INPUT",
        type=argparse.FileType("r"), help="input file")
    parser.add_argument("output_handle", metavar="OUTPUT",
        type=argparse.FileType("w"), help="output file")
    parser.add_argument("job_type", metavar="TYPE", choices=job_types,
        help="batch job type ({})".format(", ".join(job_types)))
    parser.add_argument("-b", dest="build", type=str,
        help="genome build (only valid for the Position Converter)")

    args = parser.parse_args()
    run_batchjob(args.input_handle, args.output_handle, args.job_type)


if __name__ == "__main__":
    main()
