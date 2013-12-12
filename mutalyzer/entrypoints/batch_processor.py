"""
Mutalyzer batch processor.

.. todo: Get rid of ugly exception logging.
.. todo: Reload configuration without restarting (for example, on SIGHUP).
"""


import argparse
import signal
import sys
import time
import traceback

from .. import config
from .. import Db
from .. import Scheduler


def process():
    """
    Run forever in a loop processing scheduled batch jobs.
    """
    database = Db.Batch()
    counter = Db.Counter()
    scheduler = Scheduler.Scheduler(database)

    def handle_exit(signum, stack_frame):
        if scheduler.stopped():
            sys.stderr.write('mutalyzer-batchd: Terminated\n')
            sys.exit(1)
        if signum == signal.SIGINT:
            sys.stderr.write('mutalyzer-batchd: Hitting Ctrl+C again will '
                             'terminate any running job!\n')
        scheduler.stop()

    signal.signal(signal.SIGTERM, handle_exit)
    signal.signal(signal.SIGINT, handle_exit)

    while not scheduler.stopped():
        # Process batch jobs. This process() method runs while there
        # exist jobs to run.
        try:
            scheduler.process(counter)
        except Exception as e:
            f = open('/tmp/batcherror.log', 'a+')
            f.write('Error (%s): %s\n' % (type(e), str(e)))
            f.write('%s\n\n' % repr(traceback.format_exc()))
            f.flush()
            f.close()
        if scheduler.stopped():
            break
        # Wait a bit and process any possible new jobs.
        time.sleep(1)

    sys.stderr.write('mutalyzer-batchd: Graceful shutdown\n')
    sys.exit(0)


def main():
    """
    Command line interface to the batch processor.
    """
    parser = argparse.ArgumentParser(
        description='Mutalyzer batch processor.',
        epilog='The process can be shutdown gracefully by sending a SIGINT '
        '(Ctrl+C) or SIGTERM signal.')

    parser.parse_args()
    process()


if __name__ == '__main__':
    main()
