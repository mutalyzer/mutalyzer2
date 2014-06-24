"""
WSGI interface to the Mutalyzer SOAP webservice.

Example Apache/mod_wsgi configuration:

.. code-block:: apache

    WSGIScriptAlias /soap /usr/local/bin/mutalyzer-service-soap

Be sure to have this line first if you also define a / alias, like this:

.. code-block:: apache

    WSGIScriptAlias /soap /usr/local/bin/mutalyzer-service-soap
    WSGIScriptAlias / /usr/local/bin/mutalyzer-website

You can also use the built-in HTTP server by running this file directly.
"""


import argparse
import logging
import sys

from wsgiref.simple_server import make_server
from spyne.server.wsgi import WsgiApplication

from ..config import settings
from ..services import soap


# Setup logging.
log_level = logging.INFO if settings.DEBUG else logging.ERROR
logging.basicConfig(level=log_level)
logging.getLogger('spyne.protocol.xml').setLevel(log_level)


#: WSGI application instance.
application = WsgiApplication(soap.application)


def debugserver(port):
    """
    Run the webservice with the Python built-in HTTP server.
    """
    sys.stderr.write('Listening on http://localhost:%d/\n' % port)
    sys.stderr.write('WDSL file is at http://localhost:%d/?wsdl\n' % port)
    make_server('localhost', port, application).serve_forever()


def main():
    """
    Command-line interface to the SOAP webservice.
    """
    parser = argparse.ArgumentParser(
        description='Mutalyzer SOAP webservice.')
    parser.add_argument(
        '-p', '--port', metavar='NUMBER', dest='port', type=int,
        default=8081, help='port to run the webservice on (default: 8081)')

    args = parser.parse_args()
    debugserver(args.port)


if __name__ == '__main__':
    main()
