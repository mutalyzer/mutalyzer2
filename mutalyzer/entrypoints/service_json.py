"""
WSGI interface to the Mutalyzer HTTP/RPC+JSON webservice.

Example Apache/mod_wsgi configuration:

.. code-block:: apache

    WSGIScriptAlias /json /usr/local/bin/mutalyzer-service-json

Be sure to have this line first if you also define a / alias, like this:

.. code-block:: apache

    WSGIScriptAlias /json /usr/local/bin/mutalyzer-service-json
    WSGIScriptAlias / /usr/local/bin/mutalyzer-website

You can also use the built-in HTTP server by running this file directly.
"""


import argparse
import logging
import sys

from wsgiref.simple_server import make_server
from spyne.server.wsgi import WsgiApplication

from . import _ReverseProxied
from ..config import settings
from ..services import json


# Setup logging.
log_level = logging.INFO if settings.DEBUG else logging.ERROR
logging.basicConfig(level=log_level)
logging.getLogger('spyne.protocol.xml').setLevel(log_level)


#: WSGI application instance.
application = WsgiApplication(json.application)
if settings.REVERSE_PROXIED:
    application = _ReverseProxied(application)


def debugserver(port):
    """
    Run the webservice with the Python built-in HTTP server.
    """
    sys.stderr.write('Listening on http://localhost:%d/\n' % port)
    make_server('localhost', port, application).serve_forever()


def main():
    """
    Command-line interface to the HTTP/RPC+JSON webservice.
    """
    parser = argparse.ArgumentParser(
        description='Mutalyzer HTTP/RPC+JSON webservice.')
    parser.add_argument(
        '-p', '--port', metavar='NUMBER', dest='port', type=int,
        default=8082, help='port to run the webservice on (default: 8082)')

    args = parser.parse_args()
    debugserver(args.port)


if __name__ == '__main__':
    main()
