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


from __future__ import unicode_literals

import argparse
import logging
import sys

from wsgiref.simple_server import make_server
from spyne.server.wsgi import WsgiApplication

from . import _cli_string, _ReverseProxied
from ..config import settings
from ..services import soap


# Setup logging.
log_level = logging.INFO if settings.DEBUG else logging.ERROR
logging.basicConfig(level=log_level)
logging.getLogger('spyne.protocol.xml').setLevel(log_level)


#: WSGI application instance.
application = WsgiApplication(soap.application)

# By default, the first request will trigger a build of the WSDL document,
# using the context (service location) from that request. For example, if the
# first request is on `http://localhost/` and subsequent requests are on
# `https://mutalyzer.nl/services/`, they will not have a valid WSDL document.
#
# This is actually what we do on our production infrastructure, where the
# service is tested (on localhost) after it has been started.
#
# The fix is to force a build of the WSDL document and specifying the location
# to use.
#
# http://spyne.io/docs/2.10/reference/server.html#spyne.server.wsgi.WsgiApplication
if settings.SOAP_WSDL_URL:
    application.doc.wsdl11.build_interface_document(settings.SOAP_WSDL_URL)

if settings.REVERSE_PROXIED:
    application = _ReverseProxied(application)


def debugserver(host, port):
    """
    Run the webservice with the Python built-in HTTP server.
    """
    sys.stderr.write('Listening on http://%s:%d/\n' % (host, port))
    sys.stderr.write('WDSL file is at http://%s:%d/?wsdl\n' % (host, port))
    make_server(host, port, application).serve_forever()


def main():
    """
    Command-line interface to the SOAP webservice.
    """
    parser = argparse.ArgumentParser(
        description='Mutalyzer SOAP webservice.')
    parser.add_argument(
        '-H', '--host', metavar='HOSTNAME', type=_cli_string, dest='host',
        default='127.0.0.1', help='hostname to listen on (default: '
        '127.0.0.1; specify 0.0.0.0 to listen on all hostnames)')
    parser.add_argument(
        '-p', '--port', metavar='PORT', dest='port', type=int,
        default=8081, help='port to listen on (default: 8081)')

    args = parser.parse_args()
    debugserver(args.host, args.port)


if __name__ == '__main__':
    main()
