#!/usr/bin/env python
"""
WSGI interface to the Mutalyzer HTTP/RPC+JSON webservice.

The WSGI interface is exposed through the module variable 'application'.

Example Apache/mod_wsgi configuration:

  WSGIScriptAlias /json /usr/local/bin/mutalyzer-json-service.wsgi

Be sure to have this line first if you also define a / alias, like this:

  WSGIScriptAlias /services /usr/local/bin/mutalyzer-json-service.wsgi
  WSGIScriptAlias / /usr/local/bin/mutalyzer-website.wsgi

You can also use the built-in HTTP server by running this file directly.

To start the built-in HTTP server on port 8082:

  /usr/local/bin/mutalyzer-json-service.wsgi 8082
"""


import sys
from wsgiref.simple_server import make_server
from spyne.server.wsgi import WsgiApplication
from mutalyzer.services import json


DEFAULT_PORT = 8082


application = WsgiApplication(json.application)


if __name__ == '__main__':
    port = DEFAULT_PORT
    if len(sys.argv) > 1:
        try:
            port = int(sys.argv[1])
        except ValueError:
            print 'Not a valid port number: %s' % sys.argv[1]
            sys.exit(1)
    print 'Listening at http://localhost:%d/' % port
    make_server('localhost', port, application).serve_forever()
