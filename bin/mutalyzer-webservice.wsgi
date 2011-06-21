#!/usr/bin/env python

"""
WSGI interface to the Mutalyzer SOAP webservice.

The WSGI interface is exposed through the module variable 'application'.

Example Apache/mod_wsgi configuration:

  WSGIScriptAlias /services /usr/local/bin/mutalyzer-webservice.wsgi

Be sure to have this line first if you also define a / alias, like this:

  WSGIScriptAlias /services /usr/local/bin/mutalyzer-webservice.wsgi
  WSGIScriptAlias / /usr/local/bin/mutalyzer-website.wsgi

You can also use the built-in HTTP server by running this file directly.

To start the built-in HTTP server on port 8081:

  /usr/local/bin/mutalyzer-webservice.wsgi 8081

@todo: Do we really use namespaces correctly?
@todo: For some reason, the server exposes its location including ?wsdl.
@todo: More thourough input checking. The @soap decorator does not do any
       kind of strictness checks on the input. For example, in
       transcriptInfo, the build argument must really be present. (Hint:
       use __checkBuild.)
@todo: The mutalyzer.config.Config object can just be instantiated once
       and we should not create it on every request.
"""


import sys
from wsgiref.simple_server import make_server
from soaplib.core.server import wsgi
from mutalyzer import webservice


DEFAULT_PORT = 8081


application = wsgi.Application(webservice.soap_application)


if __name__ == '__main__':
    port = DEFAULT_PORT
    if len(sys.argv) > 1:
        try:
            port = int(sys.argv[1])
        except ValueError:
            print 'Not a valid port number: %s' % sys.argv[1]
            sys.exit(1)
    print 'Listening to http://localhost:%d/' % port
    print 'WDSL file is at http://localhost:%d/?wsdl' % port
    make_server('localhost', port, application).serve_forever()
