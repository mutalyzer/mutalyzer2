"""
Mutalyzer SOAP/1.1 webservice.
"""


from rpclib.application import Application
from rpclib.protocol.soap import Soap11
from rpclib.server.wsgi import WsgiApplication

import mutalyzer
from mutalyzer.services import rpc


# SOAP/1.1 application.
application = Application([rpc.MutalyzerService], tns=mutalyzer.SOAP_NAMESPACE,
                          in_protocol=Soap11(), out_protocol=Soap11(),
                          name='Mutalyzer')


# Below we define WSGI applications for use with e.g. Apache/mod_wsgi.
# Note: We would like to create the wsgi.Application instance only in the
#     bin/mutalyzer-webservice.wsgi script, but unfortunately this breaks the
#     get_interface_document method of rpclib which we use to generate API
#     documentation in website.py.
wsgi_application = WsgiApplication(application)