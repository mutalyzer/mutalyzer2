"""
Mutalyzer webservice HTTP/RPC with JSON response payloads.
"""


from rpclib.application import Application
from rpclib.protocol.http import HttpRpc
from rpclib.protocol.json import JsonObject

import mutalyzer
from mutalyzer.services import rpc


# HTTP/RPC application.
application = Application([rpc.MutalyzerService], tns=mutalyzer.SOAP_NAMESPACE,
                          in_protocol=HttpRpc(), out_protocol=JsonObject(skip_depth=2),
                          name='Mutalyzer')
