"""
Mutalyzer webservice HTTP/RPC with JSON response payloads.
"""


from spyne.application import Application
from spyne.protocol.http import HttpRpc
from spyne.protocol.json import JsonObject

import mutalyzer
from mutalyzer.services import rpc


# HTTP/RPC application.
application = Application([rpc.MutalyzerService], tns=mutalyzer.SOAP_NAMESPACE,
                          in_protocol=HttpRpc(), out_protocol=JsonObject(skip_depth=2),
                          name='Mutalyzer')
