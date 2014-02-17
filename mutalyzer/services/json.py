"""
Mutalyzer web service HTTP/RPC with JSON response payloads.
"""


from spyne.application import Application
from spyne.protocol.http import HttpRpc
from spyne.protocol.json import JsonDocument

import mutalyzer
from mutalyzer.services import rpc


#: HTTP/RPC+JSON application.
application = Application([rpc.MutalyzerService], tns=mutalyzer.SOAP_NAMESPACE,
                          in_protocol=HttpRpc(validator='soft'),
                          out_protocol=JsonDocument(),
                          name='Mutalyzer')
