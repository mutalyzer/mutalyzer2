"""
Mutalyzer web service HTTP/RPC with JSON response payloads.
"""


from spyne.application import Application
from spyne.protocol.http import HttpRpc
from spyne.protocol.json import JsonObject

import mutalyzer
from mutalyzer.services import rpc


# HTTP/RPC application.
# Output JSON can be made more consise by specifying skip_depth=2 in the
# JsonObject constructor, but to match the API documentation better (and hence
# have more predictable behaviour) we leave it out.
application = Application([rpc.MutalyzerService], tns=mutalyzer.SOAP_NAMESPACE,
                          in_protocol=HttpRpc(),
                          out_protocol=JsonObject(),
                          name='Mutalyzer')
