"""
Mutalyzer webservice HTTP/RPC with JSON response payloads.
"""


from mutalyzer.util import monkey_patch_spyne; monkey_patch_spyne()

from spyne.application import Application
from spyne.protocol.http import HttpRpc
from spyne.protocol.json import JsonObject

import mutalyzer
from mutalyzer.services import rpc


# HTTP/RPC application.
# Note that we originally provided skip_depth=2 to the JsonObject constructor
# to get rid of some annoying wrappers around the json results. However, this
# breaks methods returning a primitive datatype (e.g. just a string). Might
# file a bug about this on spyne.
application = Application([rpc.MutalyzerService], tns=mutalyzer.SOAP_NAMESPACE,
                          in_protocol=HttpRpc(), out_protocol=JsonObject(),
                          name='Mutalyzer')
