"""
Mutalyzer SOAP/1.1 web service.
"""


from __future__ import unicode_literals

from spyne.application import Application
from spyne.protocol.soap import Soap11

import mutalyzer
from mutalyzer.services import rpc


#: SOAP/1.1 application.
application = Application([rpc.MutalyzerService], tns=mutalyzer.SOAP_NAMESPACE,
                          in_protocol=Soap11(validator='lxml'),
                          out_protocol=Soap11(),
                          name='Mutalyzer')
