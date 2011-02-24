#!/usr/bin/python

"""
Collection of Serilizable Objects used by the webservice

@requires: soaplib.serializers.primitive.String
@requires: soaplib.serializers.primitive.Integer
@requires: soaplib.serializers.primitive.Array
@requires: soaplib.serializers.clazz.ClassSerializer

@todo: documentation
"""
from soaplib.core.model.primitive import String, Integer, Boolean
from soaplib.core.model.clazz import ClassModel, Array


# Default attributes for soaplib models:
#   nillable = True
#   min_occurs = 0
#   max_occurs = 1
#
# Additional attributes values for String model:
#   min_len = 0
#   max_len = "unbounded"
#   pattern = None


class Mandatory(object):
    """
    This is soaplib.model.primitive.Mandatory, but without min_length=1 for
    the String model.
    """
    String = String(min_occurs=1, nillable=False)
    Integer = Integer(min_occurs=1, nillable=False)
    Boolean = Boolean(min_occurs=1, nillable=False)


# Todo: Use Mandatory.* models in the classmodels below?
# Todo: See if it improves client code if we use Array(_, nillable=False)


class SoapMessage(ClassModel):
    """
    Type of messages used in SOAP method return values.
    """
    __namespace__ = 'http://mutalyzer.nl/2.0/services'

    errorcode = Mandatory.String
    message = Mandatory.String
#SoapMessage


class Mapping(ClassModel) :
    """
    Return type of SOAP method mappingInfo.
    """
    __namespace__ = 'http://mutalyzer.nl/2.0/services'

    startmain = Integer
    startoffset = Integer
    endmain = Integer
    endoffset = Integer
    start_g = Integer
    end_g = Integer
    mutationType = String
    errorcode = Integer
    messages = Array(SoapMessage)
#Mapping


class Transcript(ClassModel) :
    """
    Return type of SOAP method transcriptInfo.
    """
    __namespace__ = 'http://mutalyzer.nl/2.0/services'

    trans_start = Integer
    trans_stop = Integer
    CDS_stop = Integer
#Transcript


class RawVariant(ClassModel):
    """
    Used in MutalyzerOutput data type.
    """
    __namespace__ = 'http://mutalyzer.nl/2.0/services'

    description = Mandatory.String
    visualisation = Mandatory.String
#RawVariant


class MutalyzerOutput(ClassModel) :
    """
    Return type of SOAP method runMutalyzer.
    """
    __namespace__ = 'http://mutalyzer.nl/2.0/services'

    original = String
    mutated = String

    origMRNA = String
    mutatedMRNA= String

    origCDS = String
    newCDS= String

    origProtein = String
    newProtein = String
    altProtein = String

    errors = Integer
    warnings = Integer
    summary = String

    chromDescription = String
    genomicDescription = String
    transcriptDescriptions = Array(String)
    proteinDescriptions = Array(String)

    rawVariants = Array(RawVariant)

    messages = Array(SoapMessage)
#MutalyzerOutput


class TranscriptNameInfo(ClassModel) :
    """
    Return type of SOAP method getGeneAndTranscript.
    """
    __namespace__ = 'http://mutalyzer.nl/2.0/services'

    transcriptName = Mandatory.String
    productName = Mandatory.String
#TranscriptNameInfo


class ExonInfo(ClassModel):
    """
    Used in TranscriptInfo data type.
    """
    __namespace__ = 'http://mutalyzer.nl/2.0/services'

    cStart = Mandatory.String
    gStart = Mandatory.Integer
    cStop = Mandatory.String
    gStop = Mandatory.Integer
#ExonInfo


class ProteinTranscript(ClassModel):
    """
    Used in TranscriptInfo data type.
    """
    __namespace__ = 'http://mutalyzer.nl/2.0/services'

    name = Mandatory.String
    id = Mandatory.String
    product = Mandatory.String
#ProteinTranscript


class TranscriptInfo(ClassModel):
    """
    Used in return type of SOAP method getTranscriptsAndInfo.

    @todo: Decide on 'stop' versus 'end'. Web interface uses 'stop' for
           both trans and CDS. Ivar asked for 'end'. Internally, we have
           trans 'end' and CDS 'stop'.
    """
    __namespace__ = 'http://mutalyzer.nl/2.0/services'

    name = Mandatory.String
    id = Mandatory.String
    product = Mandatory.String

    cTransStart = Mandatory.String
    gTransStart = Mandatory.Integer
    cTransEnd = Mandatory.String
    gTransEnd = Mandatory.Integer
    sortableTransEnd = Mandatory.Integer

    cCDSStart = Mandatory.String
    gCDSStart = Mandatory.Integer
    cCDSStop = Mandatory.String
    gCDSStop = Mandatory.Integer

    locusTag = Mandatory.String
    linkMethod = Mandatory.String

    exons = Array(ExonInfo)

    proteinTranscript = ProteinTranscript
#TranscriptInfo


class CheckSyntaxOutput(ClassModel) :
    """
    Return type of SOAP method checkSyntax.
    """
    __namespace__ = 'http://mutalyzer.nl/2.0/services'

    valid = Mandatory.Boolean
    messages = Array(SoapMessage)
#CheckSyntaxOutput
