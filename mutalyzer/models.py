"""
Collection of serilizable objects used by the SOAP webservice. They extend
from the soaplib ClassModel.

Default attributes for the soaplib ClassModel:
- nillable = True
- min_occurs = 0
- max_occurs = 1

Additional attributes values for the soaplib String model:
- min_len = 0
- max_len = 'unbounded'
- pattern = None

@todo: Use Mandatory.* models in the ClassModel extensions?
@todo: See if it improves client code if we use Array(_, nillable=False).
"""


from soaplib.core.model.primitive import String, Integer, Boolean, DateTime
from soaplib.core.model.clazz import ClassModel, Array

from mutalyzer import SOAP_NAMESPACE


class Mandatory(object):
    """
    This is soaplib.model.primitive.Mandatory, but without min_length=1 for
    the String model.
    """
    String = String(min_occurs=1, nillable=False)
    Integer = Integer(min_occurs=1, nillable=False)
    Boolean = Boolean(min_occurs=1, nillable=False)
    DateTime = DateTime(min_occurs=1, nillable=False)
#Mandatory


class SoapMessage(ClassModel):
    """
    Type of messages used in SOAP method return values.
    """
    __namespace__ = SOAP_NAMESPACE

    errorcode = Mandatory.String
    message = Mandatory.String
#SoapMessage


class Mapping(ClassModel):
    """
    Return type of SOAP method mappingInfo.
    """
    __namespace__ = SOAP_NAMESPACE

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


class Transcript(ClassModel):
    """
    Return type of SOAP method transcriptInfo.
    """
    __namespace__ = SOAP_NAMESPACE

    trans_start = Integer
    trans_stop = Integer
    CDS_stop = Integer
#Transcript


class RawVariant(ClassModel):
    """
    Used in MutalyzerOutput data type.
    """
    __namespace__ = SOAP_NAMESPACE

    description = Mandatory.String
    visualisation = Mandatory.String
#RawVariant


class MutalyzerOutput(ClassModel):
    """
    Return type of SOAP method runMutalyzer.
    """
    __namespace__ = SOAP_NAMESPACE

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


class TranscriptNameInfo(ClassModel):
    """
    Return type of SOAP method getGeneAndTranscript.
    """
    __namespace__ = SOAP_NAMESPACE

    transcriptName = Mandatory.String
    productName = Mandatory.String
#TranscriptNameInfo


class ExonInfo(ClassModel):
    """
    Used in TranscriptInfo data type.
    """
    __namespace__ = SOAP_NAMESPACE

    cStart = Mandatory.String
    gStart = Mandatory.Integer
    cStop = Mandatory.String
    gStop = Mandatory.Integer
#ExonInfo


class ProteinTranscript(ClassModel):
    """
    Used in TranscriptInfo data type.
    """
    __namespace__ = SOAP_NAMESPACE

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
    __namespace__ = SOAP_NAMESPACE

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


class CheckSyntaxOutput(ClassModel):
    """
    Return type of SOAP method checkSyntax.
    """
    __namespace__ = SOAP_NAMESPACE

    valid = Mandatory.Boolean
    messages = Array(SoapMessage)
#CheckSyntaxOutput


class InfoOutput(ClassModel):
    """
    Return type of SOAP method info.
    """
    __namespace__ = SOAP_NAMESPACE

    version = String
    versionParts = Array(String)
    releaseDate = String
    nomenclatureVersion = String
    nomenclatureVersionParts = Array(String)
    serverName = String
    contactEmail = String
#InfoOutput


class CacheEntry(ClassModel):
    """
    Used in getCache SOAP method.
    """
    __namespace__ = SOAP_NAMESPACE

    name = Mandatory.String
    gi = String
    hash = Mandatory.String
    chromosomeName = String
    chromosomeStart = Integer
    chromosomeStop = Integer
    chromosomeOrientation = Integer
    url = String
    created = Mandatory.DateTime
    cached = Mandatory.Boolean
#CacheEntry
