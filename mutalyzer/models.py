"""
Collection of serilizable objects used by the SOAP web service. They extend
from the spyne ClassModel.

Default attributes for the spyne ClassModel:
- nillable = True
- min_occurs = 0
- max_occurs = 1

Additional attributes values for the spyne String model:
- min_len = 0
- max_len = 'unbounded'
- pattern = None

@todo: Use Mandatory.* models in the ClassModel extensions?
@todo: See if it improves client code if we use Array(_, nillable=False).
"""


from spyne.model.primitive import String, Integer, Boolean, DateTime
from spyne.model.binary import ByteArray
from spyne.model.complex import ComplexModel, Array

from mutalyzer import SOAP_NAMESPACE


class Mandatory(object):
    """
    This is spyne.model.primitive.Mandatory, but without min_length=1 for
    the String model.
    """
    String = String(type_name='mandatory_string', min_occurs=1, nillable=False)
    Integer = Integer(type_name='mandatory_integer', min_occurs=1, nillable=False)
    Boolean = Boolean(type_name='mandatory_boolean', min_occurs=1, nillable=False)
    DateTime = DateTime(type_name='mandatory_date_time', min_occurs=1, nillable=False)
    ByteArray = ByteArray(type_name='mandatory_byte_array', min_occurs=1, nillable=False)
#Mandatory


class SoapMessage(ComplexModel):
    """
    Type of messages used in SOAP method return values.
    """
    __namespace__ = SOAP_NAMESPACE

    errorcode = Mandatory.String
    message = Mandatory.String
#SoapMessage


class Mapping(ComplexModel):
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


class Transcript(ComplexModel):
    """
    Return type of SOAP method transcriptInfo.
    """
    __namespace__ = SOAP_NAMESPACE

    trans_start = Integer
    trans_stop = Integer
    CDS_stop = Integer
#Transcript


class RawVariant(ComplexModel):
    """
    Used in MutalyzerOutput data type.
    """
    __namespace__ = SOAP_NAMESPACE

    description = Mandatory.String
    visualisation = Mandatory.String
#RawVariant


class RawVar(ComplexModel):
    """
    Used in MutalyzerOutput data type.
    """
    __namespace__ = SOAP_NAMESPACE

    DNA = Mandatory.Boolean
    start = Mandatory.Integer
    start_offset = Mandatory.Integer
    end = Mandatory.Integer
    end_offset = Mandatory.Integer
    type = Mandatory.String
    deleted = Mandatory.String
    inserted = Mandatory.String
    shift = Mandatory.Integer
    startAA = Mandatory.String
    endAA = Mandatory.String
    term = Mandatory.Integer
    hgvs = Mandatory.String
    hgvsLength = Mandatory.Integer
#RawVar


class Allele(ComplexModel):
    """
    Used in MutalyzerOutput data type.
    """
    __namespace__ = SOAP_NAMESPACE

    description = Mandatory.String
    allele = Array(RawVar)
#Allele


class ExonInfo(ComplexModel):
    """
    Used in TranscriptInfo and MutalyzerOutput data types.
    """
    __namespace__ = SOAP_NAMESPACE

    cStart = Mandatory.String
    gStart = Mandatory.Integer
    chromStart = Integer
    cStop = Mandatory.String
    gStop = Mandatory.Integer
    chromStop = Integer
#ExonInfo


class MutalyzerOutput(ComplexModel):
    """
    Return type of SOAP method runMutalyzer.
    """
    __namespace__ = SOAP_NAMESPACE

    referenceId = Mandatory.String
    sourceId = Mandatory.String
    sourceAccession = String
    sourceVersion = String
    sourceGi = String
    molecule = Mandatory.String

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

    exons = Array(ExonInfo)

    rawVariants = Array(RawVariant)

    messages = Array(SoapMessage)
#MutalyzerOutput


class TranscriptNameInfo(ComplexModel):
    """
    Return type of SOAP method getGeneAndTranscript.
    """
    __namespace__ = SOAP_NAMESPACE

    transcriptName = Mandatory.String
    productName = Mandatory.String
#TranscriptNameInfo


class ProteinTranscript(ComplexModel):
    """
    Used in TranscriptInfo data type.
    """
    __namespace__ = SOAP_NAMESPACE

    name = Mandatory.String
    id = Mandatory.String
    product = Mandatory.String
#ProteinTranscript


class TranscriptInfo(ComplexModel):
    """
    Used in return type of SOAP method getTranscriptsAndInfo.

    @todo: Decide on 'stop' versus 'end'. Web interface uses 'stop' for both
        trans and CDS. Ivar asked for 'end'. Internally, we have trans 'end'
        and CDS 'stop'.
    @todo: We should really also provide the chromosome (or its accession
        number) next to the chromosomal positions, if available.
    """
    __namespace__ = SOAP_NAMESPACE

    name = Mandatory.String
    id = Mandatory.String
    product = Mandatory.String

    cTransStart = Mandatory.String
    gTransStart = Mandatory.Integer
    chromTransStart = Integer
    cTransEnd = Mandatory.String
    gTransEnd = Mandatory.Integer
    chromTransEnd = Integer
    sortableTransEnd = Mandatory.Integer

    cCDSStart = Mandatory.String
    gCDSStart = Mandatory.Integer
    chromCDSStart = Integer
    cCDSStop = Mandatory.String
    gCDSStop = Mandatory.Integer
    chromCDSStop = Integer

    locusTag = Mandatory.String
    linkMethod = Mandatory.String

    exons = Array(ExonInfo)

    proteinTranscript = ProteinTranscript
#TranscriptInfo


class TranscriptMappingInfo(ComplexModel):
    """
    Used in return type of SOAP method getTranscriptsRange.
    """
    __namespace__ = SOAP_NAMESPACE

    name = Mandatory.String
    version = Mandatory.Integer
    gene = Mandatory.String
    protein = Mandatory.String
    orientation = Mandatory.String

    start = Mandatory.Integer
    stop = Mandatory.Integer

    cds_start = Mandatory.Integer
    cds_stop = Mandatory.Integer
#TranscriptMappingInfo


class CheckSyntaxOutput(ComplexModel):
    """
    Return type of SOAP method checkSyntax.
    """
    __namespace__ = SOAP_NAMESPACE

    valid = Mandatory.Boolean
    messages = Array(SoapMessage)
#CheckSyntaxOutput


class InfoOutput(ComplexModel):
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


class CacheEntry(ComplexModel):
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
    cached = String
#CacheEntry
