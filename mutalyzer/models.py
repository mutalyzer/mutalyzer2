"""
Collection of serilizable objects used by the SOAP web service. They extend
from the Spyne model classes.

@todo: Move all these models to the mutalyzer.services package and refactor
  all uses of them in other places. The SOAP_NAMESPACE variable can then also
  be moved there.
"""


from __future__ import unicode_literals

from spyne.model.primitive import Integer, Boolean, DateTime, Unicode
from spyne.model.binary import ByteArray
from spyne.model.complex import ComplexModel, Array

from mutalyzer import SOAP_NAMESPACE


class Mandatory(object):
    """
    This is spyne.model.primitive.Mandatory, but without min_length=1 for
    the Unicode model.
    """
    Unicode = Unicode(type_name='mandatory_unicode', min_occurs=1, nillable=False)
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

    errorcode = Mandatory.Unicode
    message = Mandatory.Unicode
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
    mutationType = Unicode
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

    description = Mandatory.Unicode
    visualisation = Mandatory.Unicode
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
    type = Mandatory.Unicode
    deleted = Mandatory.Unicode
    inserted = Mandatory.Unicode
    shift = Mandatory.Integer
    startAA = Mandatory.Unicode
    endAA = Mandatory.Unicode
    term = Mandatory.Integer
    hgvs = Mandatory.Unicode
    hgvsLength = Mandatory.Integer
#RawVar


class Allele(ComplexModel):
    """
    Used in MutalyzerOutput data type.
    """
    __namespace__ = SOAP_NAMESPACE

    description = Mandatory.Unicode
    allele = Array(RawVar)
#Allele


class ExonInfo(ComplexModel):
    """
    Used in TranscriptInfo and MutalyzerOutput data types.
    """
    __namespace__ = SOAP_NAMESPACE

    cStart = Mandatory.Unicode
    gStart = Mandatory.Integer
    chromStart = Integer
    cStop = Mandatory.Unicode
    gStop = Mandatory.Integer
    chromStop = Integer
#ExonInfo


class MutalyzerOutput(ComplexModel):
    """
    Return type of SOAP method runMutalyzer.
    """
    __namespace__ = SOAP_NAMESPACE

    referenceId = Mandatory.Unicode
    sourceId = Mandatory.Unicode
    sourceAccession = Unicode
    sourceVersion = Unicode
    sourceGi = Unicode
    molecule = Mandatory.Unicode

    original = Unicode
    mutated = Unicode

    origMRNA = Unicode
    mutatedMRNA= Unicode

    origCDS = Unicode
    newCDS= Unicode

    origProtein = Unicode
    newProtein = Unicode
    altProtein = Unicode

    errors = Integer
    warnings = Integer
    summary = Unicode

    chromDescription = Unicode
    genomicDescription = Unicode
    transcriptDescriptions = Array(Unicode)
    proteinDescriptions = Array(Unicode)

    exons = Array(ExonInfo)

    rawVariants = Array(RawVariant)

    messages = Array(SoapMessage)
#MutalyzerOutput


class TranscriptNameInfo(ComplexModel):
    """
    Return type of SOAP method getGeneAndTranscript.
    """
    __namespace__ = SOAP_NAMESPACE

    transcriptName = Mandatory.Unicode
    productName = Mandatory.Unicode
#TranscriptNameInfo


class ProteinTranscript(ComplexModel):
    """
    Used in TranscriptInfo data type.
    """
    __namespace__ = SOAP_NAMESPACE

    name = Mandatory.Unicode
    id = Mandatory.Unicode
    product = Mandatory.Unicode
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

    name = Mandatory.Unicode
    id = Mandatory.Unicode
    product = Mandatory.Unicode

    cTransStart = Mandatory.Unicode
    gTransStart = Mandatory.Integer
    chromTransStart = Integer
    cTransEnd = Mandatory.Unicode
    gTransEnd = Mandatory.Integer
    chromTransEnd = Integer
    sortableTransEnd = Mandatory.Integer

    cCDSStart = Mandatory.Unicode
    gCDSStart = Mandatory.Integer
    chromCDSStart = Integer
    cCDSStop = Mandatory.Unicode
    gCDSStop = Mandatory.Integer
    chromCDSStop = Integer

    locusTag = Mandatory.Unicode
    linkMethod = Mandatory.Unicode

    exons = Array(ExonInfo)

    proteinTranscript = ProteinTranscript
#TranscriptInfo


class TranscriptMappingInfo(ComplexModel):
    """
    Used in return type of SOAP method getTranscriptsRange.
    """
    __namespace__ = SOAP_NAMESPACE

    name = Mandatory.Unicode
    version = Mandatory.Integer
    gene = Mandatory.Unicode
    orientation = Mandatory.Unicode

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

    version = Unicode
    versionParts = Array(Unicode)
    releaseDate = Unicode
    nomenclatureVersion = Unicode
    nomenclatureVersionParts = Array(Unicode)
    serverName = Unicode
    contactEmail = Unicode
    announcement = Unicode
    announcementUrl = Unicode
#InfoOutput


class CacheEntry(ComplexModel):
    """
    Used in getCache SOAP method.
    """
    __namespace__ = SOAP_NAMESPACE

    name = Mandatory.Unicode
    gi = Unicode
    hash = Mandatory.Unicode
    chromosomeName = Unicode
    chromosomeStart = Integer
    chromosomeStop = Integer
    chromosomeOrientation = Integer
    url = Unicode
    created = Mandatory.DateTime
    cached = Unicode
#CacheEntry
