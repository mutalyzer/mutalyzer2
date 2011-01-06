#!/usr/bin/python

"""
Collection of Serilizable Objects used by the webservice

@requires: soaplib.serializers.primitive.String
@requires: soaplib.serializers.primitive.Integer
@requires: soaplib.serializers.primitive.Array
@requires: soaplib.serializers.clazz.ClassSerializer

@todo: documentation
"""
from soaplib.model.primitive import String, Integer
from soaplib.model.clazz import ClassModel


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


class List(object):
    """
    Unbounded (in number of occurrences) versions of models.
    """
    String = String(nillable=False, max_occurs="unbounded")
    Integer = Integer(nillable=False, max_occurs="unbounded")


# Todo: use Mandatory.* models in the classmodels below?


class SoapMessage(ClassModel):
    """
    Send info message over the soapline

    Attributes:
        - errorcode   ; The error code affiliated with the error message
        - message     ; The error message
    """

    __namespace__ = 'soapmessage'

    errorcode = Mandatory.String
    message = Mandatory.String
#SoapMessage

class Mapping(ClassModel) :
    """
    Extended ClassSerializer object with mixed types of attributes

    Attributes:
        - startmain    ; Define the type of startmain.
        - startoffset  ; Define the type of startoffset.
        - endmain      ; Define the type of endmain value.
        - endoffset    ; Define the type of endoffset value.
        - start_g      ; Define the type of start_g value.
        - end_g        ; Define the type of end_g value.
        - mutationType ; Define the type of mutation type
    """

    __namespace__ = 'mapping'

    startmain = Integer
    startoffset = Integer
    endmain = Integer
    endoffset = Integer
    start_g = Integer
    end_g = Integer
    mutationType = String
    errorcode = Integer
    messages = SoapMessage(max_occurs="unbounded")
#Mapping

class Transcript(ClassModel) :
    """
    Extended ClassSerializer object with mixed types of attributes

    Attributes:
        - trans_start ; Define the type of trans_start
        - trans_stop  ; Define the type of trans_stop
        - CDS_stop    ; Define the type of CDS_stop
    """

    __namespace__ = 'transcript'

    trans_start = Integer
    trans_stop = Integer
    CDS_stop = Integer
#Transcript

class MutalyzerOutput(ClassModel) :
    """
        Extended ClassSerializer object with mixed types of attributes

        Attributes:
    """

    __namespace__ = 'mutalyzeroutput'

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
#MutalyzerOutput
