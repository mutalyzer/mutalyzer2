#!/usr/bin/python

"""
Collection of Serilizable Objects used by the webservice

@requires: soaplib.serializers.primitive.String
@requires: soaplib.serializers.primitive.Integer
@requires: soaplib.serializers.primitive.Array
@requires: soaplib.serializers.clazz.ClassSerializer

@todo: documentation
"""
from soaplib.serializers.primitive import String, Integer, Array
from soaplib.serializers.clazz import ClassSerializer

class SoapMessage(ClassSerializer):
    """
    Send info message over the soapline

    Attributes:
        - errorcode   ; The error code affiliated with the error message
        - message     ; The error message
    """

    class types():
        errorcode = String
        message = String
    #types
#SoapMessage

class Mapping(ClassSerializer) :
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

    class types() :
        """
        Types are defined here for the soaplib module.
        """

        startmain = Integer
        startoffset = Integer
        endmain = Integer
        endoffset = Integer
        start_g = Integer
        end_g = Integer
        mutationType = String
        errorcode = Integer
        messages = Array(SoapMessage)
    #types
#Mapping

class Transcript(ClassSerializer) :
    """
    Extended ClassSerializer object with mixed types of attributes

    Attributes:
        - trans_start ; Define the type of trans_start
        - trans_stop  ; Define the type of trans_stop
        - CDS_stop    ; Define the type of CDS_stop
    """

    class types() :
        """
        """

        trans_start = Integer
        trans_stop = Integer
        CDS_stop = Integer
    #types
#Transcript

class MutalyzerOutput(ClassSerializer) :
    """
        Extended ClassSerializer object with mixed types of attributes

        Attributes:
    """

    class types() :
        """
        """

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
    #types
#MutalyzerOutput
