#!/usr/bin/python

"""
Module that provides general functions used by the web interfaces.

@requires: sys
@requires: re
@requires: urllib 
@requires: cStringIO.StringIO
@requires: simpletal.simpleTALES
@requires: simpletal.simpleTAL
@requires: Config
"""
# Public classes:
#     - Web ; General functions used by the web interfaces.


import sys                        # sys.stdout
import re                         # match
import urllib                     # quote
from cStringIO import StringIO    # StringIO() getvalue()
from simpletal import simpleTALES # context(), addGlobal()
from simpletal import simpleTAL   # compileHTMLTemplate,
                                  # compileXMLTemplate,
import Config


class Web() :
    """
    General functions used by the web interfaces.

    Public variables:
        - version ; This is the version that is displayed on the web pages,
                    WSDL files, etc.

    Special methods:
        - __init__() ; Initialise the class.

    Public methods:
        - run(func, *args)            ; Run func(*args) and return stdout.
        - tal(scheme, filename, args) ; Compile a TAL template to HTML or XML.
        - read(path, req)             ; Read a file and return the content.
    """

    def __init__(self) :
        """
        Initialise the class.

        Public variables (altered):
            - version ; Here the displayed version is defined.
        """

        self.version = "2.0&nbsp;&beta;-6"
        self.nomenclatureVersion = "2.0"
        self.releaseDate = "12 Jan 2011"

        C = Config.Config()
        self.email = C.Retriever.email
    #__init__

    def run(self, func, *args) :
        """
        Run any function and return standard output as a string.

        @arg func: The function that has to be called
        @type func: function
        @arg args: arguments for the function to call
        @type args: list

        @return: Everything that func(*args) writes to standard output
        @rtype: string
        """

        old_stdout = sys.stdout
        sys.stdout = StringIO()
        func(*args)
        reply = sys.stdout.getvalue()
        sys.stdout = old_stdout

        return reply
    #run

    def tal_old(self, scheme, filename, args) :
        #TODO merge this function with 'tal' (below).
        """
        Compile a TAL template to HTML or XML.
        
        @todo: merge this function with 'tal'

        @arg scheme:    Either "HTML" or "XML", output will be in this format
        @type scheme:   string
        @arg filename:  The filename of the template
        @type filename: string
        @arg args:      A dictionary with variables (whose name correspond to
                        the ones in the template) and their values
        @type args:     dictionary

        @return: An HTML or XML file
        @rtype:  string
        """

        from simpletal import simpleTALES # context(), addGlobal()
        from simpletal import simpleTAL   # compileHTMLTemplate,
                                          # compileXMLTemplate,

        context = simpleTALES.Context()

        for i in args :
            context.addGlobal(i, args[i])

        #templateFile = open("templates/menu.html", 'r')
        #macros = simpleTAL.compileHTMLTemplate(templateFile)
        #templateFile.close()
        #context.addGlobal("sitemacros", macros)

        templateFile = open(filename, 'r')

        if scheme == "HTML" :
            template = simpleTAL.compileHTMLTemplate(templateFile)
        else :
            template = simpleTAL.compileXMLTemplate(templateFile)

        templateFile.close()

        string = StringIO()
        template.expand(context, string)

        return string.getvalue()
    #tal

    def tal(self, scheme, filename, args) :
        """
        Compile a TAL template to HTML or XML.

        @arg scheme:    Either "HTML" or "XML", output will be in this format
        @type scheme:   string
        @arg filename:  The filename of the template
        @type filename: string
        @arg args:      A dictionary with variables (whose name correspond to
                        the ones in the template) and their values
        @type args:     dictionary

        @return: An HTML or XML file
        @rtype: string
        """

        context = simpleTALES.Context()

        context.addGlobal("version", self.version)
        context.addGlobal("nomenclatureVersion", self.nomenclatureVersion)
        context.addGlobal("releaseDate", self.releaseDate)
        context.addGlobal("contactEmail", self.email)
        for i in args :
            context.addGlobal(i, args[i])

        if scheme == "HTML" :
            templateFile = open(filename, 'r')
            macros = simpleTAL.compileHTMLTemplate(templateFile)
            templateFile.close()
            context.addGlobal("sitemacros", macros)

            templateFile = open("templates/menu.html", 'r')
            template = simpleTAL.compileHTMLTemplate(templateFile)
            templateFile.close()
        #if
        else :
            templateFile = open(filename, 'r')
            template = simpleTAL.compileHTMLTemplate(templateFile)
            templateFile.close()
        #else

        string = StringIO()
        template.expand(context, string)

        return string.getvalue()
    #tal

    def read(self, path, req) :
        """
        Read a file and return its content.

        @arg path: Path to the file
        @type path: string
        @arg req: HTTP request (used to extract the filename)
        @type req: string

        @return: The content of the file
        @rtype: string
        """

        handle = open(path + req.uri.split('/', 2)[2], "r")
        s = handle.read()
        handle.close

        return s
    #read

    def isEMail(self, eMail) :
        #TODO documentation
        """
        Check if argument is a valid email address.
        
        @arg eMail: email address to check
        @type eMail: string
        
        @return: True or False
        @rtype: boolean
        """

        if re.match("^[a-zA-Z0-9._%-]+@[a-zA-Z0-9._%-]+.[a-zA-Z]{2,6}$",
                    eMail) :
            return True
        return False
    #isEmail

    def urlEncode(self, descriptions) :
        #TODO documentation
        """
        @arg descriptions:
        @type descriptions: list
        
        @return: urlEncode descriptions???????????????
        @rtype: list
        """

        newDescr = []
        for i in descriptions :
            newDescr.append([i, urllib.quote(i)])
        return newDescr
    #urlEncode
#Web
