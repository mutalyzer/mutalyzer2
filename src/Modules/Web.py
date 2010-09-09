#!/usr/bin/python

"""
    Module that provides general functions used by the web interfaces.

    Public classes:
        Web ; General functions used by the web interfaces.
"""

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
            version ; This is the version that is displayed on the web pages,
                      WSDL files, etc.

        Special methods:
            __init__() ; Initialise the class.

        Public methods:
            run(func, *args)            ; Run func(*args) and return stdout.
            tal(scheme, filename, args) ; Compile a TAL template to HTML or XML.
            read(path, req)             ; Read a file and return the content.
    """

    def __init__(self) :
        """
            Initialise the class.

            Public variables (altered):
                version ; Here the displayed version is defined.
        """

        self.version = "2.0&nbsp;&beta;-2"
        self.nomenclatureVersion = "2.0"
        self.releaseDate = "9 Sep 2010"

        C = Config.Config()
        self.email = C.Retriever.email
    #__init__

    def run(self, func, *args) :
        """
            Run any function and return standard output as a string.

            Arguments:
                func  ; The function that has to be called.
                *args ; The arguments of func.

            Returns:
                string ; Everything that func(*args) writes to standard output.
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

            Arguments:
                scheme   ; Either "HTML" or "XML", output will be in this
                           format.
                filename ; The filename of the template.
                args     ; A dictionary with variables (whose name correspond
                           to the ones in the template) and their values.

            Returns:
                string ; An HTML or XML file.
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

            Arguments:
                scheme   ; Either "HTML" or "XML", output will be in this
                           format.
                filename ; The filename of the template.
                args     ; A dictionary with variables (whose name correspond
                           to the ones in the template) and their values.

            Returns:
                string ; An HTML or XML file.
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

            Arguments:
                path ; Path to the file.
                req  ; HTTP request (used to extract the filename).

            Returns:
                string ; The content of the file.
        """

        handle = open(path + req.uri.split('/', 2)[2], "r")
        s = handle.read()
        handle.close

        return s
    #read

    def isEMail(self, eMail) :
        #TODO documentation
        """
        """

        if re.match("^[a-zA-Z0-9._%-]+@[a-zA-Z0-9._%-]+.[a-zA-Z]{2,6}$",
                    eMail) :
            return True
        return False
    #isEmail

    def urlEncode(self, descriptions) :
        #TODO documentation
        """
        """

        newDescr = []
        for i in descriptions :
            newDescr.append([i, urllib.quote(i)])
        return newDescr
    #urlEncode
#Web
