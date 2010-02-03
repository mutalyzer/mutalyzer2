#!/usr/bin/python

import sys
from cStringIO import StringIO
from simpletal import simpleTAL, simpleTALES

class Web() :
    def __init__(self) :
        self.version = "2.0&alpha;"
    #__init__

    def run(self, func, *args) :
        old_stdout = sys.stdout
        sys.stdout = StringIO()
        func(*args)
        reply = sys.stdout.getvalue()
        sys.stdout = old_stdout
        
        return reply 
    #run

    def tal(self, scheme, filename, args) :
        context = simpleTALES.Context()
    
        for i in args :
            context.addGlobal(i, args[i])
    
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

    def read(self, path, req) :
        handle = open(path + req.uri.split('/')[-1], "r")
        s = handle.read()
        handle.close
    
        return s
    #read
#Web
