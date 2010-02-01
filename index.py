#!/usr/bin/python

import os
import sys
from mod_python import apache
from cStringIO import StringIO
from simpletal import simpleTAL, simpleTALES

directory = os.path.dirname(__file__)
Mutalyzer = apache.import_module('Mutalyzer', path=[directory + "/src"])
VI = apache.import_module('Variant_info', path=[directory + "/src"])
os.chdir(directory)

Mut_ver = "2.0&alpha;"

def __run(func, *args) :
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    func(*args)
    reply = sys.stdout.getvalue()
    sys.stdout = old_stdout

    return reply
#__run

"""
def __htmlread(filename) :
    handle = open("./html/" + filename, "r")
    s = handle.read()
    handle.close

    return s
#__htmlread
"""

def __tal(filename, args) :
    context = simpleTALES.Context()

    for i in args :
        context.addGlobal(i, args[i])
    
    templateFile = open(filename, 'r')
    template = simpleTAL.compileHTMLTemplate(templateFile)
    templateFile.close()

    string = StringIO()
    template.expand(context, string)

    return string.getvalue()
#__tal

def index(req) :
    name = ""
    reply = ""
    if req.form :
        name = req.form['mutationName']
        reply = __run(Mutalyzer.main, name)
    #if

    args = {
        "version"    : Mut_ver, 
        "lastpost"   : name, 
        "mut_output" : reply
    }

    return __tal("html/check.html", args)
#index

def Variant_info(req) :
    LOVD_ver = req.form['LOVD_ver']
    build = req.form['build']
    acc = req.form['acc']
    var = ""
    if req.form.has_key('var') :
        var = req.form['var']

    result = __run(VI.main, LOVD_ver, build, acc, var)

    if LOVD_ver == "2.0-23" : # Obsoleted error messages, remove when possible.
        import re

        return re.sub("^Error \(.*\):", "Error:", result)
    #if
    return result
#Variant_info
