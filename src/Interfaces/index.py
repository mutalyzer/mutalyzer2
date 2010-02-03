#!/usr/bin/python

import sys

import Mutalyzer
import Variant_info as VI
from Modules import Web

def index(req) :
    W = Web.Web()

    name = ""
    reply = ""
    if req.form :
        name = req.form['mutationName']
        reply = W.run(Mutalyzer.main, name)
    #if

    args = {
        "version"    : W.version, 
        "lastpost"   : name, 
        "mut_output" : reply
    }

    return W.tal("HTML", "templates/check.html", args)
#index

def Variant_info(req) :
    W = Web.Web()

    LOVD_ver = req.form['LOVD_ver']
    build = req.form['build']
    acc = req.form['acc']
    var = ""
    if req.form.has_key('var') :
        var = req.form['var']

    result = W.run(VI.main, LOVD_ver, build, acc, var)

    if LOVD_ver == "2.0-23" : # Obsoleted error messages, remove when possible.
        import re

        return re.sub("^Error \(.*\):", "Error:", result)
    #if
    return result
#Variant_info
