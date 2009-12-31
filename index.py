#!/usr/bin/python

import os
import sys
from mod_python import apache
from cStringIO import StringIO

directory = os.path.dirname(__file__)
Main = apache.import_module('Main', path=[directory + "/src"])
VI = apache.import_module('Variant_info', path=[directory + "/src"])
os.chdir(directory)

def __run(func, *args) :
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    func(*args)
    reply = sys.stdout.getvalue()
    sys.stdout = old_stdout

    return reply
#__run

def __htmlread(filename) :
    handle = open("./html/" + filename, "r")
    s = handle.read()
    handle.close

    return s
#__htmlread

def index(req) :
    name = ""
    reply = ""
    if req.form :
        name = req.form['mutationName']
        reply = __run(Main.main, name)
    #if

    return __htmlread("check.html") % (name, reply)
#index

def Variant_info(req) :
    LOVD_ver = req.form['LOVD_ver']
    build = req.form['build']
    acc = req.form['acc']
    var = req.form['var']

    return __run(VI.main, LOVD_ver, build, acc, var)
#Variant_info
