#!/usr/bin/python

import os
import sys
from mod_python import apache
from cStringIO import StringIO

directory = os.path.dirname(__file__)
Main = apache.import_module('Main', path=[directory + "/src"])
VI = apache.import_module('Variant_info', path=[directory + "/src"])
os.chdir(directory)

s = """
<html>
    <head>
        <title>Mutalyzer 2.0</title>
    </head>
    <body>
        <center>
            <big>Mutalyzer 2.0 nomenclature check.</big>
        </center><br>
        Please insert the mutation name using the format:<br>
        <small>
          &lt;Accession Number&gt;.&lt;version 
          number&gt;(&lt;Gene symbol&gt;):&lt;sequence 
          type&gt;.&lt;mutation&gt;
        </small><br>
        <br>
        Example: AB026906.1:c.274G&gt;T<br>
        <br>
        <form action = "" method = "post">
            <input 
                type = "text" 
                name = "mutationName" 
                value = "%s"
                size = "100%%"
            >
        </form>
        <pre>
%s
        </pre>
    </body>
</html>
"""

def index(req) :
    name = ""
    reply = ""
    if req.form :
        name = req.form['mutationName']

        #fsock = open('panic.log', 'a')
        #sys.stderr = fsock

        old_stdout = sys.stdout
        sys.stdout = StringIO()
        Main.rrr(name)
        reply = sys.stdout.getvalue() #.replace('\n', "<br>")
        sys.stdout = old_stdout
    #if

    return s % (name, reply)
#index

def Variant_info(req) :
    LOVD_ver = req.form['LOVD_ver']
    build = req.form['build']
    acc = req.form['acc']
    var = req.form['var']

    old_stdout = sys.stdout
    sys.stdout = StringIO()
    VI.main(LOVD_ver, build, acc, var)
    reply = sys.stdout.getvalue() #.replace('\n', "<br>")
    sys.stdout = old_stdout

    return reply
#Variant_info
