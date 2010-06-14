#!/usr/bin/python

"""
    The HTML publisher.

    These functions appear as HTML pages on the web server.
    
    Public methods:
        index(req)        ; The mutation checker page.
        Variant_info(req) ; The g. to c. and vice versa interface for LOVD.
        download(req)     ; The download page.
"""

#import sys

import Mutalyzer
import VarInfo
import pydoc
import webservice

from mod_python import apache

from Modules import Web
from Modules import Config
from Modules import Output
from Modules import Db
from Modules import Scheduler
from Modules import File

def index(req) :
    """
        The mutation checker page.

        If the incoming request has a form, run Mutalyzer. The output of 
        Mutalyzer is used together with a version and the last posted value
        to make an HTML page from a TAL template.

        Arguments:
            req ; The request:
                  req.form['mutationName'] ; A description of a variant.

        Returns:
            string ; An HTML page containing the results of Mutalyzer.
    """

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

    ret = W.tal("HTML", "templates/check.html", args)
    del W
    return ret
#index

def Variant_info(req) :
    """
        The g. to c. and vice versa interface for LOVD.

        Arguments:
            req ; The request:
                  req.form['LOVD_ver'] ; The version of the calling LOVD.
                  req.form['build']    ; The human genome build (hg19 assumed).
                  req.form['acc']      ; The accession number (NM number).
                  req.form['var']      ; A description of the variant.

        Returns:
            string ; An HTML page containing the results of Variant_info.
    """

    W = Web.Web()

    LOVD_ver = req.form['LOVD_ver']
    build = req.form['build']
    acc = req.form['acc']
    var = ""
    if req.form.has_key('var') :
        var = req.form['var']

    result = W.run(VarInfo.main, LOVD_ver, build, acc, var)

    if LOVD_ver == "2.0-23" : # Obsoleted error messages, remove when possible.
        import re

        return re.sub("^Error \(.*\):", "Error:", result)
    #if
    return result
#Variant_info

def download(req) :
    """
        The download page.

        Arguments:
            req ; The request.

        Returns:
            string ; An HTML page.
    """

    W = Web.Web()

    args = {"version" : W.version}
    ret = W.tal("HTML", "templates/download.html", args)
    del W
    return ret
#download

def upload(req) :
    """
    """

    C = Config.Config()
    maxUploadSize = C.Retriever.maxDldSize
    del C

    if req.method == 'POST' :
        length = req.headers_in.get('Content-Length')
        if not length :
            req.status = apache.HTTP_LENGTH_REQUIRED
            req.write("Content length required.")
            return None
        #if
        if int(length) > maxUploadSize :
            req.status = apache.HTTP_REQUEST_ENTITY_TOO_LARGE
            req.write("Upload limit exceeded.")
            return None
        #if
    #if

    W = Web.Web()
    args = {"version" : W.version}
    ret = W.tal("HTML", "templates/gbupload.html", args)
    del W
    return ret
#upload

def batch(req) :
    """
    """

    W = Web.Web()
    eMail = ""
    if req.form :
        eMail = req.form['eMail']
        fileUpload = req.form['file']
        
        if fileUpload.filename and W.isEMail(eMail) :
            C = Config.Config()
            D = Db.Batch(C.Db)
            S = Scheduler.Scheduler(C.Scheduler, D)
            O = Output.Output(__file__, C.Output)
            FileInstance = File.File(C.File, O)

            job = FileInstance.parseBatchFile(fileUpload.file)
            S.addJob("1231243", eMail, job, "http://%s%s" % (req.hostname, 
                                                             req.uri))

            del FileInstance, S, D, C
        #if            
    #if

    args = {
        "version" : W.version,
        "lastEMail" : eMail
    }

    ret = W.tal("HTML", "templates/batch.html", args)
    del W
    return ret
#batch

def documentation(req) :
    """
        Generate documentation for the webservice.

        Arguments:
            req ; The request.

        Returns:
            string ; An HTML page.
    """

    htmldoc = pydoc.HTMLDoc()
    doc = "<html><body>"
    doc += htmldoc.docmodule(webservice)
    del htmldoc
    doc += "</body></html>"

    return doc
#documentation

def menu(req) :
    W = Web.Web()

    return W.tal2("", "")
#menu
