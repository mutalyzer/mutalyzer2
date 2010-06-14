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

from Modules import Mapper
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

from SOAPpy import WSDL

def numberingConversion(req) :
    """
        The c. to g. notation or vice versa for the web interface
        
        Arguments:
            req ; The request:
                  req.form['build']    ; The human genome build.
                  req.form['variant']  ; A description of the variant.

        Returns:
            string ; An HTML page containing the results of either cTog or gToc
    """
    W = Web.Web()
    C = Config.Config() # Read the configuration file.
    availBuilds = C.Db.dbNames  # available builds in config file
    build = availBuilds[1] # default to the highest build
    if req.form :
        build = req.form['build']
    # error if build not in available builds
    if build not in availBuilds :
        req.write("This build is not available.")
        return None
    D = Db.Db("local", build, C.Db)
    L = Output.Output(__file__, C.Output)
    
    variant = ""
    reply = ""
    replyTranscripts = []
    final = []
    title = ""
    title2 = ""
        
    if req.form.has_key('variant') :
        variant = req.form['variant']
        '''
        if ">" in variant :
            import re
            matchstr = re.compile(r'.+([ATGC]>[ATGC])$')
            if not matchstr.match(variant) :
                reply = "Input does not match the pattern [ATGC]>[ATGC]"
                args = {
                    "version"     : W.version, 
                    "build"       : build,
                    "variant"     : variant,
                    "conv_output" : reply
                }
            
                ret = W.tal("HTML", "templates/convert.html", args)
                del W
                return ret
        #if  
        '''
        if "c." in variant :
            reply = Mapper.cTog(build, variant, C, D, L)
            if not reply :
                reply = "Nothing found, please check your input"
            mrnaAcc = variant.split(':')[0]
            accno = Mapper.mrnaSplit(mrnaAcc)[0]
            chrom = D.get_chromName(accno)
            chromacc = D.chromAcc(chrom)
            reply = chromacc + ":" + reply
            # to check for other transcripts
            replyTranscripts = Mapper.gToc(build, reply, C, D, L)
            for acc in replyTranscripts :
                genesDict = {}
                gene = D.get_GeneName(acc.split('.')[0])
                if not gene :
                    reply = "No gene found"
#                    req.write("No gene found")
#                    return None
                genesDict['gene'] = gene
                genesDict['mrnaAcc'] = acc
                final.append(genesDict)
            #for
            title2 = "Chromosomal position: "
            title = "Genes and transcripts found associated with this variant"\
                    " (including the provided one):"
        #if
        if "g." in variant :
            # check if build and NC_ accNo agree
            if "NC_" in variant :
                accNo = variant.split(':')[0] # the chromosome accession number
                # get the chromosome name
                chrom = D.chromName(accNo)
                if not chrom :
                    reply = "This chromosome accession number was not found in this build"
                    args = {
                        "version"     : W.version, 
                        "build"       : build,
                        "variant"     : variant,
                        "conv_output" : reply
                    }
                
                    ret = W.tal("HTML", "templates/convert.html", args)
                    del W
                    return ret
            #if
            if "chr" in variant :
                # Get the chromosome accession number (NC)
                mrnaAcc = D.chromAcc(variant.split(':')[0])
                variant = mrnaAcc + ":" + variant.split(':')[1]
            #if
            replyTranscripts = Mapper.gToc(build, variant, C, D, L)
            if not replyTranscripts :
                reply = "Nothing found, please check your input"
                args = {
                    "version"     : W.version, 
                    "build"       : build,
                    "variant"     : variant,
                    "conv_output" : reply
                }
            
                ret = W.tal("HTML", "templates/convert.html", args)
                del W
                return ret
            #if not
            for acc in replyTranscripts :
                genesDict = {}
                gene = D.get_GeneName(acc.split('.')[0])
                genesDict['gene'] = gene
                genesDict['mrnaAcc'] = acc
                final.append(genesDict)
            #for
            title = "Genes and transcripts found associated with this variant:"
        #if
    #if
    # now sort the dictionary genesDict
    final.sort()
    # use the TAL template
    args = {
        "version"     : W.version, 
        "build"       : build,
        "variant"     : variant,
        "conv_output" : reply,
#        "trans_output": replyTranscripts,
        "genes_output": final,
        "title"       : title,
        "title2"      : title2
    }
    ret = W.tal("HTML", "templates/convert.html", args)
    del W
    return ret
#numberingConversion

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
