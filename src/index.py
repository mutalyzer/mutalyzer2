#!/usr/bin/python

"""
    The HTML publisher.

    These functions appear as HTML pages on the web server.

    Public methods:
        index(req)        ; The mutation checker page.
        Variant_info(req) ; The g. to c. and vice versa interface for LOVD.
        download(req)     ; The download page.
"""

import Mutalyzer
import VarInfo
import pydoc
import webservice
import string
#import sys

from mod_python import apache

from Modules import Parser
from Modules import Mapper
from Modules import Web
from Modules import Config
from Modules import Output
from Modules import Db
from Modules import Scheduler
from Modules import Retriever
from Modules import File

def index(req) :
    W = Web.Web()
    return W.tal("HTML", "templates/index.html", {})
#index

def check(req) :
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
    C = Config.Config()
    O = Output.Output(__file__, C.Output)

    name = ""
    reply = ""
    if req.form :
        name = req.form['mutationName']
        O.addMessage(__file__, -1, "INFO", "Received variant %s" % name)
        RD = Mutalyzer.process(name, C, O)
        O.addMessage(__file__, -1, "INFO", "Finished processing variant %s" % \
                     name)
    #if
    errors, warnings, summary = O.Summary()

    altStart = O.getOutput("altstart")
    if altStart :
        altStart = altStart[0]

    genomicDescription = O.getOutput("genomicDescription")
    if genomicDescription :
        genomicDescription = genomicDescription[0]

    args = {
        "version"            : W.version,
        "lastpost"           : name,
        "messages"           : O.getMessages(),
        "summary"            : summary,
        "genomicDescription" : genomicDescription,
        "visualisation"      : O.getOutput("visualisation"),
        "descriptions"       : O.getOutput("descriptions"),
        "protDescriptions"   : O.getOutput("protDescriptions"),
        "oldProtein"         : O.getOutput("oldProteinFancy"),
        "altStart"           : altStart,
        "altProtein"         : O.getOutput("altProteinFancy"),
        "newProtein"         : O.getOutput("newProteinFancy"),
        "legends"            : O.getOutput("legends"),
        "mut_output"         : None
    }

    return W.tal("HTML", "templates/check.html", args)
#check

def syntaxCheck(req) :
    """
        Checks the syntax of a variant

        Arguments:
            req ; The request:
                  req.form['variant']  ; A description of the variant.

        Returns:
            string ; An HTML page containing the remark if the variant syntax
                     is OK or not
    """

    W = Web.Web()
    C = Config.Config() # Read the configuration file.
    O = Output.Output(__file__, C.Output)
    variant = req.form.get("variant", None)
    args = {
        "version"       : W.version,
        "variant"       : variant,
        "messages"      : [],
        "debug"         : ""
    }
    if variant:
        if variant.find(',') >= 0:
            O.addMessage(__file__, 2, "WCOMMASYNTAX",
                    "Comma's are not allowed in the syntax, autofixed")
            variant = variant.replace(',', '')
            args["variant"]=variant
        P = Parser.Nomenclatureparser(O)
        parsetree = P.parse(variant)
        if not parsetree :
            args["messages"].append("This variant does not have the right"
                " syntax. Please try again.")
            args["messages"].extend(O.getMessages())
        else :
            args["messages"].append("The syntax of this variant is OK!")
    #if
    ret = W.tal("HTML", "templates/parse.html", args)
    del W
    return ret
#checkingSyntax

def positionConverter(req):
    W = Web.Web()
    C = Config.Config()
    O = Output.Output(__file__, C.Output)

    if not req.form: req.form = {}
    build = req.form.get("build", "")
    variant = req.form.get("variant", "")

    attr = {
        "version" :         W.version,
        "avail_builds" :    C.Db.dbNames[::-1],
        "variant"  :        variant,
        "gName"     :       "",
        "cNames"    :       [],
        "messages" :        [],
        "errors" :          [],
        "debug" :           []
        }

    if build and variant:
        converter = Mapper.Converter(build, C, O)

        #Conver chr accNo to NC number
        variant = converter.correctChrVariant(variant)

        if not(":c." in variant or ":g." in variant):
            #Bad name
            P = Parser.Nomenclatureparser(O)
            parsetree = P.parse(variant)

        if ":c." in variant:
            # Do the c2chrom dance
            variant = converter.c2chrom(variant)
        if variant and ":g." in variant:
            # Do the g2c dance
            variants = converter.chrom2c(variant)
            if variants:
                attr["gName"] = variant
                output = ["%-10s:\t%s" % (key[:10], "\n\t\t".join(value))\
                        for key, value in variants.items()]
                attr["cNames"].extend(output)

        attr["errors"].extend(O.getMessages())
    return W.tal("HTML", "templates/converter.html", attr)
#positionConverter

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

    O = Output.Output(__file__, C.Output)
    D = Db.Cache(C.Db)
    R = Retriever.GenBankRetriever(C.Retriever, O, D)

    if req.method == 'POST' :
        if req.form["invoermethode"] == "file" :
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
        if req.form["invoermethode"] == "gene" :
            geneName = req.form["genesymbol"]
            organism = req.form["organism"]
            upStream = int(req.form["5utr"])
            downStream = int(req.form["3utr"])

            UD = R.retrievegene(geneName, organism, upStream, downStream)
            req.write(UD)
            return None
        #if
        if req.form["invoermethode"] == "chr" :
            accNo = req.form["chracc"]
            start = int(req.form["start"])
            stop = int(req.form["stop"])
            orientation = int(req.form["orientation"])
            UD = R.retrieveslice(accNo, start, stop, orientation)
            req.write(UD)
            return None
        #if
    #if

    W = Web.Web()
    args = {"version" : W.version}
    ret = W.tal("HTML", "templates/gbupload.html", args)
    del W
    return ret
#upload

def batch(req, batchType=None):
    """
        Batch function to add batch jobs to the Database
    """
    W = Web.Web()
    C = Config.Config()
    O = Output.Output(__file__, C.Output)

    attr = {"version"       : W.version,
            "messages"      : [],
            "errors"        : [],
            "debug"         : [],
            "batchTypes"    : ["NameChecker",
                               "SyntaxChecker",
                               "ConversionChecker"],
            "hideTypes"     : batchType and 'none' or '',
            "selected"      : "0",
            "batchType"     : batchType or "",
            "avail_builds"  : C.Db.dbNames[::-1]
         }

    # Use an empty dictionary if no form is filed
    if not(req.form): req.form = {}

    #get email and inFile
    email = req.form.get("batchEmail", None)
    inFile = req.form.get("batchFile", None)
    arg1 = req.form.get("arg1", "")

    #Make sure the correct page is displayed for an entrypoint
    batchType =  req.form.get("batchType", batchType or "NameChecker")
    if batchType in attr["batchTypes"]:
        attr["selected"] = str(attr["batchTypes"].index(batchType))

    if email and W.isEMail(email) and inFile:
        D = Db.Batch(C.Db)
        S = Scheduler.Scheduler(C.Scheduler, D)
        FileInstance = File.File(C.File, O)

        # Generate the fromhost URL from which the results can be fetched
        fromHost = "http://%s%s" % (
            req.hostname, req.uri.rsplit("/", 1)[0]+"/")

        job = FileInstance.parseBatchFile(inFile.file)
        if job is None:
            O.addMessage(__file__, 4, "PRSERR", "Could not parse input"
                " file, please check your file format.")
        else:
            #TODO: Add Binair Switches to toggle some events
            S.addJob("BINSWITHCES", email, job, fromHost, batchType, arg1)
            attr["messages"].append("Your file has been parsed and the job"
                " is scheduled, you will receive an email when the job is "
                "finished.")

        attr["errors"].extend(O.getMessages())

    return W.tal("HTML", "templates/batch.html", attr)
#batch

def batchNameChecker(req):
    return batch(req, "NameChecker")
#batchCheck

def batchConversionChecker(req):
    return batch(req, "ConversionChecker")
#batchConvert

def batchSyntaxChecker(req):
    return batch(req, "SyntaxChecker")
#batchCheckSyntaxch

def documentation(req) :
    """
        Generate documentation for the webservice.

        Arguments:
            req ; The request.

        Returns:
            string ; An HTML page.
    """

    htmldoc = pydoc.HTMLDoc()
    doc = "<html><body>%s</body></html>" % htmldoc.docmodule(webservice)
    return doc
#documentation

#TODO: taltest.html does not exist
def taltest(req) :
    W = Web.Web()
    C = Config.Config()
    variant = ""
    finalbuilds = []
    availBuilds = C.Db.dbNames  # available builds in config file
    for x in availBuilds :
        builds = {}
        builds["build"] = x
        finalbuilds.append(builds)
#    build = availBuilds[len(availBuilds)-1] # default to the highest build
    build = "hg18"
    args = {
        "version" : W.version,
        "build" : build,
        "avail_builds" : finalbuilds,
        "variant"  : variant
        }
    ret = W.tal("HTML", "templates/taltest.html", args)
    del W
    return ret
#taltest
