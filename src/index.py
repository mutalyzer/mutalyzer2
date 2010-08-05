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

from mod_python import apache, Session, util
from mod_python import Session, util

from Modules import Parser
from Modules import Mapper
from Modules import Web
from Modules import Config
from Modules import Output
from Modules import Db
from Modules import Scheduler
from Modules import Retriever
from Modules import File

def snp(req) :
    C = Config.Config()
    O = Output.Output(__file__, C.Output)
    W = Web.Web()

    rsId = None
    if req.form :
        rsId = req.form.get('rsId', None)
    if rsId :
        O.addMessage(__file__, -1, "INFO", "Received rs%s" % rsId)
        R = Retriever.Retriever(C.Retriever, O, None)
        R.snpConvert(rsId)
        O.addMessage(__file__, -1, "INFO", "Finished processing rs%s" % rsId)

    args = {
        "snp"      : O.getOutput("snp"),
        "messages" : O.getMessages(),
        "summary"  : O.Summary()[2],
        "lastpost" : rsId
    }

    return W.tal("HTML", "templates/snp.html", args)
#snp

class InputException(Exception):
    pass

def index(req) :
    W = Web.Web()
    return W.tal("HTML", "templates/index.html", {})
#index

def nameGenerator(req):
    W = Web.Web()
    return W.tal("HTML", "templates/generator.html", {})
#generator

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
        name = req.form.get('mutationName', None)
    else:
        session = Session.Session(req)
        session.load()
        name = session.get("mut", None)
        #Remove from session
        session["mut"] = None
        session.save()
    if name:
        O.addMessage(__file__, -1, "INFO", "Received variant %s" % name)
        RD = Mutalyzer.process(name, C, O)
        O.addMessage(__file__, -1, "INFO", "Finished processing variant %s" % \
                     name)
    #if
    errors, warnings, summary = O.Summary()
    recordType = O.getIndexedOutput("recordType",0)
    reference = O.getIndexedOutput("reference", 0)
    if recordType == "LRG" :
        reference += ".xml"
    else :
        reference += ".gb"

    args = {
        "lastpost"           : name,
        "messages"           : O.getMessages(),
        "summary"            : summary,
        "parseError"         : O.getOutput("parseError"),
        "errors"             : errors,
        "genomicDescription" : W.urlEncode([O.getIndexedOutput(
                                   "genomicDescription", 0)])[0],
        "chromDescription"   : O.getIndexedOutput("genomicChromDescription", 0),
        "visualisation"      : O.getOutput("visualisation"),
        "descriptions"       : W.urlEncode(O.getOutput("descriptions")),
        "protDescriptions"   : O.getOutput("protDescriptions"),
        "oldProtein"         : O.getOutput("oldProteinFancy"),
        "altStart"           : O.getIndexedOutput("altStart", 0),
        "altProtein"         : O.getOutput("altProteinFancy"),
        "newProtein"         : O.getOutput("newProteinFancy"),
        "exonInfo"           : O.getOutput("exonInfo"),
        "cdsStart_g"         : O.getIndexedOutput("cdsStart_g", 0),
        "cdsStart_c"         : O.getIndexedOutput("cdsStart_c", 0),
        "cdsStop_g"          : O.getIndexedOutput("cdsStop_g", 0),
        "cdsStop_c"          : O.getIndexedOutput("cdsStop_c", 0),
        "restrictionSites"   : O.getOutput("restrictionSites"),
        "legends"            : O.getOutput("legends"),
        "reference"          : reference
    }

    if req.method == 'GET' and req.form :
        args["interactive"] = False
        ret = W.tal_old("HTML", "templates/check.html", args)
    else :
        args["interactive"] = True
        ret = W.tal("HTML", "templates/check.html", args)
    del W
    return ret
#check

def checkForward(req) :
    session = Session.Session(req)
    session['mut'] = req.form.get("mutationName", None)
    session.save()
    util.redirect(req, "check", permanent=False)
#checkForward


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
    if variant:
        if variant.find(',') >= 0:
            O.addMessage(__file__, 2, "WCOMMASYNTAX",
                    "Comma's are not allowed in the syntax, autofixed")
            variant = variant.replace(',', '')
            args["variant"]=variant
        P = Parser.Nomenclatureparser(O)
        parsetree = P.parse(variant)
        #if not parsetree :
        #    args["messages"].append("This variant does not have the right"
        #    args["messages"].extend(O.getMessages())
        #else :
        #    args["messages"].append("The syntax of this variant is OK!")
    #if
    args = {
        "variant"       : variant,
        "messages"      : O.getMessages(),
        "parseError"    : O.getOutput("parseError"),
        "debug"         : ""
    }
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

    avail_builds = C.Db.dbNames[::-1]

    if build :
        avail_builds.remove(build)
        avail_builds.insert(0, build)
    #if

    attr = {
        "avail_builds" : avail_builds,
        "variant"      : variant,
        "gName"        : "",
        "cNames"       : [],
        "messages"     : [],
        "errors"       : [],
        "debug"        : []
        }

    if build and variant:
        converter = Mapper.Converter(build, C, O)

        #Conver chr accNo to NC number
        variant = converter.correctChrVariant(variant)

        if variant :
            if not(":c." in variant or ":g." in variant):
                #Bad name
                P = Parser.Nomenclatureparser(O)
                parsetree = P.parse(variant)
            #if

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
                #if
            #if
        #if

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
    var = req.form.get("var", "")

    result = W.run(VarInfo.main, LOVD_ver, build, acc, var)

    if LOVD_ver == "2.0-23" : # Obsoleted error messages, remove when possible.
        import re

        return re.sub("^Error \(.*\):", "Error:", result)
    #if
    return result
#Variant_info

def webservices(req) :
    """
        The download page.

        Arguments:
            req ; The request.

        Returns:
            string ; An HTML page.
    """

    W = Web.Web()

    ret = W.tal("HTML", "templates/webservices.html", {})
    del W
    return ret
#download

def __checkInt(inpv, refname):
    #remove , . and -
    inpv = inpv.replace(',','').replace('.','').replace('-','')
    try:
        return int(inpv)
    except ValueError, e:
        raise InputException("Expected an integer in field: %s" % refname)

def upload(req) :
    """
    """

    C = Config.Config()
    maxUploadSize = C.Retriever.maxDldSize

    O = Output.Output(__file__, C.Output)
    D = Db.Cache(C.Db)
    R = Retriever.GenBankRetriever(C.Retriever, O, D)

    UD, errors = "", []

    if req.method == 'POST' :
        try:
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
                UD = R.uploadrecord(req.form["bestandsveld"].file.read())
            #if
            elif req.form["invoermethode"] == "url" :
                UD = R.downloadrecord(req.form["urlveld"])
            #if
            elif req.form["invoermethode"] == "gene" :
                geneName = req.form["genesymbol"]
                organism = req.form["organism"]
                upStream = __checkInt((req.form["5utr"]),
                        "5' flanking nucleotides")
                downStream = __checkInt((req.form["3utr"]),
                        "3' flanking nucleotides")
                UD = R.retrievegene(geneName, organism, upStream, downStream)
            #if
            elif req.form["invoermethode"] == "chr" :
                accNo = req.form["chracc"]
                start = __checkInt((req.form["start"]),
                        "Start position")
                stop = __checkInt((req.form["stop"]),
                        "Stop position")
                orientation = int(req.form["orientation"])
                UD = R.retrieveslice(accNo, start, stop, orientation)
            #if
            else:
                #unknown "invoermethode"
                raise InputException("Wrong method selected")
        except InputException, e:
            #DUMB USERS
            errors.append(e)
        finally:
            if not UD:
                #Something went wrong
                errors += ["The request could not be completed"]
                errors.extend(O.getMessages())
    #if

    W = Web.Web()
    args = {
        "UD"      : UD,
        "maxSize" : float(maxUploadSize) / 1048576,
        "errors"  : errors
    }
    ret = W.tal("HTML", "templates/gbupload.html", args)
    del W
    return ret
#upload

def progress(req):
    """
        Progress page for batch runs
    """
    W = Web.Web()
    C = Config.Config()
    O = Output.Output(__file__, C.Output)

    attr = {"percentage"    : 0}

    try:
        jobID = int(req.form["jobID"])
        total = int(req.form["totalJobs"])
    except Exception, e:
        return
    D = Db.Batch(C.Db)
    left = D.entriesLeftForJob(jobID)
    percentage = int(100 - (100 * left / float(total)))
    if req.form.get("ajax", None):
        if percentage == 100:
            #download url, check if file still exists
            ret = "OK"
        else:
            ret = percentage
        return ret
    else:
        #Return progress html page
        return W.tal("HTML", "templates/progress.html", attr)


def batch(req, batchType=None):
    """
        Batch function to add batch jobs to the Database
    """
    W = Web.Web()
    C = Config.Config()
    O = Output.Output(__file__, C.Output)

    attr = {"messages"      : [],
            "errors"        : [],
            "debug"         : [],
            "batchTypes"    : ["NameChecker",
                               "SyntaxChecker",
                               "PositionConverter"],
            "hideTypes"     : batchType and 'none' or '',
            "selected"      : "0",
            "batchType"     : batchType or "",
            "avail_builds"  : C.Db.dbNames[::-1],
            "jobID"         : None,
            "totalJobs"     : None
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
            attr["jobID"] =\
                    S.addJob("BINSWITHCES", email, job, fromHost, batchType, arg1)
            attr["totalJobs"] = len(job) or 1
            attr["messages"].append("Your file has been parsed and the job"
                " is scheduled, you will receive an email when the job is "
                "finished.")

        attr["errors"].extend(O.getMessages())

    return W.tal("HTML", "templates/batch.html", attr)
#batch

def disclaimer(req) :
    W = Web.Web()
    return W.tal("HTML", "templates/disclaimer.html", [])
#disclaimer

def batchNameChecker(req):
    return batch(req, "NameChecker")
#batchCheck

def batchPositionConverter(req):
    return batch(req, "PositionConverter")
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
        "build" : build,
        "avail_builds" : finalbuilds,
        "variant"  : variant
        }
    ret = W.tal("HTML", "templates/taltest.html", args)
    del W
    return ret
#taltest
