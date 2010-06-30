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
    ret = W.tal("HTML", "templates/index.html", {})

    del W
    return ret
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

    ret = W.tal("HTML", "templates/check.html", args)
    del W
    return ret
#check

def checkingSyntax(req) :
    """
        Checks the syntax of a variant
        
        Arguments:
            req ; The request:
                  req.form['build']    ; The human genome build.
                  req.form['variant']  ; A description of the variant.

        Returns:
            string ; An HTML page containing the remark if the variant syntax
                     is OK or not
    """
    
    W = Web.Web()
    C = Config.Config() # Read the configuration file.
    L = Output.Output(__file__, C.Output)
    variant = ""
    args = {
        "version"       : W.version,
        "variant"       : variant,
        "parse_output"  : ""
    }
    if req.form.has_key('variant') :
#    if req.form :
        variant  = string.replace(req.form['variant'], ",", "")
        args["variant"] = variant
        if variant == "" :
            args["parse_output"] = "You did not fill in the variant name"\
                " field.Please enter a variant using the format as shown"\
                " above."
            ret = W.tal("HTML", "templates/parse.html", args)
            del W
            return ret
        #if
        P = Parser.Nomenclatureparser(L)
        parsetree = P.parse(variant)
        del P
        if not parsetree :
            args["parse_output"] = "This variant does not have the right"\
                " syntax. Please try again."
            ret = W.tal("HTML", "templates/parse.html", args)
            del W
            return ret
        #if
        else :
            args["parse_output"] = "The syntax of this variant is OK!"
    #if
    ret = W.tal("HTML", "templates/parse.html", args)
    del W
    return ret
#checkingSyntax

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
    build = ""
    variant = ""
    reply = ""
    replyTranscripts = []
    final = []
    title = ""
    title2 = ""
    finalbuilds = []
    availBuilds = C.Db.dbNames  # available builds in config file
    for x in availBuilds :
        builds = {}
        builds["build"] = x
        finalbuilds.append(builds)
    build = availBuilds[len(availBuilds)-1] # default to the highest build
    args = {
        "version" : W.version,
        "build" : build,
        "avail_builds" : finalbuilds,
        "variant"  : variant
        }
    """
    availBuilds = C.Db.dbNames  # available builds in config file
    args = {
        "version"  : W.version,
        "build"    : build,
        "variant"  : variant
        }
    """
    # Check the build
#    if req.form.has_key('build') :
    if req.form :
        build = req.form['build']
        args["build"] = build
        buildStr = ', '.join(availBuilds)
        if req.form.has_key('variant') :
            args["variant"] = req.form['variant']
        # error if build not in available builds
        if build == "" :
            args["conv_output"] = "You did not provide the build. "\
                "Available builds are: " + buildStr
            ret = W.tal("HTML", "templates/convert.html", args)
            del W
            return ret
        #if
        if build not in availBuilds :
            args["conv_output"] = "This build is not available."\
                " Available builds are: " + buildStr
            ret = W.tal("HTML", "templates/convert.html", args)
            del W
            return ret
        #if
    #if
        
    D = Db.Mapping(build, C.Db)
    L = Output.Output(__file__, C.Output)
    
    if req.form.has_key('variant') :
        variant = req.form['variant']
        if variant == "" :
            args["conv_output"] = "You did not fill in the variant name field."\
                " Please enter a variant using the format as shown above."
            ret = W.tal("HTML", "templates/convert.html", args)
            del W
            return ret
        #if
        if "c." in variant :
            reply = Mapper.cTog(build, variant, C, D, L)
            if not reply :
                args["conv_output"] = "Nothing found. Either no transcripts"\
                " are known at that position or your input is incorrect. "\
                " Please check your input, e.g. did you provide the version"\
                " number or the correct change description?"
                ret = W.tal("HTML", "templates/convert.html", args)
                del W
                return ret
            #if
            mrnaAcc = variant.split(':')[0]
            accno = Mapper.mrnaSplit(mrnaAcc)[0]
            chrom = D.get_chromName(accno)
            chromacc = D.chromAcc(chrom)
            if not chromacc :
                chromacc = "NC_0000.0"
            reply = chromacc + ":" + reply
            # to check for other transcripts
            replyTranscripts = Mapper.gToc(build, reply, C, D, L)
            if not replyTranscripts :
                args["conv_output"] = "No transcripts found for: " + reply + " (the chromosome found was: " + chrom + " and the accession number was: " + accno
                ret = W.tal("HTML", "templates/convert.html", args)
                del W
                return ret
            #if
            for acc in replyTranscripts :
                genesDict = {}
                gene = D.get_GeneName(acc.split('.')[0])
                if not gene : # will this ever happen??
                    args["conv_output"] = "No gene found"
                    ret = W.tal("HTML", "templates/convert.html", args)
                    del W
                    return ret
                #if
                genesDict['gene'] = gene
                genesDict['mrnaAcc'] = acc
                final.append(genesDict)
            #for
            title2 = "Chromosomal position: "
            title = "Genes and transcripts found associated with this position"\
                    " (fully and partially overlapping, and including"\
                    " the provided one), using build %s:" % (build)
        #if
        if "g." in variant :
            #check if there is a minus (-) sign
            if "-" in variant :
                args["conv_output"] = "Invalid character (-) in your input, "\
                    "please correct."
                ret = W.tal("HTML", "templates/convert.html", args)
                del W
                return ret
            #if
                
            # check if build and NC_ accNo agree
            if "NC_" in variant :
                accNo = variant.split(':')[0] # the chromosome accession number
                # get the chromosome name
                chrom = D.chromName(accNo)
                if not chrom :
                    args["conv_output"] = "This chromosome accession number"\
                        " was not found in this build. You could either use"\
                        " the chromosome name (e.g. chr1) or look for the"\
                        " correct accession number at "
                    args["link_output"] = {
                                        "title" : "NCBI",
                                        "page" : "http://www.ncbi.nlm.nih.gov/"\
                                        "genome/guide/human/release_notes.html"
                                        }
                    ret = W.tal("HTML", "templates/convert.html", args)
                    del W
                    return ret
                #if
            #if
            if "chr" in variant :
                # Get the chromosome accession number (NC)
                mrnaAcc = D.chromAcc(variant.split(':')[0])
                variant = mrnaAcc + ":" + variant.split(':')[1]
            #if
            replyTranscripts = Mapper.gToc(build, variant, C, D, L)
            if not replyTranscripts :
                args["conv_output"] = "Nothing found. This can mean that"\
                        " either no transcripts were found at that position,"\
                        " or your input is incorrect (e.g. the position is"\
                        " not on the chromosome or the variant description"\
                        " is incomplete)."
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
            title = "Genes and transcripts found associated with this position"\
                 " (fully and partially overlapping), using build %s:" % (build)
        #if
    #if
    # now sort the dictionary genesDict
    final.sort()
    # use the TAL template
    args["conv_output"] = reply
    args["genes_output"] = final
    args["title"] = title
    args["title2"] = title2
    ret = W.tal("HTML", "templates/convert.html", args)
    del W
    return ret
#numberingConversion

def batchConversion(req) :
    """
        The c. to g. notation or vice versa for the batch web interface
        
        Arguments:
            req ; The request:
                  req.form['build']    ; The human genome build.
                  req.form['fileName']  ; A file name.

        Returns:
            
    """

    W = Web.Web()
    C = Config.Config() # Read the configuration file.
    build = ""
    fileName = ""
    eMail = ""
    reply = ""
    message = ""
    replyTranscripts = []
    final = []
    finalbuilds = []
    availBuilds = C.Db.dbNames  # available builds in config file
    for x in availBuilds :
        builds = {}
        builds["build"] = x
        finalbuilds.append(builds)
#    build = availBuilds[len(availBuilds)-1] # default to the highest build
    build = ""
    args = {
        "version" : W.version,
        "build" : build,
        "avail_builds" : finalbuilds,
        "fileName" : fileName,
        "eMail" : eMail,
        "message" : message
        }
    # Check the build
    if req.form.has_key('build') :
        build = req.form['build']
        args["build"] = build
        buildStr = ', '.join(availBuilds)
        # error if build not in available builds
        if build == "" :
            args["message"] = "You did not provide the build. "\
                                  "Available builds are: " + buildStr
            ret = W.tal("HTML", "templates/batch_convert.html", args)
            del W
            return ret
        #if
        if build not in availBuilds :
            args["message"] = "This build is not available."\
                " Available builds are: " + buildStr
            ret = W.tal("HTML", "templates/batch_convert.html", args)
            del W
            return ret
        #if
    #if
    D = Db.Mapping(build, C.Db)
    L = Output.Output(__file__, C.Output)
    
    if req.form.has_key('fileName') :
        fileName = req.form['fileName']
        if fileName == "" :
            args["message"] = "You did not provide a file name."\
                " Please enter a file name in the appropriate field."
            ret = W.tal("HTML", "templates/batch_convert.html", args)
            del W
            return ret
        #if
        
    #if
    if req.form.has_key('eMail') :
        eMail = req.form['eMail']
        if eMail == "" :
            args["message"] = "You did not provide an email"\
                " address. Please enter an email address."
            ret = W.tal("HTML", "templates/batch_convert.html", args)
            del W
            return ret
        #if
        fileName = req.form['fileName']
        
#        if fileName.filename and W.isEMail(eMail) :
#            C = Config.Config()
#            D = Db.Batch(C.Db)
#            S = Scheduler.Scheduler(C.Scheduler, D)
#            O = Output.Output(__file__, C.Output)
#            FileInstance = File.File(C.File, O)

#            job = FileInstance.parseBatchFile(fileUpload.file)
#            S.addJob("1231243", eMail, job, "http://%s%s" % (req.hostname, 
#                                                             req.uri))

#            del FileInstance, S, D, C
        inputFile = open(fileName, 'r')
        c_outputFile = open("batchconversion_to_cnum.txt", 'w')
        g_outputFile = open("batchconversion_to_gnum.txt", 'w')
        for line in inputFile :
            variant = line.strip()
            if "g." in variant :
                info = variant.split(':')
                chrom = info[0]
                if "chr" in chrom :
                    accNo = D.chromAcc(chrom)
                    var = accNo + ":" + info[1]
                replyTranscripts = Mapper.gToc(build, var, C, D, L)
    #            final.append(var)
                for acc in replyTranscripts :
#                    genesDict = {}
#                    gene = D.get_GeneName(acc.split('.')[0])
#                    genesDict['var'] = var
#                    genesDict['gene'] = gene
#                    genesDict['mrnaAcc'] = acc
#                    final.append(genesDict)
    #                outputFile.write(var + "\t")
    #                outputFile.write(gene + "\t")
                    c_outputFile.write(acc + "\n")
                #for
#                c_outputFile.close
            #if
            if "c." in variant :
                reply = Mapper.cTog(build, variant, C, D, L)
                if not reply :
                    var = ""
                else :
                    mrnaAcc = variant.split(':')[0]
                    accno = Mapper.mrnaSplit(mrnaAcc)[0]
                    chrom = D.get_chromName(accno)
                    chromacc = D.chromAcc(chrom)
                    if not chromacc :
                        chromacc = "NC_0000.0"
                    var = chromacc + ":" + reply
#                genesDict = {}
#                genesDict['var'] = variant
#                genesDict['mrnaAcc'] = var
#                final.append(genesDict)
#                outputFile.write(var + "\t")
#                outputFile.write(gene + "\t")
                g_outputFile.write(var + "\n")
#                g_outputFile.close
                
        #for
        inputFile.close
        c_outputFile.close
        g_outputFile.close
        args["message"] = "Conversion completed! You can find the resulting"\
        " file here:"
        args["link_output"] = {
                            "title" : "output",
                            "page" : "http://localhost/batchconversion_to_cnum.txt"
                            }
    #if
    args["batch_output"] = final
    ret = W.tal("HTML", "templates/batch_convert.html", args)
    del W
    return ret
    
#batchConversion
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
    R = Retriever.Retriever(C.Retriever, O, D)

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
