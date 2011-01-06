#!/usr/bin/python

"""
The HTML publisher.

These functions appear as HTML pages on the web server.

Public methods:
    - index(req)        ; The mutation checker page.
    - Variant_info(req) ; The I{g.} to I{c.} and vice versa interface for LOVD.
    - download(req)     ; The download page.

@requires: Mutalyzer
@requires: VarInfo
@requires: pydoc
@requires: webservice
@requires: string

@requires: mod_python import apache
@requires: mod_python import Session
@requires: mod_python import util

@requires: Modules.Parser
@requires: Modules.Mapper
@requires: Modules.Web
@requires: Modules.Config
@requires: Modules.Output
@requires: Modules.Db
@requires: Modules.Scheduler
@requires: Modules.Retriever
@requires: Modules.File
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

class InputException(Exception):
    """
    @todo: documentation
    """
    pass

def snp(req) :
    """
    @todo: documentation
    
    @arg req: the HTTP request
    @type req: object
    @return: compiled TAL template
    @rtype: object
    """
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
    #if

    args = {
        "snp"      : O.getOutput("snp"),
        "messages" : O.getMessages(),
        "summary"  : O.Summary()[2],
        "lastpost" : rsId
    }

    return W.tal("HTML", "templates/snp.html", args)
#snp

def getGS(req):
    """
    LOVD bypass to get the correct GeneSymbol incl Transcript variant.

    Used by LOVD to get the correct transcript variant out of a genomic
    record. LOVD uses a genomic reference (NC_?) in combination with a gene
    symbol to pass variant info to mutalyzer. Mutalyzer 1.0 was only using
    the first transcript. LOVD supplies the NM of the transcript needed but
    this was ignored. This helper allows LOVD to get the requested
    transcript variant from a genomic reference.

    @arg req: The request:
        - req.form['mutationName'] ; the mutationname without gene symbol
        - re.form['variantRecord'] ; the NM reference of the variant
        - re.form['forward']       ; if set this forwards the request to the name
                                     checker
    @type req:
    
    @return:
        - string ; The GeneSymbol with the variant notation
        - web    ; If forward is set the request is forwarded to check
    """
    W = Web.Web()
    C = Config.Config()
    O = Output.Output(__file__, C.Output)

    if not req.form:
        return "Error in input"
    mutationName = req.form.get("mutationName", None)
    variantRecord = req.form.get("variantRecord", None)
    forward = req.form.get("forward", None)

    # We are only interested in the legend
    Mutalyzer.process(mutationName, C, O)

    legends = O.getOutput("legends")

    # Filter the transcript from the legend
    legends = [l for l in legends if "_v" in l[0]]
    for l in legends:
        if l[1] == variantRecord:
            if forward:
                p,a = mutationName.split(':')
                req.form["mutationName"] = p+'('+l[0]+'):'+a
                return check(req)
            else:
                return l[0]
    return "Transcript not found"#+`legends`
#getGS

def Variant_info(req) :
    """
    The I{g.} to I{c.} and vice versa interface for LOVD.

    @arg req: The request:
      - req.form['LOVD_ver'] ; The version of the calling LOVD
      - req.form['build']    ; The human genome build (hg19 assumed)
      - req.form['acc']      ; The accession number (NM number)
      - req.form['var']      ; A description of the variant
    @type req: object

    @return: An HTML page containing the results of Variant_info
    @rtype: string
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


def __checkInt(inpv, refname):
    """
    @arg inpv:
    @type inpv:
    @arg refname:
    @type refname:
    
    @todo: documentation
    """
    #remove , . and -
    inpv = inpv.replace(',','').replace('.','').replace('-','')
    try:
        return int(inpv)
    except ValueError, e:
        raise InputException("Expected an integer in field: %s" % refname)

def upload(req) :
    """
    @arg req:
    @type req:
    
    @return:
    @rtype:
    
    @todo: documentation
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
