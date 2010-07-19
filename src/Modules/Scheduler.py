#!/usr/bin/python

"""
    Public classes:
        Scheduler ;
"""

#from __future__ import with_statement

import time

import subprocess
import psutil
import os
import smtplib
from email.mime.text import MIMEText
import daemon

from Modules import Config
from Modules import Output
from Modules import Db
from Modules import Parser
from Modules import Mapper

import Mutalyzer
#import BatchChecker

def debug(f):
    """Debug Wrapper for functions called from within the daemon"""
    def _tempf(*args):
        of = open("/tmp/daemon.out", "a+")
        try:
            of.write("\nFunction %s\n\targs: %s\n\t" % (`f`, `args`))
            ret = f(*args)
            of.write("Returns: %s" % `ret`)
            return ret
        except Exception, e:
            import traceback
            of.write("\nEXCEPTION:\n")
            traceback.print_exc(file=of)
    return _tempf


class Scheduler() :
    """
        Special methods:
            __init__(config, database) ;
        
        Public methods:
            isDaemonRunning() ;
            process()         ;
            addJob(outputFilter, eMail, queue, fromHost) ;

    """

    def __init__(self, config, database) :
        """
            Arguments:
                config   ;
                database ;
        """

        self.__config = config
        self.__database = database
    #__init__

    def __sendMail(self, mailTo, url) :
        """
            Send an e-mail containing an url to a batch job submitter.

            Arguments:
                mailTo ; The batch job submitter.
                url    ; The url containing the results.

            Private variables:
                __config ; The variables mailMessage, mailSubject and mailFrom
                           are used.
        """
        #TODO: Handle Connection errors in a try, except fassion
        #socket.error 

        handle = open(self.__config.mailMessage)
        message = MIMEText(handle.read() % url)
        handle.close()

        message["Subject"] = self.__config.mailSubject
        message["From"] = self.__config.mailFrom
        message["To"] = mailTo

        smtpInstance = smtplib.SMTP()
        smtpInstance.connect()
        smtpInstance.sendmail(self.__config.mailFrom, mailTo,
                              message.as_string())
        smtpInstance.quit()
    #__sendMail

    def isDaemonRunning(self) :
        """
            Returns:
                True if an other scheduler is already running, False otherwise.
        """

        myPid = os.getpid()
        for i in psutil.get_process_list() :
            if i.cmdline and i.cmdline[0] == self.__config.processName  and \
               i.pid != myPid :
                return True
        return False
    #isDaemonRunning

    def process(self) :
        """
        """
        jobList = self.__database.getJobs()

        while jobList :
            for i, jobType, arg1 in jobList :
                inputl = self.__database.getFromQueue(i)
                if inputl:
                    if jobType == "NameChecker":
                        self._processNameBatch(inputl, i)
                    elif jobType == "SyntaxChecker":
                        self._processSyntaxCheck(inputl, i)
                    elif jobType == "ConversionChecker":
                        self._processConversion(inputl, i, arg1)
                    else: #unknown jobType
                        pass #TODO: Scream burning water
                else :
                    eMail, stuff, fromHost = self.__database.removeJob(i)
                    print "Job %s finished, email %s file %s" % (i, eMail, i)
                    self.__sendMail(eMail, "%sResults_%s.csv" % (fromHost, i))
                #else
            #for
            jobList = self.__database.getJobs()
        #while
    #process

    def _processNameBatch(self, cmd, i):
        """
            Process an entry from the Name Batch,

            write the results to the job-file
        """
        #FIXME: This output is generated in an UUUUUUGLY way, see Mutalyzer.py
        # Preferably receive the 12 output fields in tact
        C = Config.Config()
        O = Output.Output(__file__, C.Output)


        #Run mutalyzer and get values from Output Object 'O'
        Mutalyzer.process(cmd, C, O)
        batchData = O.getOutput("batch")
        if batchData:
            bO = batchData[0]
        else:
            bO = {}

        #AccNo, GeneSymbol, Mutation FIXME for oneliners
        outputline =  "%s\t%s\t%s\t" % (cmd, "", "") #split up the cmd

        if bO and bO["exitcode"] == 0:
            g = bO["gName"]                     #Genomic Plain Description
            c = bO["cName"]                     #Coding  Plain Description
            p = bO["pName"]                     #Protein Plain Description
            gc = "%s:%s" % (bO["cSymbol"], c)   #Coding  GeneSymbol Description
            gp = "%s:%s" % (bO["pSymbol"], p)   #Protein GeneSymbol Description
            ag = cmd.split(':')[0].split('(')[0]#Genomic Accesion Description
            ac = bO["cAcc"] or ag               #Coding  Accesion Description
            ap = bO["pAcc"] or ag               #Protein Accesion Description

            #plain          gName, cName, pName
            outputline += "%s\t%s\t%s\t" % (g,c,p)

            #genesymbol     cName, pName
            outputline += "%s\t%s\t" % (gc, gp)

            #AccNo          gName, cName, pName
            outputline += "%s\t%s\t%s\t" % (ag, ac, ap)

        #if
        else:
            #Something went wrong, skip fields
            outputline += "\t"*8

        outputline += "%s\n" % "|".join(O.getBatchMessages(3)+\
                bO.get("messages",[]))

        #Output
        filename = "%s/Results_%s.csv" % (self.__config.resultsDir, i)
        if not os.path.exists(filename):
            handle = open(filename, 'a')
            handle.write("%s\n" % "\t".join(i for i in
                self.__config.nameCheckOutHeader))
        else:
            handle = open(filename, 'a')

        handle.write(outputline)
        handle.close()

    #_processNameBatch

    def _processSyntaxCheck(self, cmd, i):
        """
            _processSyntaxCheck docstring
        """
        C = Config.Config()
        O = Output.Output(__file__, C.Output)
        P = Parser.Nomenclatureparser(O)

        #Process
        parsetree = P.parse(cmd)
        if parsetree:
            result = "OK"
        else:
            result = "|".join(O.getBatchMessages(3))

        #Output
        filename = "%s/Results_%s.csv" % (self.__config.resultsDir, i)
        if not os.path.exists(filename):
            handle = open(filename, 'a')
            handle.write("%s\n" % "\t".join(i for i in
                self.__config.syntaxCheckOutHeader))
        else:
            handle = open(filename, 'a')

        handle.write("%s\t%s\n" % (cmd, result))
        handle.close()
    #_processSyntaxCheck

    def _processConversion(self, cmd, i, build):
        """
            _processConversion docstring
        """
        C = Config.Config()
        O = Output.Output(__file__, C.Output)
        variant = cmd

        #process
        converter = Mapper.Converter(build, C, O)

        variants = None
        gName = ""
        cNames = [""]

        #Also accept chr accNo
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
                gName = variant
                cNames = [cName for cName2 in variants.values() for cName in
                        cName2]

        error = "%s" % "|".join(O.getBatchMessages(3))

        #Output
        filename = "%s/Results_%s.csv" % (self.__config.resultsDir, i)
        if not os.path.exists(filename):
            handle = open(filename, 'a')
            handle.write("%s\n" % "\t".join(i for i in
                self.__config.positionConverterOutHeader))
        else:
            handle = open(filename, 'a')

        handle.write("%s\t%s\t%s\t%s\n" % (cmd, error, gName, "\t".join(cNames)))
        handle.close()
    #_processConversion



    def addJob(self, outputFilter, eMail, queue, fromHost, jobType, Arg1) :
        """
            Arguments:
                outputFilter ; Filter the output of Mutalyzer
                eMail        ; e-mail address of batch supplier
                queue        ; A list of jobs
                fromHost     ; From where is the request made
                jobType      ; The type of Batch Job that should be run
        """
        # Add jobs to the database
        jobID = self.__database.addJob(outputFilter, eMail,
                fromHost, jobType, Arg1)
        for inputl in queue :
            self.__database.addToQueue(jobID, inputl)

        # Spawn child
        p = subprocess.Popen(["MutalyzerBatch",
            "src/BatchChecker.py"], executable="python")

        #Wait for the BatchChecker to fork of the Daemon
        p.communicate()
    #addJob
#Scheduler

#
# Unit test.
#
if __name__ == "__main__" :
    pass
#if
