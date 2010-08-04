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
    #TODO documentation
    """
        Debug Wrapper for functions called from within the daemon
    """

    def _tempf(*args):
        #TODO documentation
        """
        """

        of = open("/tmp/daemon.out", "a+")
        try:
            of.write("\nFunction %s\n\targs: %s\n\t" % (`f`, `args`))
            ret = f(*args)
            of.write("Returns: %s" % `ret`)
            return ret
        #try
        except Exception, e:
            import traceback
            of.write("\nEXCEPTION:\n")
            traceback.print_exc(file=of)
        #except
    return _tempf
#debug


class Scheduler() :
    #TODO documentation
    """
        Special methods:
            __init__(config, database) ;
        
        Public methods:
            isDaemonRunning() ;
            process()         ;
            addJob(outputFilter, eMail, queue, fromHost) ;

    """

    def __init__(self, config, database) :
        #TODO documentation
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

    def __processFlags(self, O, flags):
        #TODO documentation
        """
        """

        #TODO: Add more info for the user 
        if not flags: return
        if 'S' in flags: #This entry is going to be skipped
            #Add a usefull message to the Output object
            O.addMessage(__file__, 4, "EBSKIP", "Skipping entry")
            return True #skip
        #if
        if 'A' in flags: #This entry is altered before execution
            O.addMessage(__file__, 3, "WEALTER", "Entry altered before "
                    "execution")
    #__processFlags                    

    def __alterBatchEntries(self, jobID, old, new, flag, nselector):
        #TODO documentation
        """
        """

        self.__database.updateBatchDb(jobID, old, new, flag, nselector)
    #__alterBatchEntries

    def __skipBatchEntries(self, jobID, flag, selector):
        #TODO documentation
        """
        """

        self.__database.skipBatchDb(jobID, selector, flag)
    #__skipBatchEntries

    def _updateDbFlags(self, O, jobID):
        #TODO documentation
        """
        """

        flags = O.getOutput("BatchFlags")
        if not flags: return
        #First check if we need to skip
        for flag, args in flags:
            if 'S' in flag:
                selector = args
                O.addMessage(__file__, 3, "WBSKIP",
                        "All further occurrences with '%s' will be "
                        "skipped" % selector)
                self.__skipBatchEntries(jobID, flag, selector)
                return
            #if
        #for
        #If not skipflags, check if we need to alter
        for flag, args in flags:
            if 'A' in flag:
                old, new, nselector = args  #nselector = negative selector
                O.addMessage(__file__, 3, "WBSUBST",
                        "All further occurrences of %s will be substituted "
                        "by %s" % (old, new))
                self.__alterBatchEntries(jobID, old, new, flag, nselector)
            #if
        #for
    #_updateDbFlags

    def process(self) :
        #TODO documentation
        """
            Start the mutalyzer Batch Processing
        """

        jobList = self.__database.getJobs()

        while jobList :
            for i, jobType, arg1 in jobList :
                inputl, flags = self.__database.getFromQueue(i)
                if inputl:
                    if jobType == "NameChecker":
                        self._processNameBatch(inputl, i, flags)
                    elif jobType == "SyntaxChecker":
                        self._processSyntaxCheck(inputl, i)
                    elif jobType == "ConversionChecker":
                        self._processConversion(inputl, i, arg1)
                    else: #unknown jobType
                        pass #TODO: Scream burning water
                else :
                    eMail, stuff, fromHost = self.__database.removeJob(i)
                    print "Job %s finished, email %s file %s" % (i, eMail, i)
                    self.__sendMail(eMail, "%sResults_%s.txt" % (fromHost, i))
                #else
            #for
            jobList = self.__database.getJobs()
        #while
    #process

    def _processNameBatch(self, cmd, i, flags):
        #TODO documentation
        """
            Process an entry from the Name Batch,

            write the results to the job-file
        """

        C = Config.Config()
        O = Output.Output(__file__, C.Output)

        #Read out the flags
        skip = self.__processFlags(O, flags)

        if not skip:
            #Run mutalyzer and get values from Output Object 'O'
            try:
                Mutalyzer.process(cmd, C, O)
            except Exception, e:
                #Catch all exceptions related to the processing of cmd
                O.addMessage(__file__, -1, "EBATCH",
                        "Error during NameChecker Batch. Input: %s" % `cmd`)
                O.addMessage(__file__, 4, "EBATCHU",
                        "Unexpected error occurred, dev-team notified")
                #import traceback
                #O.addMessage(__file__, 4, "DEBUG", `traceback.format_exc()`)
            #except
            finally:
                #check if we need to update the database
                self._updateDbFlags(O, i)
        #if

        batchOutput = O.getOutput("batchDone")

        outputline =  "%s\t" % cmd
        outputline += "%s\t" % "|".join(O.getBatchMessages(3))

        if batchOutput:
            outputline += batchOutput[0]

        outputline += "\n"

        #Output
        filename = "%s/Results_%s.txt" % (self.__config.resultsDir, i)
        if not os.path.exists(filename):
            handle = open(filename, 'a')
            handle.write("%s\n" % "\t".join(i for i in #TODO clarify
                self.__config.nameCheckOutHeader))
        #if
        else:
            handle = open(filename, 'a')

        handle.write(outputline)
        handle.close()
    #_processNameBatch

    def _processSyntaxCheck(self, cmd, i):
        #TODO documentation
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
        filename = "%s/Results_%s.txt" % (self.__config.resultsDir, i)
        if not os.path.exists(filename):
            handle = open(filename, 'a')
            handle.write("%s\n" % "\t".join(i for i in #TODO clarify
                self.__config.syntaxCheckOutHeader))
        #if
        else:
            handle = open(filename, 'a')

        handle.write("%s\t%s\n" % (cmd, result))
        handle.close()
    #_processSyntaxCheck

    def _processConversion(self, cmd, i, build):
        #TODO documentation
        """
            _processConversion docstring
        """

        C = Config.Config()
        O = Output.Output(__file__, C.Output)
        variant = cmd
        variants = None
        gName = ""
        cNames = [""]

        try:
            #process
            converter = Mapper.Converter(build, C, O)

            #Also accept chr accNo
            variant = converter.correctChrVariant(variant)

            if not (":c." in variant or ":g." in variant):
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
                    gName = variant #TODO clarify
                    cNames = [cName for cName2 in variants.values() for cName in
                            cName2]
                #if
            #if
        #try
        except Exception, e:
            #Catch all exceptions related to the processing of cmd
            O.addMessage(__file__, -1, "EBATCH",
                    "Error during ConversionBatch. Input: %s" % `cmd`)
            O.addMessage(__file__, 4, "EBATCHU",
                    "Unexpected error occurred, dev-team notified")
        #except

        error = "%s" % "|".join(O.getBatchMessages(3))

        #Output
        filename = "%s/Results_%s.txt" % (self.__config.resultsDir, i)
        if not os.path.exists(filename):
            handle = open(filename, 'a')
            handle.write("%s\n" % "\t".join(i for i in
                self.__config.positionConverterOutHeader))
        #if
        else:
            handle = open(filename, 'a')

        handle.write("%s\t%s\t%s\t%s\n" % (cmd, error, gName, "\t".join(cNames)))
        handle.close()
    #_processConversion



    def addJob(self, outputFilter, eMail, queue, fromHost, jobType, Arg1) :
        #TODO documentation
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
        return jobID
    #addJob
#Scheduler

#
# Unit test.
#
if __name__ == "__main__" :
    pass
#if
