#!/usr/bin/python

"""
    Public classes:
        Scheduler ;
"""

import time

import subprocess
import psutil
import os
import smtplib
from email.mime.text import MIMEText

from Modules import Config
from Modules import Output
from Modules import Db

import Mutalyzer

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
        #TODO: Handle Connection errors in a try, except fasion
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
            for i, jobType in jobList :
                results = self.__database.getFromQueue(i)
                if results:
                    if jobType == "NameChecker":
                        self._processNameBatch(results, i)
                    elif jobType == "SyntaxChecker":
                        self._processSyntaxCheck(results, i)
                    elif jobType == "Conversions":
                        self._processConversion(results, i)
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

    def _processNameBatch(self, results, i):
        """
            Process an entry from the Name Batch,

            write the results to the job-file
        """
        C = Config.Config()
        O = Output.Output(__file__, C.Output)

        #create NameCheck cmd
        if results[1] :
            if results[0].startswith("LRG"):
                cmd = "%s%s:%s" % results
            else:
                cmd = "%s(%s):%s" % results
        else :
            cmd = "%s:%s" % (results[0], results[2])

        #Run mutalyzer and get values from Output Object 'O'
        Mutalyzer.process(cmd, C, O)
        batchData = O.getOutput("batch")    #what if batch does not exist
        if batchData:
            bO = batchData[0]
        else:
            bO = None

        #AccNo, GeneSymbol, Mutation
        outputline =  "%s\t%s\t%s\t" % results[:3]

        if bO:
            g = bO["gName"]                     #Genomic Plain Description
            c = bO["cName"]                     #Coding  Plain Description
            p = bO["pName"]                     #Protein Plain Description
            gc = "%s:%s" % (bO["cSymbol"], c)   #Coding  GeneSymbol Description
            gp = "%s:%s" % (bO["pSymbol"], p)   #Protein GeneSymbol Description
            ag = results[0]                     #Genomic Accesion Description
            ac = bO["cAcc"] or ag               #Coding  Accesion Description
            ap = bO["pAcc"] or ag               #Protein Accesion Description

            #plain          gName, cName, pName
            outputline += "%s\t%s\t%s\t" % (g,c,p)

            #genesymbol     cName, pName
            outputline += "%s\t%s\t" % (gc, gp)

            #AccNo          gName, cName, pName
            outputline += "%s\t%s\t%s\t" % (ag, ac, ap)

            #Messages & Newline
            outputline += "{%s}\n" % "|".join(bO["messages"])

       #if
        else:
            #Parsing went wrong, skip fields
            outputline += "\t"*9

        #Add Debug Data
        outputline += "DEBUG:\n\t%s\n" % "\n\t".join(O.getMessages())


        filename = "%s/Results_%s.txt" % (self.__config.resultsDir, i)
        #TODO: If path not exists create file with correct header
        handle = open(filname, 'a')
        handle.write(outputline)
        handle.close()

    #_processNameBatch

    def _processSyntaxCheck(self, results, i):
        """
            _processSyntaxCheck docstring
        """
        pass
    #_processSyntaxCheck

    def _processConversion(self, results, i):
        """
            _processConversion docstring
        """
        pass
    #_processConversion



    def addJob(self, outputFilter, eMail, queue, fromHost, jobType) :
        """
            Arguments:
                outputFilter ; Filter the output of Mutalyzer
                eMail        ; e-mail address of batch supplier
                queue        ; A list of jobs
                fromHost     ; From where is the request made
                jobType      ; The type of Batch Job that should be run
        """

        print "called addjob"
        jobList = self.__database.getJobs()
        jobID = self.__database.addJob(outputFilter, eMail, fromHost, jobType)
        for i in queue :
            self.__database.addToQueue(jobID, *i)
        #FIXME!:    This will not work, because Popen doesn't detach from the
        #           parent, so on return, the process is interrupted. BAD!
        #           Use multiprocessing, check if PID catch still works
        subprocess.Popen([self.__config.processName, "src/BatchChecker.py"],
                         executable = "python")
        time.sleep(1)
    #addJob
#Scheduler

#
# Unit test.
#
if __name__ == "__main__" :
    pass
#if
