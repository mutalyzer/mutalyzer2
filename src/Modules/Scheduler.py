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
            for i in jobList :
                results = self.__database.getFromQueue(i)
                if results :
                    if results[1] :
                        cmd = "%s(%s):%s" % results
                    else :
                        cmd = "%s:%s" % (results[0], results[2])
                    C = Config.Config()
                    O = Output.Output(__file__, C.Output)
                    Mutalyzer.process(cmd, C, O)
                    handle = open("%s/Results_%s.txt" % (
                        self.__config.resultsDir, i), "a")
                    handle.write(str(O.getOutput("variantdescription")))
                    handle.close()
                    del O, C
                #if
                else :
                    eMail, stuff, fromHost = self.__database.removeJob(i)
                    #print "Job %s finished, email %s file %s" % (i, eMail, i)
                    self.__sendMail(eMail, "%sResults_%s.txt" % (fromHost, i))
                #else
            #for
            jobList = self.__database.getJobs()
        #while
    #process

    def addJob(self, outputFilter, eMail, queue, fromHost) :
        """
            Arguments:
                outputFilter ;
                eMail        ;
                queue        ;
        """

        print "called addjob"
        jobList = self.__database.getJobs()
        jobID = self.__database.addJob(outputFilter, eMail, fromHost)
        for i in queue :
            self.__database.addToQueue(jobID, *i)
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
