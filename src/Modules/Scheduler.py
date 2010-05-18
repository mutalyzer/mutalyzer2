#!/usr/bin/python

import time

import Config
import Db
import subprocess
import psutil
import os

class Scheduler() :
    """
    """

    def __init__(self, config, database) :
        """
        """

        self.__config = config
        self.__database = database
    #__init__

    def isDaemonRunning(self) :
        """
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
                    print i, results
                else :
                    eMail = self.__database.removeJob(i)
                    print "Job %s finished, email %s file %s" % (i, eMail, i)
                time.sleep(1)
            #for
            jobList = self.__database.getJobs()
        #while
    #process

    def addJob(self, outputFilter, eMail, queue) :
        """
        """

        print "called addjob"
        jobList = self.__database.getJobs()
        jobID = self.__database.addJob(outputFilter, eMail)
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
