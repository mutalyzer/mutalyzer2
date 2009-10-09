#!/usr/bin/python

class batchJob(object) :
    """
        A batchJob object, to store a list of jobs and callback information for
        when the job is finished. 

        Public variables:
            queue ; The batch queue for this job.
            hook  ; The callback information for this job.
    """

    def __init__(self, queue, hook) :
        """
            Initialise the class.

            Public variables (altered):
                queue ; The batch queue for this job.
                hook  ; The callback information for this job.
        """

        self.queue = queue
        self.hook = hook
    #__init()
#batchJob

class Scheduler() :
    """
        The job scheduler consists of one queue for interactive jobs and
        a queue of queues for batch jobs. 
        - Each time an interactive job is scheduled, it is appended to the
          interactive queue, each time a batch job is scheduled, the batch
          (which is a list of jobs) is placed as a whole in the queue for batch
          jobs.
        - Each time a job is requested for processing, each other time an
          interactive job is selected. If a batch job is selected, the
          particular batch queue is selected in a round robin way.
        - Is a batch job is finished, the provided callback information is 
          used.
        This way, interactive jobs have a high priority and all batch jobs
        have a priority linear to their size, so a short batch job will finish
        after a short time, even is there is a large batch job running as well.

        Special methods:
            __init__() ; Initialise the class.

        Private variables:
            __interactive ; A fifo queue for interactive jobs.
            __batch       ; A round-robin queue of fifo queues for batch jobs.
            __turn        ; A switch to see who's turn it is, "interactive" or
                            "batch".
            __batchturn   ; A switch to see which batch queue's turn it is.
            __sleeptime   ; Time to sleep if the scheduler is idle.
            
        Private methods:
            __switch()      ; Switch between interactive jobs and batch jobs.
            __getfifo(list) ; Get a job from a queue.

        Public methods:
            queue(job)       ; Put an interactive job in the queue.
            batchqueue(list) ; Put a list of jobs in the batchqueue.
            getjob()         ; Get the next job.
        
    """

    def __init__(self) :
        """
            Initialise the class.

            Private variables (altered):
                __interactive ; The fifo queue for interactive jobs.
                __batch       ; The round-robin queue of fifo queues for batch
                                jobs.
                __turn        ; A switch to see who's turn it is, "interactive"
                                or "batch".
                __batchturn   ; A switch to see which batch queue's turn it is.
                __sleeptime   ; Time to sleep if the scheduler is idle.
        """

        self.__interactive = []
        self.__batch = []
        self.__turn = "batch"
        self.__batchturn = 0
        self.__lastret = None
        self.__sleeptime = 1
    #__init__

    def __switch(self) :
        """
            Select the next job to be processed.
            If __turn is "interactive", switch to "batch" and select the next
            batch queue in a round robin way.
            If __turn is "batch", switch to "interactive".

            Private variables (altered):
                __turn      ; Flipped from "interactive" to "batch" or vice
                              versa.
                __batchturn ; Is increased (modulo the number of batch queues)
                              if __turn is switched to "batch".
        """

        if self.__turn == "interactive" and len(self.__batch):
            self.__turn = "batch"
            self.__batchturn = (self.__batchturn + 1) % len(self.__batch)
        #if
        else :
            self.__turn = "interactive"
    #__switch

    def __getfifo(self, list) :
        """
            Get the first element in a queue and remove it from the queue.

            Arguments:
                list ; A queue of jobs, this list will be altered (the first
                       element is deleted).

            Returns:
                job ; The first element of the list.
        """

        ret = None

        if list :
            ret = list[0]
            list.pop(0)
        #if

        return ret
    #__getfifo

    def queue(self, job) :
        """
            Schedule an interactive job.

            Arguments:
                job ; The job to be placed in the interactive queue.

            Private variables (altered):
                __interactive ; The job is appended to this list.
        """

        self.__interactive.append(job)
    #queue

    def batchqueue(self, jobList) :
        """
            Schedule a list of jobs (a batch job).

            Arguments:
                jobList ; The list of jobs to be placed in the batch queue.

            Private variables (altered):
                __batch ; The list of jobs is appended to this list.
        """

        #if (jobList) :
        self.__batch.append(jobList)
    #batchqueue

    def getjob(self) :
        """
            Request a job for processing. If both __lastret (the previous 
            result) and ret (the current result) are None, all queues are
            empty, so we can sleep for a while.
            Make sure to call this function from a separate thread, otherwise
            it may block (for __sleeptime seconds).
           
            Returns:
                job ; The job which is first in line.

            Private variables:
                __sleeptime ; The time to sleep if there is nothing to do.

            Private variables (altered):
                __interactive ; May be altered by __getfifo().
                __batch       ; If __getfifo() returns None, then the empty
                                list is removed from the queue of batch queues
                                (the batch job is finished).
                __lastret     ; Set to the current result.
                __sleeptime   ; Seconds to sleep.
        """

        import time

        self.__switch() # Select the next job in line.
        if self.__turn == "interactive" :
            ret =  self.__getfifo(self.__interactive)
        else :
            ret = self.__getfifo(self.__batch[self.__batchturn - 1].queue)
            # Remove empty (finished) batch queues.
            if not self.__batch[self.__batchturn - 1].queue :
                # Of course, do more than printing the hook...
                print self.__batch[self.__batchturn - 1].hook
                self.__batch.pop(self.__batchturn - 1)
        #else

        # Sleep for __sleeptime if we have nothing to do.
        if not self.__lastret and not ret :
            time.sleep(self.__sleeptime)
        self.__lastret = ret

        return ret
    #getjob
#Scheduler
"""
#!/usr/bin/python

import os
import time
import threading

shell_command = "grep -c \"^processor\" /proc/cpuinfo"

cpus = int(os.popen(shell_command).read())
inuse = 0

class hmmm(threading.Thread) :
    def __init__(self, i) :
        self.i = i
        threading.Thread.__init__(self)
    #__init__

    def run(self) :
        global inuse
    
        print time.strftime("%H:%M:%S") + ' ' + str(self.i) + " starting"
        time.sleep(i + 1)
        print time.strftime("%H:%M:%S") + ' ' + str(self.i) + " ends"
        inuse -= 1
    #__run__

for i in range(10) :
    while inuse >= cpus :
        time.sleep(0.01)
    inuse += 1
    bla = hmmm(i)
    bla.start()
#for
"""
