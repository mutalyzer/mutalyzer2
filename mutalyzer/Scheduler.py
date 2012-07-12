"""
Module used to add and manage the Batch Jobs.

@requires: subprocess
@requires: os
@requires: smtplib
@requires: email.mime.text.MIMEText
@requires: Modules.Config
@requires: Modules.Mapper
@requires: Mutalyzer
"""
# Public classes:
#     - Scheduler ; Manages the batch jobs and contains the methods for
#             - Batch Name Checker
#             - Batch Syntax Checker
#             - Batch Position Converter

import os                               # os.path.exists
import smtplib                          # smtplib.STMP
from email.mime.text import MIMEText    # MIMEText

import mutalyzer
from mutalyzer import config
from mutalyzer import variantchecker
from mutalyzer.grammar import Grammar
from mutalyzer.output import Output
from mutalyzer.mapping import Converter
from mutalyzer import Retriever           # Retriever.Retriever


__all__ = ["Scheduler"]


class Scheduler() :
    """
    Special methods:
        - __init__(database) ;

    Public methods:
        - addJob(outputFilter, eMail, queue, fromHost, jobType, Arg1) ; Add a
          job to the database jobqueue and start the batchChecker daemon.
        - process() ; Iterate over & process the jobs in the jobqueue.

    @summary: Manages the batch jobs and contains the methods for
        - Batch Name Checker
        - Batch Syntax Checker
        - Batch Position Converter
    """

    def __init__(self, database) :
        #TODO: documentation
        """
        Initialize the Scheduler, which requires a database connection.

        @todo: documentation

        @arg database:
        @type database:
        """
        self.__database = database
        self.__run = True
    #__init__

    def stop(self):
        """
        If the {process} method is running, the current job item will be
        processed and {process} will return.
        """
        self.__run = False
    #stop

    def stopped(self):
        """
        Test if the scheduler instance is stopped (i.e. the {stop} method is
        called).

        @return: True if {stop} was called, False otherwise.
        @rtype: bool
        """
        return not self.__run
    #stopped

    def __sendMail(self, mailTo, url) :
        """
        Send an e-mail containing an url to a batch job submitter.

        @todo: Handle Connection errors in a try, except clause

        @arg mailTo: The batch job submitter
        @type mailTo: string
        @arg url: The url containing the results
        @type url: string
        """
        if mutalyzer.is_test():
            return

        # Note: The above check with mutalyzer.is_test is bogus, since during
        # a normal unit test, the batch checker process is not started in the
        # environment of the unit tests.
        # As sort of a hack, we therefore check for the patented address
        # 'test@test.test', used in the unit tests.
        if mailTo == 'test@test.test':
            return

        #TODO: Handle Connection errors in a try, except clause
        #Expected errors: socket.error

        message = MIMEText("""Dear submitter,

The batch operation you have submitted, has been processed successfully.

Your results can be found here:
%s

Thanks for using Mutalyzer.


With kind regards,
Mutalyzer batch checker.""" % url)

        message["Subject"] = config.get('mailSubject')
        message["From"] = config.get('mailFrom')
        message["To"] = mailTo

        smtpInstance = smtplib.SMTP()
        smtpInstance.connect()
        smtpInstance.sendmail(config.get('mailFrom'), mailTo,
                              message.as_string())
        smtpInstance.quit()
    #__sendMail

    def __processFlags(self, O, flags) :
        """
        Translate the flags to error & info messages.

        Side-effect:
            - Added messages to the Output object.

        @arg O: Output object of the current batchentry
        @arg flags: Flags of the current batchentry

        @return: skip ; True if the entry must be skipped
        @rtype: boolean
        """

        if not flags :
            return False
        if 'S' in flags : #This entry is going to be skipped
            #Add a usefull message to the Output object
            if "S0" in flags :
                message = "Entry could not be formatted correctly, check "\
                        "batch input file help for details"
            elif "S9" in flags :
                message = "Empty Line"
            else :
                message = "Skipping entry"
            O.addMessage(__file__, 4, "EBSKIP", message)
            return True #skip
        #if
        if 'A' in flags : #This entry is altered before execution
            O.addMessage(__file__, 3, "WEALTER", "Entry altered before "
                    "execution")
        return False
    #__processFlags

    def __alterBatchEntries(self, jobID, old, new, flag, nselector) :
        """
        Alias for the database.updateBatchDb method.

        Replace within one JobID all entries matching old with new, if
        they do not match the negative selector.

        Example:
        NM_002001(FCER1A_v001):c.1A>C ; this would result in the continuous
        fetching of the reference because no version number is given.
        In this case the arguments would be:
            - old         ;   NM_002001
            - new         ;   NM_002001.2
            - nselector   ;   NM_002001[[.period.]]

        The nselector is used to prevent the replacement of
        false positives. e.g. NM_002001.1(FCER1A_v001):c.1A>C should not
        be replaced. The double bracket notation is the MySQL escape char
        for a regular expression.

        @arg jobID:
        @type jobID:
        @arg old:
        @type old:
        @arg new:
        @type new:
        @arg flag:
        @type flag:
        @arg nselector:
        @type nselector:
        """

        self.__database.updateBatchDb(jobID, old, new, flag, nselector)
    #__alterBatchEntries

    def __skipBatchEntries(self, jobID, flag, selector) :
        """
        Alias for the database.skipBatchDb method.

        Skip all batch entries that match a certain selector.

        @arg jobID:
        @type jobID:
        @arg flag:
        @type flag:
        @arg selector:
        @type selector:
        """

        self.__database.skipBatchDb(jobID, selector, flag)
    #__skipBatchEntries

    def _updateDbFlags(self, O, jobID) :
        """
            Check and set the flags for other entries of jobID.

            After each entry is ran, the Output object can contain BatchFlags.
            If these are set, this means that identical entries need to be
            skipped / altered.

            Side-effect:
               -  Added flags to entries in the database

            @arg O:     Output object of the current batchentry
            @type O:    object
            @arg jobID: ID of job, so that the altering is only done within one
            job
            @type jobID:
        """

        flags = O.getOutput("BatchFlags")
        # NOTE:
        # Flags is a list of tuples. Each tuple consists of a flag and its
        # arguments. A skipped entry has only one argument, the selector
        # E.g. ("S1", "NM_002001.$")
        # An altered entry has three arguments,
        #               old,           new          negative selector
        # E.g.("A2",("NM_002001", "NM_002001.2", "NM_002001[[.period.]]"))

        # Flags are set when an entry could be sped up. This is either the
        # case for the Retriever as for the Mutalyzer module

        if not flags: return
        #First check if we need to skip
        for flag, args in flags :
            if 'S' in flag :
                selector = args     # Strip argument
                O.addMessage(__file__, 3, "WBSKIP",
                        "All further occurrences with '%s' will be "
                        "skipped" % selector)
                self.__skipBatchEntries(jobID, flag, selector)
                return
            #if
        #for
        #If not skipflags, check if we need to alter
        for flag, args in flags :
            if 'A' in flag :
                old, new, nselector = args  #Strip arguments
                O.addMessage(__file__, 3, "WBSUBST",
                        "All further occurrences of %s will be substituted "
                        "by %s" % (old, new))
                self.__alterBatchEntries(jobID, old, new, flag, nselector)
            #if
        #for
    #_updateDbFlags

    def process(self) :
        """
        Start the mutalyzer Batch Processing. This method retrieves
        all jobs from the database and processes them in a roundrobin
        fashion. If all jobs are done the process checks if new jobs are
        added during the last processing round.

        If during this process the {stop} method is called, the current
        job item is completed and we return.

        This method uses two database tables, BatchJob and BatchQueue.

        The jobList is an array of tuples with three elements
            - jobID       ;   The ID of the job
            - jobType     ;   The type of the job
            - argument1   ;   Currently only used for the ConversionChecker
                            to send the build version.

        If the jobList is not empty, the method will iterate once over the
        list and fetch the first entry of a job from the database table
        BatchQueue. This request returns both the input for the batch and
        the flags for the job.

        #Flags
        A job can be flagged in three ways:
          - A       ;   Altered - this means that the input is altered
                        before execution. This could be the case if an
                        entry uses an accession number without a version.
                        If a version is retrieved from the NCBI, all
                        further occurences of that accession will be
                        replaced by the accession with version number.
          - S       ;   Skipped - this means that this batchentry will be
                        skipped by the batchprocess. This could be the
                        case if the user made a mistake that could not be
                        auto fixed and henceforth all occurences of the
                        mistake will be skipped.
          - C       ;   Continue - this means the input does not end the
                        current row, so no new row in the output should
                        be started.

        A Flag consists of either an A, S or C followed by a digit, which
        refers to the reason of alteration / skip.
        """
        jobList = self.__database.getJobs()

        while jobList and self.__run:
            for i, jobType, arg1 in jobList :
                inputl, flags = self.__database.getFromQueue(i)
                if not (inputl is None) :
                    if jobType == "NameChecker" :
                        self._processNameBatch(inputl, i, flags)
                    elif jobType == "SyntaxChecker" :
                        self._processSyntaxCheck(inputl, i, flags)
                    elif jobType == "PositionConverter" :
                        self._processConversion(inputl, i, arg1, flags)
                    elif jobType == "SnpConverter" :
                        self._processSNP(inputl, i, flags)
                    else: #unknown jobType
                        pass #TODO: Scream burning water and remove from Queue
                else :
                    eMail, stuff, fromHost = self.__database.removeJob(i)
                    print "Job %s finished, email %s file %s" % (i, eMail, i)
                    self.__sendMail(eMail, "%sResults_%s.txt" % (fromHost, i))
                #else
                if not self.__run:
                    break
            #for
            jobList = self.__database.getJobs()
        #while
    #process

    def _processNameBatch(self, cmd, i, flags) :
        """
        Process an entry from the Name Batch, write the results
        to the job-file. If an Exception is raised, catch and continue.

        Side-effect:
            - Output written to outputfile.

        @arg cmd: The NameChecker input
        @type cmd:
        @arg i: The JobID
        @type i:
        @arg flags: Flags of the current entry
        @type flags:
        """
        O = Output(__file__)
        O.addMessage(__file__, -1, "INFO",
            "Received NameChecker batchvariant " + cmd)

        #Read out the flags
        skip = self.__processFlags(O, flags)

        if not skip :
            #Run mutalyzer and get values from Output Object 'O'
            try :
                variantchecker.check_variant(cmd, O)
            except Exception:
                #Catch all exceptions related to the processing of cmd
                O.addMessage(__file__, 4, "EBATCHU",
                        "Unexpected error occurred, dev-team notified")
                import traceback
                O.addMessage(__file__, 4, "DEBUG", repr(traceback.format_exc()))
            #except
            finally :
                #check if we need to update the database
                self._updateDbFlags(O, i)
        #if

        batchOutput = O.getOutput("batchDone")

        outputline =  "%s\t" % cmd
        outputline += "%s\t" % "|".join(O.getBatchMessages(3))

        if batchOutput :
            outputline += batchOutput[0]

        #Output
        filename = "%s/Results_%s.txt" % (config.get('resultsDir'), i)
        if not os.path.exists(filename) :
            # If the file does not yet exist, create it with the correct
            # header above it. The header is read from the config file as
            # a list. We need a tab delimited string.
            header = config.get('nameCheckOutHeader')
            handle = open(filename, 'a')
            handle.write("%s\n" % "\t".join(header))
        #if
        else :
            handle = open(filename, 'a')

        if flags and 'C' in flags:
            separator = '\t'
        else:
            separator = '\n'

        handle.write("%s%s" % (outputline, separator))
        handle.close()
        O.addMessage(__file__, -1, "INFO",
            "Finished NameChecker batchvariant " + cmd)
    #_processNameBatch

    def _processSyntaxCheck(self, cmd, i, flags) :
        """
        Process an entry from the Syntax Check, write the results
        to the job-file.

        Side-effect:
            - Output written to outputfile

        @arg cmd:   The Syntax Checker input
        @type cmd:
        @arg i:     The JobID
        @type i:
        @arg flags: Flags of the current entry
        @type flags:
        """
        output = Output(__file__)
        grammar = Grammar(output)

        output.addMessage(__file__, -1, "INFO",
                           "Received SyntaxChecker batchvariant " + cmd)

        skip = self.__processFlags(output, flags)
        #Process
        if not skip :
            parsetree = grammar.parse(cmd)
        else :
            parsetree = None

        if parsetree :
            result = "OK"
        else :
            result = "|".join(output.getBatchMessages(3))

        #Output
        filename = "%s/Results_%s.txt" % (config.get('resultsDir'), i)
        if not os.path.exists(filename) :
            # If the file does not yet exist, create it with the correct
            # header above it. The header is read from the config file as
            # a list. We need a tab delimited string.
            header = config.get('syntaxCheckOutHeader')
            handle = open(filename, 'a')
            handle.write("%s\n" % "\t".join(header))
        #if
        else :
            handle = open(filename, 'a')

        if flags and 'C' in flags:
            separator = '\t'
        else:
            separator = '\n'

        handle.write("%s\t%s%s" % (cmd, result, separator))
        handle.close()
        output.addMessage(__file__, -1, "INFO",
                          "Finished SyntaxChecker batchvariant " + cmd)
    #_processSyntaxCheck

    def _processConversion(self, cmd, i, build, flags) :
        """
        Process an entry from the Position Converter, write the results
        to the job-file. The Position Converter is wrapped in a try except
        block which ensures that he Batch Process keeps running. Errors
        are caught and the user will be notified.

        Side-effect:
            - Output written to outputfile.

        @arg cmd: The Syntax Checker input
        @type cmd: string
        @arg i: The JobID
        @type i: integer
        @arg build: The build to use for the converter
        @type build: string
        @arg flags: Flags of the current entry
        @type flags:
        """
        O = Output(__file__)
        variant = cmd
        variants = None
        gName = ""
        cNames = [""]

        O.addMessage(__file__, -1, "INFO",
            "Received PositionConverter batchvariant " + cmd)

        skip = self.__processFlags(O, flags)
        if not skip :
            try :
                #process
                converter = Converter(build, O)

                #Also accept chr accNo
                variant = converter.correctChrVariant(variant)

                #TODO: Parse the variant and check for c or g. This is ugly
                if not(":c." in variant or ":n." in variant or ":g." in variant) :
                    #Bad name
                    grammar = Grammar(O)
                    grammar.parse(variant)
                #if

                if ":c." in variant or ":n." in variant :
                    # Do the c2chrom dance
                    variant = converter.c2chrom(variant)
                    # NOTE:
                    # If we received a coding reference convert that to the
                    # genomic position variant. Use that variant as the input
                    # of the chrom2c.

                # If the input is a genomic variant or if we converted a
                # coding variant to a genomic variant we try to find all
                # other affected coding variants.
                if variant and ":g." in variant :
                    # Do the chrom2c dance
                    variants = converter.chrom2c(variant, "dict")
                    if variants :
                        gName = variant
                        # Due to the cyclic behavior of the Position Converter
                        # we know for a fact that if a correct chrom name is
                        # generated by the converter.c2chrom that we will at
                        # least find one variant with chrom2c. Collect the
                        # variants from a nested lists and store them.
                        cNames = [cName for cName2 in variants.values() \
                                for cName in cName2]
            except Exception:
                #Catch all exceptions related to the processing of cmd
                O.addMessage(__file__, 4, "EBATCHU",
                        "Unexpected error occurred, dev-team notified")
            #except
        #if

        error = "%s" % "|".join(O.getBatchMessages(3))

        #Output
        filename = "%s/Results_%s.txt" % (config.get('resultsDir'), i)
        if not os.path.exists(filename) :
            # If the file does not yet exist, create it with the correct
            # header above it. The header is read from the config file as
            # a list. We need a tab delimited string.
            header = config.get('positionConverterOutHeader')
            handle = open(filename, 'a')
            handle.write("%s\n" % "\t".join(header))
        #if
        else :
            handle = open(filename, 'a')

        if flags and 'C' in flags:
            separator = '\t'
        else:
            separator = '\n'

        handle.write("%s\t%s\t%s\t%s%s" % (cmd, error, gName, "\t".join(cNames), separator))
        handle.close()
        O.addMessage(__file__, -1, "INFO",
            "Finisehd PositionConverter batchvariant " + cmd)
    #_processConversion


    def _processSNP(self, cmd, i, flags) :
        """
        Process an entry from the SNP converter Batch, write the results
        to the job-file. If an Exception is raised, catch and continue.

        Side-effect:
            - Output written to outputfile.

        @arg cmd: The SNP converter input
        @type cmd:
        @arg i: The JobID
        @type i:
        @arg flags: Flags of the current entry
        @type flags:
        """
        O = Output(__file__)
        O.addMessage(__file__, -1, "INFO",
            "Received SNP converter batch rs" + cmd)

        #Read out the flags
        # Todo: Do something with the flags?
        skip = self.__processFlags(O, flags)

        descriptions = []
        if not skip :
            R = Retriever.Retriever(O, None)
            descriptions = R.snpConvert(cmd)

        # Todo: Is output ok?
        outputline =  "%s\t" % cmd
        outputline += "%s\t" % "|".join(descriptions)
        outputline += "%s\t" % "|".join(O.getBatchMessages(3))

        #Output
        filename = "%s/Results_%s.txt" % (config.get('resultsDir'), i)
        if not os.path.exists(filename) :
            # If the file does not yet exist, create it with the correct
            # header above it. The header is read from the config file as
            # a list. We need a tab delimited string.
            header = config.get('snpConverterOutHeader')
            handle = open(filename, 'a')
            handle.write("%s\n" % "\t".join(header))
        #if
        else :
            handle = open(filename, 'a')

        if flags and 'C' in flags:
            separator = '\t'
        else:
            separator = '\n'

        handle.write("%s%s" % (outputline, separator))
        handle.close()
        O.addMessage(__file__, -1, "INFO",
                     "Finished SNP converter batch rs%s" % cmd)
    #_processSNP


    def addJob(self, outputFilter, eMail, queue, columns, fromHost, jobType,
               Arg1):
        """
        Add a job to the Database and start the BatchChecker.

        @arg outputFilter:  Filter the output of Mutalyzer
        @type outputFilter:
        @arg eMail:         e-mail address of batch supplier
        @type eMail:        string
        @arg queue:         A list of jobs
        @type queue:        list
        @arg columns:       The number of columns.
        @type columns:      int
        @arg fromHost:      From where is the request made
        @type fromHost:
        @arg jobType:       The type of Batch Job that should be run
        @type jobType:
        @arg Arg1:          Batch Arguments, for now only build info
        @type Arg1:

        @return: jobID
        @rtype:

        @todo: outputFilter is not used
        """
        #TODO: outputFilter is not used

        # Add jobs to the database
        jobID = self.__database.addJob(outputFilter, eMail, fromHost, jobType,
                                       Arg1)

        for i, inputl in enumerate(queue):
            # NOTE:
            # This is a very dirty way to skip entries before they are fed
            # to the batch processes. This is needed for e.g. an empty line
            # or because the File Module noticed wrong formatting. These lines
            # used to be discarded but are now preserved by the escape string.
            # The benefit of this is that the users input will match the
            # output in terms of input line and outputline.
            if inputl.startswith("~!"): #Dirty Escape
                inputl = inputl[2:]
                if inputl:
                    flag = "S0"     # Flag for wrong format
                else:
                    flag = "S9"     # Flag for empty line
                    inputl = " " #Database doesn't like an empty inputfield
            else:
                flag = None
            if (i + 1) % columns:
                # Add flag for continuing the current row
                flag = '%s%s' % (flag if flag else '', 'C0')
            self.__database.addToQueue(jobID, inputl, flag)

        # Spawn child
        # Todo: Executable should be in bin/ directory.
        #p = subprocess.Popen(["MutalyzerBatch",
        #    "bin/batch_daemon"], executable="python")

        #Wait for the BatchChecker to fork of the Daemon
        #p.communicate()
        return jobID
    #addJob
#Scheduler
