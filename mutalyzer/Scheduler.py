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

from __future__ import unicode_literals

import io
import os                               # os.path.exists
import smtplib                          # smtplib.STMP
from email.mime.text import MIMEText    # MIMEText
from sqlalchemy import func
from sqlalchemy.orm.exc import NoResultFound

from mutalyzer.config import settings
from mutalyzer.db import queries, session
from mutalyzer.db.models import Assembly, BatchJob, BatchQueueItem
from mutalyzer import ncbi
from mutalyzer import stats
from mutalyzer import variantchecker
from mutalyzer.grammar import Grammar
from mutalyzer.output import Output
from mutalyzer.mapping import Converter
from mutalyzer import website


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

    def __init__(self) :
        #TODO: documentation
        """
        Initialize the Scheduler, which requires a database connection.

        @todo: documentation
        """
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

    def __sendMail(self, mailTo, result_id):
        """
        Send an e-mail containing an url to a batch job submitter.

        @todo: Handle Connection errors in a try, except clause

        @arg mailTo: The batch job submitter
        @type mailTo: unicode
        @arg result_id: Identifier for the job result.
        @type result_id: unicode
        """
        if settings.TESTING:
            return

        # Mail is set to '<IP ADDRESS>@<INTERFACE>.mutalyzer' if the batch job
        # was submitted without specifying an email address.
        if mailTo.endswith('.mutalyzer'):
            return

        #TODO: Handle Connection errors in a try, except clause
        #Expected errors: socket.error

        def encode_address(address):
            """
            Unfortunately, smtplib can only handle ASCII email addresses.

            The domain name part should be punycode-encoded. The part before
            the '@' should be encoded as ASCII if possible or as UTF8 if the
            first-hop mail server supports this.

            This is way to complicated, so we just encode the entire address
            as ASCII or punycode if that fails.

            https://bugs.python.org/issue20084
            https://bugs.python.org/issue20083
            """
            try:
                return address.encode('ascii')
            except UnicodeEncodeError:
                return address.encode('idna')

        from_address = 'Mutalyzer batch job <%s>' % (
            settings.BATCH_NOTIFICATION_EMAIL or settings.EMAIL)
        download_url = website.url_for('batch_job_result',
                                       result_id=result_id)

        message = MIMEText("""Dear submitter,

The batch operation you have submitted, has been processed successfully.

Your results can be found here:
%s

Thanks for using Mutalyzer.


With kind regards,
Mutalyzer batch scheduler""" % download_url)

        message["Subject"] = "Result of your Mutalyzer batch job"
        message["From"] = from_address
        message["To"] = encode_address(mailTo)

        try:
            smtpInstance = smtplib.SMTP('localhost')
        except smtplib.SMTPConnectError as e:
            print 'Could not send email to %s: (%s) %s' % (mailTo,
                                                           unicode(e.smtp_code),
                                                           unicode(e.smtp_error))
            return
        except IOError as e:
            print 'Could not send email to %s: %s' % (mailTo, unicode(e))
            return

        try:
            smtpInstance.sendmail(from_address, encode_address(mailTo),
                                  message.as_string())
        except smtplib.SMTPRecipientsRefused as e:
            for address, (code, error) in e.recipients.items():
                print 'Could not send email to %s: (%s) %s' % (address,
                                                               code,
                                                               error)
        except smtplib.SMTPResponseException as e:
            print 'Could not send email to %s: (%s) %s' % (mailTo,
                                                           unicode(e.smtp_code),
                                                           unicode(e.smtp_error))
        finally:
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
            O.addMessage(__file__, 3, "EBSKIP", message)
            return True #skip
        #if
        if 'A' in flags : #This entry is altered before execution
            O.addMessage(__file__, 2, "WEALTER", "Entry altered before "
                    "execution")
        return False
    #__processFlags

    def __alterBatchEntries(self, jobID, old, new, flag, nselector, O) :
        """
        Replace within one JobID all entries matching old with new, if they do
        not match the negative selector.

        This is used to alter batch entries that would otherwise take a long
        time to process. E.g. a batch job with a lot of the same accession
        numbers without version numbers would take a long time because
        mutalyzer would fetch the file from the NCBI for each entry. A
        database update over all entries with the same accession number speeds
        up the job considerably.

        Example:
        NM_002001(FCER1A_v001):c.1A>C ; this would result in the continuous
        fetching of the reference because no version number is given.
        In this case the arguments would be:
            - old         ;   NM_002001
            - new         ;   NM_002001.2
            - nselector   ;   NM_002001.

        The nselector is used to prevent the replacement of
        false positives. e.g. NM_002001.1(FCER1A_v001):c.1A>C should not
        be replaced. For this reason, any items starting with the nselector
        value are ignored.

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
        #query = '''UPDATE batch_queue_items
        #             SET item = REPLACE(item, :old, :new),
        #                 flags = flags || :flag
        #             WHERE batch_job_id = :batch_job_id
        #                   AND NOT item LIKE :nselector%'''
        #parameters = {'batch_job_id': jobID,
        #              'old': old,
        #              'new': new,
        #              'flag': flag,
        #              'nselector': nselector}
        #session.execute(query, parameters)
        try:
            BatchQueueItem.query \
                .filter_by(batch_job_id=jobID) \
                .filter(BatchQueueItem.item.startswith(old),
                        ~BatchQueueItem.item.startswith(nselector)) \
                .update({'item': func.replace(BatchQueueItem.item, old, new),
                         'flags': BatchQueueItem.flags + flag},
                            synchronize_session=False)
        except Exception as ex:
            message = ("An exception of type '%s' occurred in __alterBatchEntries() "
                       "with the following arguments: %s. "
                       "Other info: old=%s, new=%s, flag=%s, nselector=%s"
                       % (type(ex).__name__, ex.args, old, new, flag, nselector))
            O.addMessage(__file__, 4, "ABATCHE", message)
        session.commit()
    #__alterBatchEntries

    def __skipBatchEntries(self, jobID, flag, selector) :
        """
        Skip all batch entries that match a certain selector.

        We flag batch entries to be skipped. This is used if it is certain
        that an entry will cause an error, or that its output is ambiguous.

        @arg jobID:
        @type jobID:
        @arg flag:
        @type flag:
        @arg selector:
        @type selector:
        """
        #update `BatchQueue` set
        #  `Flags` = CONCAT(IFNULL(`Flags`, ""), %s)
        #  where `JobID` = %s AND
        #  `Input` RLIKE %s;
        BatchQueueItem.query \
            .filter_by(batch_job_id=jobID) \
            .filter(BatchQueueItem.item.startswith(selector)) \
            .update({'flags': BatchQueueItem.flags + flag},
                    synchronize_session=False)
        session.commit()
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
        # E.g. ("S1", "NM_002001.")
        # An altered entry has three arguments,
        #               old,           new          negative selector
        # E.g.("A2",("NM_002001", "NM_002001.2", "NM_002001."))

        # Flags are set when an entry could be sped up. This is either the
        # case for the Retriever as for the Mutalyzer module

        if not flags: return
        #First check if we need to skip
        for flag, args in flags :
            if 'S' in flag :
                selector = args     # Strip argument
                O.addMessage(__file__, 2, "WBSKIP",
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
                O.addMessage(__file__, 2, "WBSUBST",
                        "All further occurrences of %s will be substituted "
                        "by %s" % (old, new))
                self.__alterBatchEntries(jobID, old, new, flag, nselector, O)
            #if
        #for
    #_updateDbFlags

    def process(self):
        """
        Start the mutalyzer Batch Processing. This method retrieves all jobs
        jobs from the database and processes them in a roundrobin fashion.
        After each round, the process checks if new jobs are added during the
        last processing round and repeats. This continue until no jobs are
        left to process.

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
        while not self.stopped():
            # Group batch jobs by email address and retrieve the oldest for
            # each address. This improves fairness when certain users have
            # many jobs.
            batch_jobs = BatchJob.query.filter(BatchJob.id.in_(
                session.query(func.min(BatchJob.id)).group_by(BatchJob.email))
            ).all()

            if len(batch_jobs) == 0:
                break

            for batch_job in batch_jobs:
                if self.stopped():
                    break

                batch_queue_item = queries.pop_batch_queue_item(batch_job)

                if batch_queue_item is not None:
                    item, flags = batch_queue_item

                    if batch_job.job_type == 'name-checker':
                        self._processNameBatch(batch_job, item, flags)
                    elif batch_job.job_type == 'syntax-checker':
                        self._processSyntaxCheck(batch_job, item, flags)
                    elif batch_job.job_type == 'position-converter':
                        self._processConversion(batch_job, item, flags)
                    elif batch_job.job_type == 'snp-converter':
                        self._processSNP(batch_job, item, flags)
                    else:
                        # Unknown job type, should never happen.
                        # Todo: Log some screaming message.
                        pass

                else:
                    print ('Job %s finished, email %s file %s' %
                           (batch_job.id, batch_job.email, batch_job.result_id))
                    self.__sendMail(batch_job.email, batch_job.result_id)
                    session.delete(batch_job)
                    session.commit()
    #process

    def _processNameBatch(self, batch_job, cmd, flags):
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

        stats.increment_counter('name-checker/batch')

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
                O.addMessage(__file__, -1, "DEBUG", unicode(repr(traceback.format_exc())))
                # We don't know what caused the exception, but it might be
                # necessary to do a rollback on the session transaction in
                # order to continue using the session (which is important in
                # a batch setup).
                # TODO: It makes sense to do this not only for name checker
                # jobs, but on a higher level in the batch processing stack.
                # The reason for doing it here is (I think) so there's still
                # a batch result entry written.
                session.rollback()
            #except
            finally :
                #check if we need to update the database
                self._updateDbFlags(O, batch_job.id)
        #if

        batchOutput = O.getOutput("batchDone")

        outputline =  "%s\t" % cmd
        outputline += "%s\t" % "|".join(O.getBatchMessages(2))

        if batchOutput :
            outputline += batchOutput[0]

        #Output
        filename = "%s/batch-job-%s.txt" % (settings.CACHE_DIR, batch_job.result_id)
        if not os.path.exists(filename) :
            # If the file does not yet exist, create it with the correct
            # header above it. The header is read from the config file as
            # a list. We need a tab delimited string.
            header = ['Input',
                      'Errors and warnings',
                      'AccNo',
                      'Genesymbol',
                      'Variant',
                      'Reference Sequence Start Descr.',
                      'Coding DNA Descr.',
                      'Protein Descr.',
                      'GeneSymbol Coding DNA Descr.',
                      'GeneSymbol Protein Descr.',
                      'Genomic Reference',
                      'Coding Reference',
                      'Protein Reference',
                      'Affected Transcripts',
                      'Affected Proteins',
                      'Restriction Sites Created',
                      'Restriction Sites Deleted']
            handle = io.open(filename, mode='a', encoding='utf-8')
            handle.write("%s\n" % "\t".join(header))
        #if
        else :
            handle = io.open(filename, mode='a', encoding='utf-8')

        if flags and 'C' in flags:
            separator = '\t'
        else:
            separator = '\n'

        handle.write("%s%s" % (outputline, separator))
        handle.close()
        O.addMessage(__file__, -1, "INFO",
            "Finished NameChecker batchvariant " + cmd)
    #_processNameBatch

    def _processSyntaxCheck(self, batch_job, cmd, flags):
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

        stats.increment_counter('syntax-checker/batch')

        skip = self.__processFlags(output, flags)
        #Process
        if not skip :
            parsetree = grammar.parse(cmd)
        else :
            parsetree = None

        if parsetree :
            result = "OK"
        else :
            result = "|".join(output.getBatchMessages(2))

        #Output
        filename = "%s/batch-job-%s.txt" % (settings.CACHE_DIR, batch_job.result_id)
        if not os.path.exists(filename) :
            # If the file does not yet exist, create it with the correct
            # header above it. The header is read from the config file as
            # a list. We need a tab delimited string.
            header = ['Input', 'Status']
            handle = io.open(filename, mode='a', encoding='utf-8')
            handle.write("%s\n" % "\t".join(header))
        #if
        else :
            handle = io.open(filename, mode='a', encoding='utf-8')

        if flags and 'C' in flags:
            separator = '\t'
        else:
            separator = '\n'

        handle.write("%s\t%s%s" % (cmd, result, separator))
        handle.close()
        output.addMessage(__file__, -1, "INFO",
                          "Finished SyntaxChecker batchvariant " + cmd)
    #_processSyntaxCheck

    def _processConversion(self, batch_job, cmd, flags):
        """
        Process an entry from the Position Converter, write the results
        to the job-file. The Position Converter is wrapped in a try except
        block which ensures that he Batch Process keeps running. Errors
        are caught and the user will be notified.

        Side-effect:
            - Output written to outputfile.

        @arg cmd: The Syntax Checker input
        @type cmd: unicode
        @arg i: The JobID
        @type i: integer
        @arg build: The build to use for the converter
        @type build: unicode
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

        stats.increment_counter('position-converter/batch')

        skip = self.__processFlags(O, flags)
        if not skip :
            try :
                #process
                try:
                    assembly = Assembly.by_name_or_alias(batch_job.argument)
                except NoResultFound:
                    O.addMessage(__file__, 3, 'ENOASSEMBLY',
                                 'Not a valid assembly: ' + batch_job.argument)
                    raise

                converter = Converter(assembly, O)

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

        error = "%s" % "|".join(O.getBatchMessages(2))

        #Output
        filename = "%s/batch-job-%s.txt" % (settings.CACHE_DIR, batch_job.result_id)
        if not os.path.exists(filename) :
            # If the file does not yet exist, create it with the correct
            # header above it. The header is read from the config file as
            # a list. We need a tab delimited string.
            header = ['Input Variant',
                      'Errors',
                      'Chromosomal Variant',
                      'Coding Variant(s)']
            handle = io.open(filename, mode='a', encoding='utf-8')
            handle.write("%s\n" % "\t".join(header))
        #if
        else :
            handle = io.open(filename, mode='a', encoding='utf-8')

        if flags and 'C' in flags:
            separator = '\t'
        else:
            separator = '\n'

        handle.write("%s\t%s\t%s\t%s%s" % (cmd, error, gName, "\t".join(cNames), separator))
        handle.close()
        O.addMessage(__file__, -1, "INFO",
            "Finisehd PositionConverter batchvariant " + cmd)
    #_processConversion


    def _processSNP(self, batch_job, cmd, flags):
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

        stats.increment_counter('snp-converter/batch')

        #Read out the flags
        # Todo: Do something with the flags?
        skip = self.__processFlags(O, flags)

        descriptions = []
        if not skip:
            descriptions = ncbi.rsid_to_descriptions(cmd, O)

        # Todo: Is output ok?
        outputline =  "%s\t" % cmd
        outputline += "%s\t" % "|".join(descriptions)
        outputline += "%s\t" % "|".join(O.getBatchMessages(2))

        #Output
        filename = "%s/batch-job-%s.txt" % (settings.CACHE_DIR, batch_job.result_id)
        if not os.path.exists(filename) :
            # If the file does not yet exist, create it with the correct
            # header above it. The header is read from the config file as
            # a list. We need a tab delimited string.
            header = ['Input Variant',
                      'HGVS description(s)',
                      'Errors and warnings']
            handle = io.open(filename, mode='a', encoding='utf-8')
            handle.write("%s\n" % "\t".join(header))
        #if
        else :
            handle = io.open(filename, mode='a', encoding='utf-8')

        if flags and 'C' in flags:
            separator = '\t'
        else:
            separator = '\n'

        handle.write("%s%s" % (outputline, separator))
        handle.close()
        O.addMessage(__file__, -1, "INFO",
                     "Finished SNP converter batch rs%s" % cmd)
    #_processSNP

    def addJob(self, email, queue, columns, job_type, argument=None):
        """
        Add a job to the Database and start the BatchChecker.

        @arg email:         e-mail address of batch supplier
        @type email:        unicode
        @arg queue:         A list of jobs
        @type queue:        list
        @arg columns:       The number of columns.
        @type columns:      int
        @arg job_type:       The type of Batch Job that should be run
        @type job_type:
        @arg argument:          Batch Arguments, for now only build info
        @type argument:

        @return: result_id
        @rtype:
        """
        # Add jobs to the database
        batch_job = BatchJob(job_type, email=email, argument=argument)
        session.add(batch_job)

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

            item = BatchQueueItem(batch_job, inputl, flags=flag)
            session.add(item)

        session.commit()
        return batch_job.result_id
    #addJob
#Scheduler
