#!/usr/bin/python


from soaplib.wsgi_soap import SimpleWSGISoapApp
from soaplib.service import soapmethod
from soaplib.serializers.primitive import String, Integer, Array
from soaplib.serializers.clazz import ClassSerializer
from ZSI import TC
from ZSI.fault import Fault

from Modules import Web
from Modules import Db
from Modules import Output
from Modules import Config
from Modules import Parser
from Modules import Mapper
    
class MutalyzerService(SimpleWSGISoapApp) :
    """
        Mutalyzer webservices.
    
        These methods are made public via a SOAP interface.
    
        Private methods:
            __checkBuild(L, D, build) ; Check if the build is supported.
            __checkChrom(L, D, chrom) ; Check if the chromosome is in our 
                                        database.
            __checkPos(L, pos)        ; Check if the posision is valid.

        Public methods:
            getTranscripts(build, chrom, pos)    ; Get all transcripts that 
                                                   overlap with a chromosomal 
                                                   position.
            getTranscriptsRange(build, chrom,    ; Get all transcripts that 
                                pos1, pos2,        overlap with a range on a 
                                method)            chromosome.                 
            getGeneName(build, accno)            ; Find the gene name associated
                                                   with a transcript.
            mappingInfo(LOVD_ver, build, accNo,  ; FIXME
                        variant) ; 
            transcriptInfo(LOVD_ver, build, ; Find transcription start and
                           accNo)             end, and CDS end (in c. 
                                              notation) for a given transcript 
            __extractChange(self, variant, O)  ; Extracts the part of a variant 
                                                 description after the
                                                 coordinates (positions)
            cTogConversion(self, build, variant) ; Convert c. to g.
            gTocConversion(self, build, variant) ; Convert g. to c.
    """

    def __checkBuild(self, L, D, build) :
        """
            Check if the build is supported (hg18 or hg19).

            Arguments:
                L     ; An output object for logging.
                D     ; A handle to the database.
                build ; The build name that needs to be checked.

            Returns:
                Nothing (but raises an EARG exception).
        """

        if not D.opened :
            L.addMessage(__file__, 4, "EARG", "EARG %s" % build)
            raise Fault(Fault.Client, "EARG", 
                detail = "The build argument (%s) was not a valid " \
                         "build name." % build)
        #if                         
    #__checkBuild

    def __checkChrom(self, L, D, chrom) :
        """
            Check if the chromosome is in our database.

            Arguments:
                L     ; An output object for logging.
                D     ; A handle to the database.
                chrom ; The name of the chromosome.

            Returns:
                Nothing (but raises an EARG exception).
        """

        if not D.isChrom(chrom) :
            L.addMessage(__file__, 4, "EARG", "EARG %s" % chrom)
            raise Fault(Fault.Client, "EARG", 
                detail = "The chrom argument (%s) was not a valid " \
                         "chromosome name." % chrom)
        #if                         
    #__checkChrom

    def __checkPos(self, L, pos) :
        """
            Check if the position is valid.

            Arguments:
                L   ; An output object for logging.
                pos ; The position.

            Returns:
                Nothing (but raises an ERANGE exception).
        """

        if pos < 1 :
            L.addMessage(__file__, 4, "ERANGE", "ERANGE %i" % pos)
            raise Fault(Fault.Client, "ERANGE", 
                detail = "The pos argument (%i) is out of range." % pos)
        #if                         
    #__checkPos

    #@soapmethod(String, _returns = String)
    def __extractChange(self, variant, O) :
        '''
            Extracts.

            
            Argument:
                string variant ; The variant in complete HGVS notation
             
            Returns:
                tuple: The position and the part behind the positions
                    in HGVS notation
        '''
        
        P = Parser.Nomenclatureparser(O)
        parsetree = P.parse(variant) # 
        del P

        position = parsetree.RawVar.StartLoc.PtLoc.Main # Start position.
        
        if parsetree.RawVar.Arg1 :
            change = parsetree.RawVar.Arg1
            if parsetree.RawVar.Arg2 :
                change += '>' + parsetree.RawVar.Arg2
        #if
        else :
            change = ''

        if parsetree.RawVar.MutationType not in ('subst', 'del') :
            changeSuffix = parsetree.RawVar.MutationType + change
        else :
            if parsetree.RawVar.MutationType == 'del' :
                changeSuffix = parsetree.RawVar.MutationType
            else :
                changeSuffix = change
        #else

        del parsetree, change

        return position, changeSuffix
    #__extractChange
    
    @soapmethod(String, String, Integer, _returns = Array(String))
    def getTranscripts(self, build, chrom, pos) :
        """
            Get all the transcripts that overlap with a chromosomal position.
    
            Arguments:
                string build ; The build name encoded as "hg18" or "hg19".
                string chrom ; A chromosome encoded as "chr1", ..., "chrY".
                int    pos   ; A postion on the chromosome.
    
            Returns:
                string ; A list of transcripts.

            On error an exception is raised:
                detail       ; Human readable description of the error.
                faultstring: ; A code to indicate the type of error.
                    EARG   ; The argument was not valid.
                    ERANGE ; An invalid range was given.
        """
    
        C = Config.Config()
        L = Output.Output(__file__, C.Output)
    
        L.addMessage(__file__, -1, "INFO", 
                     "Received request getTranscripts(%s %s %s)" % (build,
                     chrom, pos))

        D = Db.Db("local", build, C.Db)
        self.__checkBuild(L, D, build)

        self.__checkChrom(L, D, chrom)
        self.__checkPos(L, pos)
        
        ret = D.get_Transcripts(chrom, pos, pos, True)
        L.addMessage(__file__, -1, "INFO", 
                     "Finished processing getTranscripts(%s %s %s)" % (build,
                     chrom, pos))
    
        del D, L, C
        return [ret]
    #getTranscripts
    
    @soapmethod(String, String, Integer, Integer, Integer, 
                _returns = Array(String))
    def getTranscriptsRange(self, build, chrom, pos1, pos2, method) :
        """
            Get all the transcripts that overlap with a range on a chromosome.
    
            Arguments:
                string build ; The build name encoded as "hg18" or "hg19".
                string chrom  ; A chromosome encoded as "chr1", ..., "chrY".
                int    pos1   ; The first postion of the range.
                int    pos2   ; The last postion of the range.
                int    method ; The method of determining overlap:
                                0 ; Return only the transcripts that completely
                                    fall in the range [pos1, pos2].
                                1 ; Return all hit transcripts.
    
            Returns:
                string ; A list of transcripts.
        """
    
        C = Config.Config()
        L = Output.Output(__file__, C.Output)
    
        L.addMessage(__file__, -1, "INFO", 
            "Received request getTranscriptsRange(%s %s %s %s %s)" % (build,
            chrom, pos1, pos2, method))

        D = Db.Db("local", build, C.Db)
        self.__checkBuild(L, D, build)

        ret = D.get_Transcripts(chrom, pos1, pos2, method)
        L.addMessage(__file__, -1, "INFO", 
            "Finished processing getTranscriptsRange(%s %s %s %s %s)" % (
            build, chrom, pos1, pos2, method))
    
        del D
        del L
        del C
        return [ret]
    #getTranscriptsRange

    @soapmethod(String, String, _returns = String)
    def getGeneName(self, build, accno) :
        """
            Find the gene name associated with a transcript.
    
            Arguments:
                string build ; The build name encoded as "hg18" or "hg19".
                string accno ; The identifier of a transcript.
    
            Returns:
                string ; The name of the associated gene.
        """
    
        C = Config.Config()
        L = Output.Output(__file__, C.Output)
        L.addMessage(__file__, -1, "INFO", 
                     "Received request getGeneName(%s %s)" % (build, accno))

        D = Db.Db("local", build, C.Db)
        self.__checkBuild(L, D, build)
    
        ret = D.get_GeneName(accno.split('.')[0])
        L.addMessage(__file__, -1, "INFO", 
                     "Finished processing getGeneName(%s %s)" % (build, accno))
    
        del D
        del L
        del C
        return ret
    #getGeneName

    @soapmethod(String, String, String, String, _returns = Mapper.Mapping)
    def mappingInfo(self, LOVD_ver, build, accNo, variant) :
        """
            Search for an NM number in the MySQL database, if the version
            number matches, get the start and end positions in a variant and
            translate these positions to g. notation if the variant is in c.
            notation and vice versa.

            - If no end position is present, the start position is assumed to
              be the end position. 
            - If the version number is not found in the database, an error
              message is generated and a suggestion for an other version is
              given.
            - If the reference sequence is not found at all, an error is
              returned.
            - If no variant is present, an error is returned.
            - If the variant is not accepted by the nomenclature parser, a
              parse error will be printed.

            
            Arguments (all strings):
                LOVD_ver ; The LOVD version.
                build ; The human genome build (ignored for now, hg19 assumed).
                accNo ; The NM accession number and version.
                variant ; The variant.
             
            Returns:
                complex object:
                    start_main   ; The main coordinate of the start position 
                                   in c. (non-star) notation.
                    start_offset ; The offset coordinate of the start position
                                   in c. notation (intronic position).
                    end_main     ; The main coordinate of the end position in 
                                   c. (non-star) notation.
                    end_offset   ; The offset coordinate of the end position in
                                   c. notation (intronic position).
                    start_g      ; The g. notation of the start position.
                    end_g        ; The g. notation of the end position.
                    type         ; The mutation type.

        """

        Conf = Config.Config()
        O = Output.Output(__file__, Conf.Output)
    
        O.addMessage(__file__, -1, "INFO", 
                     "Reveived request mappingInfo(%s %s %s %s)" % (LOVD_ver,
                     build, accNo, variant))

        result = Mapper.mainMapping(LOVD_ver, build, accNo, variant,
                                             Conf, O)

        O.addMessage(__file__, -1, "INFO", 
                     "Finished processing mappingInfo(%s %s %s %s)" % (
                     LOVD_ver, build, accNo, variant))

        del O, Conf
        return result
    #mappingInfo

    @soapmethod(String, String, String, _returns = Mapper.Transcript)
    def transcriptInfo(self, LOVD_ver, build, accNo) :
        """
            Search for an NM number in the MySQL database, if the version
            number matches, the transcription start and end and CDS end 
            in c. notation is returned.

            
            Arguments (all strings:
                LOVD_ver ; The LOVD version.
                build ; The human genome build (ignored for now, hg19 assumed).
                accNo ; The NM accession number and version.
             
            Returns:
                complex object:
                    trans_start  ; Transcription start in c. notation.
                    trans_stop   ; Transcription stop in c. notation.
                    CDS_stop     ; CDS stop in c. notation.
        """

        Conf = Config.Config()
        O = Output.Output(__file__, Conf.Output)
    
        O.addMessage(__file__, -1, "INFO", 
                     "Received request transcriptInfo(%s %s %s)" % (LOVD_ver,
                     build, accNo))

        result = Mapper.mainTranscript(LOVD_ver, build, accNo, Conf, O)

        O.addMessage(__file__, -1, "INFO",
                     "Finished processing transcriptInfo(%s %s %s)" % (
                     LOVD_ver, build, accNo))

        del O, Conf
        return result
    #transcriptInfo
    
    @soapmethod(String, String, _returns = String)
    def cTogConversion(self, build, variant) :
        """
            Converts a complete HGVS c. notation into a g. notation

            
            Arguments (all strings:
                build  ; The human genome build (ignored for now, hg19 assumed).
                variant ; The variant in c.notation.
             
            Returns:
                string var_in_g ; The variant in g. notation.
        
        """
    
        Conf = Config.Config() # Read the configuration file.
        O = Output.Output(__file__, Conf.Output)
    
        O.addMessage(__file__, -1, "INFO",
                     "Received request cTogConversion(%s %s)" % (
                     build, variant))
        
        accno = variant.split(':')[0] # the transcript (NM_..) accession number 
        var = variant.split(':')[1]  # the variant
        changeSuffix = self.__extractChange(variant, O)[1]
        mapped_mrnaAcc = Mapper.mainMapping("123", build, accno, var, Conf,
                                            O)
        if str(mapped_mrnaAcc.start_g) != str(mapped_mrnaAcc.end_g) :
            var_in_g = "g." + str(mapped_mrnaAcc.start_g) \
                        + "_" + str(mapped_mrnaAcc.end_g) + changeSuffix
        #if
        else :
            var_in_g = "g." + str(mapped_mrnaAcc.start_g) + changeSuffix
        #else
        
        del changeSuffix, mapped_mrnaAcc

        O.addMessage(__file__, -1, "INFO",
                     "Finished processing cTogConversion(%s %s)" % (
                     build, variant))
        del O

        return var_in_g
    #cTogConversion
        
    @soapmethod(String, String, _returns = String)
    def gTocConversion(self, build, variant) :
        """
            Converts a complete HGVS g. notation into a c. notation

            
            Arguments (all strings):
                build  ; The human genome build (ignored for now, hg19 assumed).
                variant ; The variant in g. notation.
             
            Returns:
                string var_in_g ; The variant in c. notation.
        
        
        """
    
        Conf = Config.Config() # Read the configuration file.
        Database = Db.Db("local", build, Conf.Db)

        O = Output.Output(__file__, Conf.Output)
    
        O.addMessage(__file__, -1, "INFO",
                     "Received request gTocConversion(%s %s)" % (
                     build, variant))

        var = variant.split(':')[1]  # the variant
        pos1 = self.__extractChange(variant, O)[0]
        changeSuffix = self.__extractChange(variant, O)[1]
        upPos = int(pos1) - 5000
        downPos = int(pos1) + 5000
        chrom = variant.split(':')[0]    # the chromosome (for now: e.g. chr1)
        transcripts = eval(self.getTranscriptsRange('chrX', upPos, downPos, 1))
        HGVS_notations = '' # To contain >= 1 HGVS notations (transcripts)
        for mrnaAcc in transcripts:
            if '.' in mrnaAcc :
                accno = mrnaAcc.split('.')[0] # The NM accession number.

            # Retrieve info on the NM accession number.
            result = Database.get_NM_info(accno)
            strand = result[4] # The orientation.
            del result

            info_mrnaAcc = Mapper.mainTranscript("123", "hg19", mrnaAcc, Conf,
                                                 O)
            mapped_mrnaAcc = Mapper.mainMapping("123", "hg19", mrnaAcc, var, 
                                                Conf, O)

            vStart = Mapper.conversionToCoding(str(mapped_mrnaAcc.startoffset), 
                                          mapped_mrnaAcc.startmain, 
                                          info_mrnaAcc.trans_start, 
                                          info_mrnaAcc.trans_stop, 
                                          info_mrnaAcc.CDS_stop)
            # end position of the variant
            vEnd = Mapper.conversionToCoding(str(mapped_mrnaAcc.endoffset), 
                                        mapped_mrnaAcc.endmain, 
                                        info_mrnaAcc.trans_start, 
                                        info_mrnaAcc.trans_stop, 
                                        info_mrnaAcc.CDS_stop)
    
            # Use n. or c. numbering
            if info_mrnaAcc.trans_start == "1" and \
                info_mrnaAcc.trans_stop == info_mrnaAcc.CDS_stop :
                numbType = ":n."
            #if
            else :
                numbType = ":c."
            #else
    
            if strand == '+' :
                var_in_c = mrnaAcc + numbType + str(vStart[1]) + str(vStart[0])
                if (vStart[0] and vEnd[0] and vStart[0] == vEnd[0]) \
                    or mapped_mrnaAcc.mutationType == 'subst' :
                    var_in_c += changeSuffix
                else :
                    var_in_c += "_" + str(vEnd[1]) + str(vEnd[0]) \
                                + changeSuffix
            else :
                # reversed orientation
                var_in_c = mrnaAcc + numbType + str(vEnd[1]) + str(vEnd[0])
                if vStart[0] and vEnd[0] and vStart[0] == vEnd[0] \
                or mapped_mrnaAcc.mutationType == 'subst' :
                    var_in_c += changeSuffix
                else :
                    var_in_c += "_" + str(vStart[1]) + str(vStart[0]) \
                                + changeSuffix
            #if
            del vStart, vEnd, info_mrnaAcc, mapped_mrnaAcc
            
            if HGVS_notations == '' :
                HGVS_notations = var_in_c
            else :
                HGVS_notations += ";" + var_in_c
        #for
        return HGVS_notations

        O.addMessage(__file__, -1, "INFO",
                     "Finished processing gTocConversion(%s %s)" % (
                     build, variant))
        del O
    #gTocConversion
#MutalyzerService
