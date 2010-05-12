#!/usr/bin/python

from Services.Mapper import mappingObj, transcriptObj

from soaplib.wsgi_soap import SimpleWSGISoapApp
from soaplib.service import soapmethod
from soaplib.serializers.primitive import String, Integer, Array

class MutalyzerService(SimpleWSGISoapApp) :
    """
        Mutalyzer webservices.
    
        These methods are made public via a SOAP interface.
    
        Public methods:
          getTranscripts(chrom, pos)  ; Get all transcripts that overlap
            with a chromosomal position.
          getTranscriptsRange(chrom, posFirst, posLast, method) ; Get all
            transcripts that overlap with a range on a chromosome.
          getGeneName(transcript)  ; Find the gene name associated with a transcript.
          mappingInfo(LOVD_ver, build, accNo, variant) ; 
          transcriptInfo(LOVD_ver, build, accNo) ; Find transcription start and
            end, and CDS end (in c. notation) for a given transcript 
          extractChange(self, variant) : Extracts
          cTogConversion(self, LOVD_ver, build, variant) ; Convert c. to g.
          gTocConversion(self, LOVD_ver, build, variant) ; Convert g. to c.
            
    """
    @soapmethod(String, Integer, _returns = Array(String))
    def getTranscripts(self, chrom, pos) :
        """
            Get all the transcripts that overlap with a chromosomal position.
    
            Arguments:
                string chrom ; A chromosome encoded as "chr1", ..., "chrY".
                int    pos ; A position on the chromosome.
    
            Returns:
                string ; A list of transcripts.
        """
    
        from Modules import Config
        from Modules import Db
        from Modules import Output
    
        Conf = Config.Config()
        Database = Db.Db(Conf, "local")
        O = Output.Output(Conf, __file__)
    
        O.LogMsg(__file__, "Reveived request getTranscripts(%s %s)" % (chrom, pos))
        ret = str(Database.get_Transcripts(chrom, pos, pos, True))
        O.LogMsg(__file__, "Finished processing getTranscripts(%s %s)" % (
                 chrom, pos))
    
        del O, Database, Conf
        return ret
    #getTranscripts
    
    @soapmethod(String, Integer, Integer, Integer, _returns = Array(String))
    def getTranscriptsRange(self, chrom, posFirst, posLast, method) :
        """
            Get all the transcripts that overlap with a range on a chromosome.
    
            Arguments:
                string chrom ; A chromosome encoded as "chr1", ..., "chrY".
                int    posFirst ; The first position of the range.
                int    posLast ; The last position of the range.
                int    method ; The method of determining overlap:
                            0 ; Return all hit transcripts.
                            1 ; Return only the transcripts that completely fall
                                in the range [posFirst, posLast].
    
            Returns:
                string ; A list of transcripts.
        """
    
        from Modules import Config
        from Modules import Db
        from Modules import Output
    
        Conf = Config.Config()
        Database = Db.Db(Conf, "local")
        O = Output.Output(Conf, __file__)
    
        O.LogMsg(__file__, 
            "Reveived request getTranscriptsRange(%s %s %s %s)" % (
            chrom, posFirst, posLast, method))
        ret = str(Database.get_Transcripts(chrom, posFirst, posLast, method))
        O.LogMsg(__file__, 
            "Finished processing getTranscriptsRange(%s %s %s %s)" % (
            chrom, posFirst, posLast, method))
    
        del O, Database, Conf
        return ret
    #getTranscriptsRange
    
    @soapmethod(String, _returns = String)
    def getGeneName(self, transcript) :
        """
            Find the gene name associated with a transcript.
    
            Arguments:
                string transcript ; The identifier of a transcript.
    
            Returns:
                string ; The name of the associated gene.
        """
    
        from Modules import Config
        from Modules import Db
        from Modules import Output
    
        Conf = Config.Config()
        Database = Db.Db(Conf, "local")
        O = Output.Output(Conf, __file__)
    
        O.LogMsg(__file__, "Received request getGeneName(%s)" % transcript)
        ret = str(Database.get_GeneName(transcript.split('.')[0]))
        O.LogMsg(__file__, "Finished processing getGeneName(%s)" % transcript)
    
        del O, Database, Conf
        return ret
    #getGeneName

    @soapmethod(String, String, String, String, _returns = mappingObj)
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

        import Services.Mapper
        from Modules import Config
        from Modules import Output
    
        Conf = Config.Config()
        O = Output.Output(Conf, __file__)
    
        O.LogMsg(__file__, "Reveived request mappingInfo(%s %s %s %s)" % (
                 LOVD_ver, build, accNo, variant))

        result = Services.Mapper.mainMapping(LOVD_ver, build, accNo, variant)

        O.LogMsg(__file__, "Finished processing mappingInfo(%s %s %s %s)" % (
                 LOVD_ver, build, accNo, variant))

        del O, Conf
        return result
    #mappingInfo

    @soapmethod(String, String, String, _returns = transcriptObj)
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

        import Services.Mapper
        from Modules import Config
        from Modules import Output
    
        Conf = Config.Config()
        O = Output.Output(Conf, __file__)
    
        O.LogMsg(__file__, "Received request transcriptInfo(%s %s %s)" % (
                 LOVD_ver, build, accNo))

        result = Services.Mapper.mainTranscript(LOVD_ver, build, accNo)

        O.LogMsg(__file__, "Finished processing transcriptInfo(%s %s %s)" % (
                 LOVD_ver, build, accNo))

        del O, Conf
        return result
    #transcriptInfo
    
    @soapmethod(String, _returns = String)
    def extractChange(self, variant) :
        '''
            Extracts.

            
            Argument:
                string variant ; The variant in complete HGVS notation
             
            Returns:
                tuple: The position and the part behind the positions
                    in HGVS notation
        '''
        
        from Modules import Parser
        from Modules import Config
        from Modules import Output
    
        Conf = Config.Config() # Read the configuration file.
        O = Output.Output(Conf, __file__)
    
        O.LogMsg(__file__, "Received request extractChange(%s)" % (variant))
        

        # Make a parsetree for the given variation.
        P = Parser.Nomenclatureparser(O)
        parsetree = P.parse(variant) # 
        del P
        position = parsetree.RawVar.StartLoc.PtLoc.Main# start position variant
        
        if parsetree.RawVar.Arg1 :
            change = parsetree.RawVar.Arg1
            if parsetree.RawVar.Arg2 :
                change += '>' + parsetree.RawVar.Arg2
        #if
        else :
            change = ''
        #else
        if parsetree.RawVar.MutationType not in ('subst', 'del') :
            changeSuffix = parsetree.RawVar.MutationType + change
        #if
        else :
            if parsetree.RawVar.MutationType == 'del' :
                changeSuffix = parsetree.RawVar.MutationType
            else :
                changeSuffix = change
        #else
        del parsetree, change

        O.LogMsg(__file__, "Finished processing extractChange(%s)" % (variant))
        del O

        return position, changeSuffix
    #extractChange

    
    @soapmethod(String, String, String, _returns = String)
    def cTogConversion(self, LOVD_ver, build, variant) :
        """
            Converts a complete HGVS c. notation into a g. notation

            
            Arguments (all strings:
                LOVD_ver ; The LOVD version.
                build  ; The human genome build (ignored for now, hg19 assumed).
                variant ; The variant in c.notation.
             
            Returns:
                string var_in_g ; The variant in g. notation.
        
        """
        from Modules import Config
        from Modules import Output
        from Services import Mapper
    
        Conf = Config.Config() # Read the configuration file.
        O = Output.Output(Conf, __file__)
        del Conf
    
        O.LogMsg(__file__, "Received request cTogConversion(%s %s %s)" % (
                 LOVD_ver, build, variant))
        
        accno = variant.split(':')[0] # the transcript (NM_..) accession number 
        var = variant.split(':')[1]  # the variant
        changeSuffix = self.extractChange(variant)[1]
        mapped_mrnaAcc = Mapper.mainMapping(LOVD_ver, build, accno, var)
        if str(mapped_mrnaAcc.start_g) != str(mapped_mrnaAcc.end_g) :
            var_in_g = "g." + str(mapped_mrnaAcc.start_g) \
                        + "_" + str(mapped_mrnaAcc.end_g) + changeSuffix
        #if
        else :
            var_in_g = "g." + str(mapped_mrnaAcc.start_g) + changeSuffix
        #else
        
        del changeSuffix, mapped_mrnaAcc

        O.LogMsg(__file__, "Finished processing cTogConversion(%s %s %s)" % (
                 LOVD_ver, build, variant))
        del O

        return var_in_g
    #cTogConversion
        
    @soapmethod(String, String, String, _returns = String)
    def gTocConversion(self, LOVD_ver, build, variant) :
        """
            Converts a complete HGVS g. notation into a c. notation

            
            Arguments (all strings):
                LOVD_ver ; The LOVD version.
                build  ; The human genome build (ignored for now, hg19 assumed).
                variant ; The variant in g. notation.
             
            Returns:
                string var_in_g ; The variant in c. notation.
        
        
        """
        from Modules import Config
        from Modules import Output
        from Modules import Parser
        from Modules import Db
        from Services import Mapper
    
        Conf = Config.Config() # Read the configuration file.
        Database = Db.Db(Conf, "local") # Open the database.
        O = Output.Output(Conf, __file__)
    
        O.LogMsg(__file__, "Received request gTocConversion(%s %s %s)" % (
                 LOVD_ver, build, variant))
        

        var = variant.split(':')[1]  # the variant
        pos1 = self.extractChange(variant)[0]
        changeSuffix = self.extractChange(variant)[1]
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
#            info_mrnaAcc = self.transcriptInfo("123", "hg19", mrnaAcc)
#            mapped_mrnaAcc = self.mappingInfo("123", "hg19", mrnaAcc, variant)
            info_mrnaAcc = Mapper.mainTranscript("123", "hg19", mrnaAcc)
            mapped_mrnaAcc = Mapper.mainMapping("123", "hg19", mrnaAcc, var)

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

        O.LogMsg(__file__, "Finished processing gTocConversion(%s %s %s)" % (
                 LOVD_ver, build, variant))
        del O
    #gTocConversion

#MutalyzerService
