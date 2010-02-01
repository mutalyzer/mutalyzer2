#!/usr/bin/python

import Config
import Db

def main(chrom, position) :
    C = Config.Config()
    D = Db.Db(C)
    
    print D.get_Transcripts(chrom, position)
    
    del D
    del C
#main

if __name__ == "__main__" :
    main("chr1", 159272155)
#if
