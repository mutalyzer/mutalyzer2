#!/usr/bin/python

import Config
import Db

def main(name) :
    C = Config.Config()
    D = Db.Db(C)
    
    print D.get_GeneName(name.split('.')[0]), # No idea why there must be a ,
    
    del D
    del C
#main

if __name__ == "__main__" :
    main("NM_002001.2")
#if
