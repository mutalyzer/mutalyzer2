#!/usr/bin/python

"""
    Get updates on mapping information from the UCSC.

    This program is intended to be run dayly from cron. 
"""

import sys
import os
os.chdir(sys.argv[0].rsplit('/', 2)[0])

from Modules import Config
from Modules import Output
from Modules import Db

C = Config.Config()
O = Output.Output(__file__, C.Output)
O.addMessage(__file__, -1, "INFO", "Starting UCSC mapping data update")

for i in C.Db.dbNames :
    RemoteDb = Db.Db("remote", i, C.Db)
    LocalDb = Db.Db("local", i, C.Db)
    
    RemoteDb.get_Update()
    LocalDb.load_Update()
    
    count_Updates = LocalDb.count_Updates()
    if count_Updates :
        O.addMessage(__file__, -1, "INFO", "%i updates found" % count_Updates)
        LocalDb.backup_cdsUpdates()
        cds_Updates = LocalDb.count_cdsUpdates()
        if cds_Updates :
            O.addMessage(__file__, -1, 
                         "%i CDS updates found, backing up" % cds_Updates)
        LocalDb.merge_cdsUpdates()
    #if
    LocalDb.merge_Update()
#for    

O.addMessage(__file__, -1, "INFO", "UCSC mapping data update end")

del LocalDb
del RemoteDb
del O
del C
