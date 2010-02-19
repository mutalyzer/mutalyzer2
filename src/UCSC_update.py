#!/usr/bin/python

"""
    Get updates on mapping information from the UCSC.

    This program is intended to be run dayly from cron. 
"""

import sys
import os
os.chdir(sys.argv[0].rsplit('/', 2)[0])

from Modules import Config
from Modules import Db
from Modules import Output

C = Config.Config()
O = Output.Output(C, __file__)
O.LogMsg(__file__, "Starting UCSC mapping data update")

RemoteDb = Db.Db(C, "remote")
LocalDb = Db.Db(C, "local")

RemoteDb.get_Update()
LocalDb.load_Update()

count_Updates = LocalDb.count_Updates()
if count_Updates :
    O.LogMsg(__file__, "%i updates found" % count_Updates)
    LocalDb.backup_cdsUpdates()
    cds_Updates = LocalDb.count_cdsUpdates()
    if cds_Updates :
        O.LogMsg(__file__, "%i CDS updates found, backing up" % cds_Updates)
    LocalDb.merge_cdsUpdates()
#if
LocalDb.merge_Update()

O.LogMsg(__file__, "UCSC mapping data update end")

del LocalDb
del RemoteDb
del O
del C
