#!/usr/bin/python

"""
Get updates on mapping information from the UCSC.

This program is intended to be run daily from cron.

@requires: sys
@requires: os

@requires: Modules.Config
@requires: Modules.Output
@requires: Modules.Remote
@requires: Modules.Update
"""

import sys # sys.argv
import os  # os.chdir()

from mutalyzer.config import Config
from mutalyzer import Output
from mutalyzer.Db import Remote
from mutalyzer.Db import Update

os.chdir(sys.argv[0].rsplit('/', 2)[0])

C = Config()
O = Output.Output(__file__, C.Output)
O.addMessage(__file__, -1, "INFO", "Starting UCSC mapping data update")

for i in C.Db.dbNames :
    RemoteDb = Remote(i, C.Db)
    tempfile = RemoteDb.get_Update()
    del RemoteDb

    LocalDb = Update(i, C.Db)
    LocalDb.load_Update(tempfile)

    count_Updates = LocalDb.count_Updates()
    if count_Updates :
        O.addMessage(__file__, -1, "INFO", "%i updates found" % count_Updates)
        LocalDb.backup_cdsUpdates()
        cds_Updates = LocalDb.count_cdsUpdates()
        if cds_Updates :
            O.addMessage(__file__, -1, "INFO",
                         "%i CDS updates found, backing up" % cds_Updates)
        LocalDb.merge_cdsUpdates()
    #if
    LocalDb.merge_Update()

    del LocalDb
#for

O.addMessage(__file__, -1, "INFO", "UCSC mapping data update end")

del O, C
