#!/usr/bin/python

"""
@todo: documentation
"""

import time

class Misc() :
    """
    @todo: documentation
    """

    def ID(self) :
        """
        Generates an ID using time()
        @todo: documentation
        
        @return:
        @rtype: 
        """

        IDsPerSec = 100

        time.sleep(1.0 / IDsPerSec)
        return int(time.time() * IDsPerSec)
    #ID
#Misc
