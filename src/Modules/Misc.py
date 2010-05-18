import time

class Misc() :
    """
    """

    def ID(self) :
        """
        """

        IDsPerSec = 100

        time.sleep(1.0 / IDsPerSec)
        return int(time.time() * IDsPerSec)
    #ID
#Misc
