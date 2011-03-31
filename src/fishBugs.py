#!/usr/bin/python

from Modules import Config

myConfig = Config.Config()
handle = open(myConfig.Output.log, "r")

scanning = False
line = handle.readline()
while line :
    if not scanning :
        if " Received " in line :
            message = line
            scanning = True
        #if
    #if
    else :
        if " Received " in line :
            print message,
            scanning = False
        #if
        if " Finished " in line :
            scanning = False
    #else
    line = handle.readline()
#while
handle.close()
