#!/usr/bin/env python
#
#############################################################################
# Ideally you should put this file in a central location, such as ${HOME}/bin,
# and add that location to your PYTHONPATH environment variable. This will 
# allow you to have only one copy of this file and be able to import it from
# any python code you write. It will still run if you do not do this, but 
# everytime you want to use the append_path function, this file will need to 
# be in whatever working directory you started from (not a huge deal for the 
# Unit-Test directory).
#
# To do this, use the following set of instructions:
#    [user@host ~]$ cd           <-- change to HOME directory
#    [user@host ~]$ mkdir bin    <-- make new directory (if it does not exist)
#
# in the proper rc file, such as .bashrc, add the following:
#    export PYTHONPATH=${HOME}/bin:$PYTHONPATH
#
# now update the rc file:
#    [user@host ~]$ source ~/.bashrc
#############################################################################
#
# import modules from a different directory by temporarily adding the 
# module's location to the PYTHONPATH. The modification to PYTHONPATH 
# lives only as long as the program, i.e. the change is not permanent
#
# R. Orvedahl 10-14-2014

import os
import sys

def append_path(location):

    if (os.path.isabs(location)):
        # location is absolute
        p = location
    else:
        # location is relative, make it absolute
        p = os.path.join(os.path.dirname(__file__), location)

    # add location of module to front of list of paths that python searches
    sys.path.insert(0, p)

    return

