#!/usr/bin/env python
#
# return average quantity of file and over records
#
# Orvedahl R. 6-18-2014
#

import sys
import numpy
from read_az_avg import *

#####################################################################
# average over files
#####################################################################
def average_over_files(files, nfiles, nr, nth, nq):

    result = numpy.zeros((nth, nr, nq))

    print "Averaging"

    if (nfiles > 1):

        print "\n\tPercent Complete:"

        for i in range(nfiles):

            data, nrecs, ierr = read_az_avg(files[i])
            if (ierr):
                print "\n\tFile Does Not Appear To Be Binary:"
                print "\t\t"+files[i]
                sys.exit(2)

            for j in range(nq):

                result[:,:,j] += average_quantity_over_records(data, nr, 
                                                               nth, j, nrecs)

            print "\t\t"+str((100.*(i+1))/nfiles)

        result = result / float(nfiles)

        print "Finished averaging and reading"

    else:

        if (type(files) == str):
            data, nrecs, ierr = read_az_avg(files)
            if (ierr):
                print "\n\tFile Does Not Appear To Be Binary:"
                print "\t\t"+files
                sys.exit(2)

        elif (type(files) == list):
            data, nrecs, ierr = read_az_avg(files[0])
            if (ierr):
                print "\n\tFile Does Not Appear To Be Binary:"
                print "\t\t"+files[0]
                sys.exit(2)

        else:
            print "\n\tLook at data type of 'files' in average_over_files\n"
            sys.exit(2)

        for j in range(nq):

            result[:,:,j] = average_quantity_over_records(data, nr, nth, 
                                                          j, nrecs)

    return result


#####################################################################
# average a quantity over all records
#####################################################################
def average_quantity_over_records(data, nr, nth, iq, nrecs):

    result = numpy.zeros((nth, nr))

    for i in range(nrecs):

        result[:,:] += data[:,:,iq,i]

    result = result / float(nrecs)

    return result

