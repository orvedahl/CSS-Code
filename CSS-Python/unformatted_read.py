#!/usr/bin/env python
#
# read an unformatted fortran data file that was generated with:
#       open(unit, file, form='unformatted', access='direct', recl=irec)
#	...
#	write(unit, rec=rec_number) some_array
#
# Orvedahl R. 6-16-2014
#

import numpy

#####################################################################
# read unformatted f90 file, return data as 1D array
#####################################################################
def read_f90(filename):

    binary = is_binary(filename)

    # want file to be binary, if it is not, throw error
    if (binary):
        ierr = False
    else:
        ierr = True

    if (not ierr):
        f = open(filename, 'rb')

        data1D = numpy.fromfile(f, dtype='float64')

        f.close()

        return data1D, ierr

    else:

        return None, ierr

#####################################################################
# check to see if file truly is a binary file
#####################################################################
def is_binary(filename):

    # return True if filename appears to be binary
    # binary is considered True if it contains NULL byte = '\0'
    # FIXME: this apparently incorrectly reports UTF-16 as binary

    f = open(filename, 'rb')
    for block in f:
        if ('\0' in block):
            f.close()
            return True
    f.close()
    return False

