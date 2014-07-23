#!/usr/bin/env python
#
# read Shell_Avg files
#
# Orvedahl R. 6-16-2014
#

import unformatted_read
import get_quantity_name as get_name
import numpy

#####################################################################
# read Shell_Avgs files
#####################################################################
def read_shell_avg(filename, rQuantities=False, oldfile=False):

    data, ierr = unformatted_read.read_f90(filename)

    if (not ierr):

        ldata = len(data)

        # extract shape of data cube and reshape data cube
        nr = int(data[0])
        nq = int(data[1])

        # print nr, nq
        print "nr, nq = ",nr, nq

        if (oldfile):
            iold = 1
        else:
            iold = 2

        # number of records
        nrecs = ldata / (nr*(nq+iold))
        if (nrecs == 0):
            nrecs = 1

        print "number of records = ", nrecs

        # the order='F' preserves the Fortran ordering
        data = numpy.reshape(data, (nr, nq+iold, nrecs), order='F')

        # extract special quantities
        r1  = data[2, 0, 0]
        r2  = data[3, 0, 0]

        # allocate shell averages
        shell_avg_values = numpy.empty((nr, nq, nrecs))

        if (oldfile):
            shell_avg_values[:,:,:] = data[:,1:nq+1,:]
        else:
            shell_avg_values[:,:,:] = data[:,2:nq+1+1,:]

        if (rQuantities):

            # get radius
            if (oldfile):
                radius = (r2-r1)*numpy.array(range(nr))/float(nr-1) + r1
            else:
                radius = data[:, 1, 0]

            # calculate quantities and radius
            quantities = data[4:nq+3+1, 0, 0].astype(int)

            print quantities
            q_names = get_name.get_quantity_name(quantities)
            print q_names

            return shell_avg_values, nrecs, radius, quantities, q_names, ierr

        else:

            return shell_avg_values, nrecs, ierr

    else:

        if (rQuantities):
            return None, None, None, None, None, ierr
        else:
            return None, None, ierr

