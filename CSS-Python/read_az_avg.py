#!/usr/bin/env python
#
# read AZ_Avgs files
#
# Orvedahl R. 6-16-2014
#

import unformatted_read
import get_quantity_name as get_name
import numpy

#####################################################################
# read AZ_Avgs files
#####################################################################
def read_az_avg(filename, rThetaQuantities=False):

    data, ierr = unformatted_read.read_f90(filename)

    if (not ierr):

        ldata = len(data)

        # extract shape of data cube and reshape data cube
        nth = int(data[0])    # = nth
        nr = int(data[1])     # = nr + 2, i.e. data[1] = nr
        nq = int(data[2])     # = nq

        # print nth, nr, nq
        print "nth, nr, nq = ",nth, nr, nq

        # number of records
        nrecs = ldata / (nth*(nr+2)*nq)
        if (nrecs == 0):
            nrecs = 1

        print "number of records = ", nrecs

        # the order='F' preserves the Fortran ordering
        data = numpy.reshape(data, (nth, nr+2, nq, nrecs), order='F')

        # extract special quantities
        r1  = data[3, 0, 0, 0]
        r2  = data[4, 0, 0, 0]
        th1 = data[5, 0, 0, 0]
        th2 = data[6, 0, 0, 0]

        # extract the az averages
        az_avg_values = numpy.empty((nth, nr, nq, nrecs))
        az_avg_values[:,:,:,:] = data[:,2:,:,:]

        if (rThetaQuantities):
            # calculate theta, quantities and radius
            theta = (th2 - th1)*numpy.array(range(nth))/(nth - 1.0) + th1
            quantities = data[7:nq+6+1, 0, 0, 0].astype(int)
            radius = data[0:nr-1+1, 1, 0, 0]

            print quantities
            q_names = get_name.get_quantity_name(quantities)
            print q_names

            return az_avg_values, nrecs, theta, radius, \
                   quantities, q_names, ierr

        else:

            return az_avg_values, nrecs, ierr

    # just return with ierr=True
    else:

        if (rThetaQuantities):
            return None, None, None, None, None, None, ierr
        else:
            return None, None, ierr
