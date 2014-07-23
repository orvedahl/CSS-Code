#!/usr/bin/env python
#
# read Shell_Slice files
#
# Orvedahl R. 7-3-2014
#

import unformatted_read
import get_quantity_name as get_name
import numpy

#####################################################################
# read AZ_Avgs files
#####################################################################
def read_shell_slice(filename, rThetaPhiQuantities=False, nonuni=False, 
                     recmax=None):

    data, ierr = unformatted_read.read_f90(filename)

    if (not ierr):

        ldata = len(data)

        # extract shape of data cube and reshape data cube
        nradii = int(data[0])
        nr     = int(data[1])
        nth    = int(data[2])
        nphi   = int(data[3])
        nq     = int(data[4])

        # print nth, nr, nq
        print "nth, nphi, nradii, nq = ",nth, nphi, nradii, nq

        # number of records
        nrecs = ldata / (nth*(nphi+1)*nradii*nq)
        if (nrecs == 0):
            nrecs = 1

        # the order='F' preserves the Fortran ordering
        data = numpy.reshape(data, (nth, nphi+1, nradii, nq, nrecs), order='F')

        # extract special quantities
        r1   = data[5, 0, 0, 0, 0]
        r2   = data[6, 0, 0, 0, 0]
        th1  = data[7, 0, 0, 0, 0]
        th2  = data[8, 0, 0, 0, 0]
        phi1 = data[9, 0, 0, 0, 0]
        phi2 = data[10,0, 0, 0, 0]

        if (rThetaPhiQuantities):

            # get phi and theta
            phi = (phi2-phi1)*numpy.array(range(nphi))/(nphi-1.) + phi1
            theta = (th2-th1)*numpy.array(range(nth ))/(nth-1. ) + th1

            # calculate theta, quantities and radius
            quantities = data[11:nq+10+1, 0, 0, 0, 0].astype(int)

            print quantities
            q_names = get_name.get_quantity_name(quantities)
            print q_names

            r = numpy.array(range(nr))/float(nr-1)
            if (nonuni):
                radius = (r2-r1)*(1.-numpy.exp(-2.*r))/(1.-numpy.exp(-2.)) + r1
            else:
                radius = (r2-r1)*r + r1

            r_ind = data[nq+11:nradii+nq+10+1, 0, 0, 0, 0].astype(int)
            radii = radius[r_ind - 1]

        if (recmax != None and recmax >= 1):
            record = int(recmax)
        else:
            record = nrecs

        if (record == 0):
            record = 1

        print "number of records = ", record

        # allocate shell_slice
        shell_slice_values = numpy.empty((nth, nphi, nradii, nq, record))

        # store values
        shell_slice_values[:,:,:,:,0:record] = data[:,1:nphi+1, :, :, 0:record]

        if (rThetaPhiQuantities):
            return shell_slice_values, record, theta, phi, radii, \
                   quantities, q_names, ierr

        else:

            return shell_slice_values, record, ierr

    else:

        if (rThetaPhiQuantities):
            return None, None, None, None, None, None, None, ierr
        else:
            return None, None, ierr

