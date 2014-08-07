#!/usr/bin/env python
#
# python module for calculating the vector curl
#
# 8-6-2014 Orvedahl R.
#

import sys
import string
import numpy
from derivatives import compact_fd6
import read_checkpoint as read_ck
import defaults


#######################################################################
# vector curl
#######################################################################
def vectorcurl(iter, case, mag=False):

    # read data from checkpoint
    data, rad, theta, phi, time, header, ierr = read_ck.read_checkpoint(\
                                                   iter, case, mag_curl=mag, 
                                                   vel_curl=True)
    if (ierr):
        print "\nERROR: could not extract data from checkpoint"
        print "\theader file: "+header
        print
        sys.exit(2)

    # get inverses
    dri = 1./(rad[-1] - rad[0])
    dti = 1./(theta[-1] - theta[0])
    dpi = 1./(phi[-1] - phi[0])

    # get trig funcs
    sines = numpy.sin(theta)
    cosines = numpy.cos(theta)
    cosec = 1./sines
    cotan = cosines/sines
    rinv = 1./rad

    nth = len(theta)
    nphi = len(phi)
    nr = len(rad)

    curl = numpy.empty((nth, nphi, nr, 3))

    # Theta
    print "\nBegin Theta..."
    for ir in range(nr):
        ri = rinv[ir]
        for ip in range(nphi):
            tmp = sines*data[:,ip,ir,1]
            bc1 = -cosines[0]*data[0,ip,ir,2]/dti
            bc2 = -cosines[-1]*data[-1,ip,ir,2]/dti
            curl[:,ip,ir,2] = ri*cosec*compact_fd6(dti, tmp, bc1, bc2, 4)
            tmp = data[:,ip,ir,2]
            bc1 = 0.
            bc2 = 0.
            curl[:,ip,ir,1] = -ri*compact_fd6(dti, tmp, bc1, bc2, 0)

        per = str(round(100.*ir/float(nr-1),1))
        # the "\r" returns to the beginning of the line
        # the extra space at the end is to place the cursor out of the way
        sys.stdout.write("\r" + "\t% Completed: "+per+20*" ")
        # flush removes the previous line
        sys.stdout.flush()

    # Phi
    print "\nBegin Phi..."
    for ir in range(nr):
        ri = rinv[ir]
        for it in range(nth):
            tmp = data[it,:,ir,2]
            bc1 = data[it,nphi-4:nphi-1+1,ir,2]
            bc2 = data[it,0:3+1,ir,2]
            curl[it,:,ir,0] = ri*cosec[it]*compact_fd6(dpi, tmp, bc1, bc2, 
                                                       4, dtype=3)
            tmp = data[it,:,ir,0]
            bc1 = data[it,nphi-4:nphi-1+1,ir,0]
            bc2 = data[it,0:3+1,ir,0]
            curl[it,:,ir,2] += -ri*cosec[it]*compact_fd6(dpi, tmp, bc1, bc2,
                                                         4, dtype=3)
        per = str(round(100.*ir/float(nr-1),1))
        sys.stdout.write("\r" + "\t% Completed: "+per+20*" ")
        sys.stdout.flush()

    # Radial
    print "\nBegin Radial..."
    for ip in range(nphi):
        for it in range(nth):
            tmp = rad*data[it,ip,:,2]
            bc1 = -data[it,ip,0,2]/dri
            bc2 = -data[it,ip,-1,2]/dri
            curl[it,ip,:,0] += -rinv*compact_fd6(dri, tmp, bc1, bc2, 4)
            tmp = rad*data[it,ip,:,0]
            bc1 = -data[it,ip,0,0]/dri
            bc2 = -data[it,ip,-1,0]/dri
            curl[it,ip,:,1] += rinv*compact_fd6(dri, tmp, bc1, bc2, 0)

        per = str(round(100.*ip/float(nphi-1),1))
        sys.stdout.write("\r" + "\t% Completed: "+per+20*" ")
        sys.stdout.flush()

    print "\nComplete\n"

    return curl, rad, theta, phi


