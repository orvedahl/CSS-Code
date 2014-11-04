#!/usr/bin/env python
#
# python module for calculating the vector curl in spherical coordinates
#      v(r, th, phi) --> curl_v(r, th, phi)
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
def vectorcurl(iter, case, mag=False, test=None):

    # read data from checkpoint
    # data is of size (nth, nphi, nr, 3) where last indices are
    #    0, w --> u_phi
    #    1, v --> u_theta
    #    2, u --> u_r
    iphi = 0
    ith  = 1
    irad = 2
    if (test == None):
        data, rad, theta, phi, time, header, ierr = read_ck.read_checkpoint(\
                                                   iter, case, mag_curl=mag, 
                                                   vel_curl=True)
        if (ierr):
            print "\nERROR: could not extract data from checkpoint"
            print "\theader file: "+header
            print
            sys.exit(2)
    else:
        # test is list of data, radius, theta & phi arrays
        data  = test[0]
        rad   = test[1]
        theta = test[2]
        phi   = test[3]

    # get inverses
    dri = 1./(rad[-1] - rad[0])
    dti = 1./(theta[-1] - theta[0])
    dpi = 1./(phi[-1] - phi[0])

    # get trig funcs
    sines = numpy.sin(theta)
    cosines = numpy.cos(theta)
    cosec = 1./sines
    rinv = 1./rad

    nth = len(theta)
    nphi = len(phi)
    nr = len(rad)

    curl = numpy.empty((nth, nphi, nr, 3))

    # Theta derivatives
    print "\nBegin Theta..."
    for ir in range(nr):
        ri = rinv[ir]
        for ip in range(nphi):
            # tmp = sin(th)*Vphi
            tmp = sines*data[:,ip,ir,iphi]

            # fill boundary conditions: Neumann on x1 & x2 (ibc=4)
            bc1 = -cosines[0]*data[0,ip,ir,irad]/dti
            bc2 = -cosines[-1]*data[-1,ip,ir,irad]/dti

            # r-component curl = 1/(r*sin(th)) * d(Vphi*sin(th))/dtheta
            curl[:,ip,ir,irad] = ri*cosec*compact_fd6(dti, tmp, bc1, bc2, 4)

            # tmp = Vr
            tmp = data[:,ip,ir,irad]

            # boundaries: dtmp/dtheta = 0 (ibc=0)
            bc1 = 0.
            bc2 = 0.

            # phi-component curl = -1/r * dVr/dtheta
            curl[:,ip,ir,iphi] = -ri*compact_fd6(dti, tmp, bc1, bc2, 0)

        per = str(round(100.*ir/float(nr-1),1))
        # the "\r" returns to the beginning of the line
        # the extra space at the end is to place the cursor out of the way
        sys.stdout.write("\r" + "\t% Completed: "+per+20*" ")
        # flush removes the previous line
        sys.stdout.flush()

    # Phi derivatives
    print "\nBegin Phi..."
    for ir in range(nr):
        ri = rinv[ir]
        for it in range(nth):
            # tmp = Vr
            tmp = data[it,:,ir,irad]

            # boundaries: Neumann on x1 & x2 (ibc=4)
            # domain type is internal (dtype=3)
            bc1 = data[it,nphi-4:nphi-1+1,ir,irad]
            bc2 = data[it,0:3+1,ir,irad]

            # theta-component curl = 1/(r*sin(th)) * dVr/dphi
            curl[it,:,ir,ith] = ri*cosec[it]*compact_fd6(dpi, tmp, bc1, bc2, 
                                                       4, dtype=3)

            # tmp = Vth
            tmp = data[it,:,ir,ith]

            # fill boundary: Neumann on x1 & x2 (ibc=4)
            # internal domain type (dtype=3)
            bc1 = data[it,nphi-4:nphi-1+1,ir,ith]
            bc2 = data[it,0:3+1,ir,ith]

            # r-component curl = -1/(r*sin(th)) * dVth/dphi
            # r-component has been partially filled already so need "+="
            curl[it,:,ir,irad] += -ri*cosec[it]*compact_fd6(dpi, tmp, bc1, bc2,
                                                         4, dtype=3)

        per = str(round(100.*ir/float(nr-1),1))
        sys.stdout.write("\r" + "\t% Completed: "+per+20*" ")
        sys.stdout.flush()

    # Radial derivatives
    print "\nBegin Radial..."
    for ip in range(nphi):
        for it in range(nth):
            # tmp = r*Vphi
            tmp = rad*data[it,ip,:,iphi]

            # boundaries: Neumann on x1 & x2 (ibc=4)
            bc1 = -data[it,ip,0,irad]/dri
            bc2 = -data[it,ip,-1,irad]/dri

            # theta-component curl = -1/r * d(r*Vphi)/dr
            # theta-component has been partially filled already so need "+="
            curl[it,ip,:,ith] += -rinv*compact_fd6(dri, tmp, bc1, bc2, 4)

            # tmp = r*Vth
            tmp = rad*data[it,ip,:,ith]

            # boundaries: dtmp/dr = 0 (ibc=0)
            bc1 = -data[it,ip,0,ith]/dri
            bc2 = -data[it,ip,-1,ith]/dri

            # phi-component curl = -1/r * d(rVth)/dr
            # phi-component has been partially filled already so need "+="
            curl[it,ip,:,iphi] += rinv*compact_fd6(dri, tmp, bc1, bc2, 0)

        per = str(round(100.*ip/float(nphi-1),1))
        sys.stdout.write("\r" + "\t% Completed: "+per+20*" ")
        sys.stdout.flush()

    print "\nComplete\n"

    return curl, rad, theta, phi


