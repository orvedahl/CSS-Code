#!/usr/bin/env python
#
# calculate Rossby number from checkpoint files
#
# Orvedahl R. 8-7-2014
#

import sys
import numpy
import vectorcurl
from read_checkpoint import *

def rossby_from_chk(iter, case, omega0, method='adv', return_rad=False):

    if (method == 'adv'):
        # Use the advective term:
        #      Ro =  U / (2*L*Omega) where U = magnitude of vel

        # get velocity, dont need rho or s --> vel_curl=True
        data, radius, theta, phi, time, header, ierr = read_checkpoint(\
                                                  iter, case, vel_curl=True)
        if (ierr):
            print "\n---ERROR: could not read checkpoint data---"
            print "\theader file: "+header+"\n"
            sys.exit(2)

        # average velocities over theta and phi
        # data is of size (nth, nphi, nr, 3) where last indices are
        #    0, w --> u_phi
        #    1, v --> u_theta
        #    2, u --> u_r
        vphi = numpy.mean(numpy.mean(data[:,:,:,0], axis=0), axis=0)
        vth  = numpy.mean(numpy.mean(data[:,:,:,1], axis=0), axis=0)
        vr   = numpy.mean(numpy.mean(data[:,:,:,2], axis=0), axis=0)

        # magnitude of velocity
        U = numpy.sqrt(vr*vr + vth*vth + vphi*vphi)

        nr = len(radius)
        rossby = numpy.zeros((nr))

        radial_domain = radius[-1] - radius[0]

        # get Omega(r)
        omega = numpy.zeros((nr))
        omega[:] = omega0 + vphi[:]/radius[:]

        # calculate Rossby number
        rossby[:] = U / (2. * radial_domain * omega)

    elif (method == 'ens'):
        # Use the Enstrophy
        #      Ro =  |curl(v)| / Omega

        # get curl
        curl, radius, theta, phi = vectorcurl.vectorcurl(iter, case)

        # get magnitude of curl --> sqrt(enstrophy)
        ens = numpy.sqrt(curl[:,:,:,0]*curl[:,:,:,0] +
                         curl[:,:,:,1]*curl[:,:,:,1] +
                         curl[:,:,:,2]*curl[:,:,:,2])

        nr = len(radius)
        rossby = numpy.zeros((nr))
        for ir in range(nr):
            rossby[ir] = numpy.average(ens[:,:,ir])

        # divide by rotation rate of the frame
        rossby[:] = rossby[:]/omega0

    else:
        print "\n---ERROR: unknown method in rossby_from_chk\n"
        sys.exit(2)

    if (return_rad):
        return rossby, radius
    else:
        return rossby


