#!/usr/bin/env python
#
# calculate Rossby number from checkpoint files
#
# Orvedahl R. 8-7-2014
#

import numpy
import vectorcurl

def rossby_from_chk(iter, case, omega0, return_rad=False):

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

    if (return_rad):
        return rossby, radius
    else:
        return rossby


