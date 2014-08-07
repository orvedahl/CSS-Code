#!/usr/bin/env python
#
# calculate Rossby number from checkpoint files
#
# Orvedahl R. 8-7-2014
#

import numpy
import vectorcurl

def rossby_from_chk(iter, case, omega):

    # get curl
    curl, rad, theta, phi = vectorcurl.vectorcurl(iter, case)

    # get magnitude of curl
    ens = numpy.sqrt(curl[:,:,:,0]*curl[:,:,:,0] +
                     curl[:,:,:,1]*curl[:,:,:,1] +
                     curl[:,:,:,2]*curl[:,:,:,2])

    nr = len(rad)
    rossby = numpy.zeros((nr))
    for ir in range(nr):
        rossby[ir] = numpy.average(ens[:,:,ir])

    if (type(omega) == float):
        rossby[:] = rossby[:]/omega
    else:
        rossby[:] = rossby[:]/omega[:]

    return rossby


