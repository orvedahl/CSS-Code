#!/usr/bin/env python
#
# calculate a shell average of data
#
# i.e. data(theta, phi, r, nq) ---> data(r, nq)_rms
#
# Orvedahl R. 9-9-2014
#

import sys
import numpy
import integrate

def shell_avg_rms(data, theta, phi, method="simp"):

    # idea:
    #                        integral( data**2 sin(th) dth dphi )
    #      avgdata(r)**2 = ----------------------------------------
    #                           integral( sin(th) dth dphi )
    # 
    #      the r**2 in the area elements cancel each other
    #
    #      integral( sin(th) dth dphi ) = (phi2 - phi1)*(cos(th1) - cos(th2))
    #         cos(th) is reversed since int(sin) = -cos: -----/\---------/\

    # get integral bounds
    phihi = phi[-1]
    philo = phi[0]
    thhi  = theta[-1]
    thlo  = theta[0]
    dphi  = phihi - philo
    dth   = thhi - thlo

    stheta = numpy.sin(theta)

    # calculate denominator integral
    area = (phihi - philo)*(numpy.cos(thlo) - numpy.cos(thhi))
    invrootarea = numpy.sqrt(1./area)

    nr = numpy.size(data, axis=2)
    nq = numpy.size(data, axis=3)
    nth = len(theta)
    nphi = len(phi)

    avgdata = numpy.zeros((nr, nq))
    tmp = numpy.zeros((nth, nr, nq))

    # first integrate over phi
    nslabs = nphi - 1
    for k in range(nq):
        for j in range(nr):
            for i in range(nth):
                val = integrate.integrate1D(data[i,:,j,k]*data[i,:,j,k], phi,
                                            dphi, nslabs, method=method)
                tmp[i,j,k] = val

    # square tmp and multiply by sin(theta)
    for j in range(nq):
        for i in range(nr):
            tmp[:,i,j] = stheta*tmp[:,i,j]*tmp[:,i,j]

    # second integrate over theta
    nslabs = nth - 1
    for j in range(nq):
        for i in range(nr):
            # this is an RMS value so take sqrt
            val = integrate.integrate1D(tmp[:,i,j], theta, dth, nslabs, 
                                      method=method)
            avgdata[i,j] = numpy.sqrt(val)

    # normalize by the area
    avgdata = invrootarea*avgdata

    return avgdata


