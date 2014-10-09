#!/usr/bin/env python
#
# calculate a shell average of data
#
# i.e. data(theta, phi, r, nq) ---> data(r, nq)
#
# Orvedahl R. 9-9-2014
#

import sys
import numpy
import integrate

def shell_avg(data, theta, phi, numq=None,  method="simp"):

    # idea:
    #                     integral( data sin(th) dth dphi )
    #      avgdata(r) = -------------------------------------
    #                        integral( sin(th) dth dphi )
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
    dphi  = phi[1] - phi[0]
    dth   = theta[1] - theta[0]

    stheta = numpy.sin(theta)

    # calculate denominator integral
    area = (phihi - philo)*(numpy.cos(thlo) - numpy.cos(thhi))
    invarea = 1./area

    nth = len(theta)
    nphi = len(phi)
    nr = numpy.size(data, axis=2)

    # sometimes we want an average of only one quantity
    if (numq==None):
        nq = numpy.size(data, axis=3)
    else:
        nq = 1
        data = numpy.reshape(data, (nth, nphi, nr, nq))

    avgdata = numpy.zeros((nr, nq))
    tmp = numpy.zeros((nth, nr, nq))

    # first: multiply by sin(theta)
    for k in range(nq):
        for j in range(nr):
            for i in range(nphi):
                data[:,i,j,k] = stheta*data[:,i,j,k]

    # second: integrate over phi
    for k in range(nq):
        for j in range(nr):
            for i in range(nth):
                val = integrate.integrate1D(data[i,:,j,k], phi, dphi, 
                                            method=method)
                tmp[i,j,k] = val

    # third: integrate over theta
    for j in range(nq):
        for i in range(nr):
            val = integrate.integrate1D(tmp[:,i,j], theta, dth, method=method)
            avgdata[i,j] = val

    # normalize by the area
    avgdata = invarea*avgdata

    return avgdata


