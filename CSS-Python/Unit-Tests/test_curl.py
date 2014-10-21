#!/usr/bin/env python
#
# python module for testing the vector curl in spherical coordinates
#
# 10-21-2014 Orvedahl R.
#

import os
import numpy
import pylab
from numpy import sin, cos
import relative_import
unit_test_dir = os.path.dirname(os.path.realpath(__file__))
relative_import.append_path(unit_test_dir+"/..") # import modules from ./.. dir
import vectorcurl as curl

# indices from read_checkpoint file
iphi = 0
ith  = 1
ir   = 2

# function to average
def func(rad, theta, phi):
    nr = len(rad)
    nth = len(theta)
    nphi = len(phi)
    nq = 3

    #                i    j     k   l
    f = numpy.empty((nth, nphi, nr, nq))

    for k in range(nr):
        r = rad[k]
        for j in range(nphi):
            p = phi[j]
            for i in range(nth):
                t = theta[i]

                f[i,j,k,ir]   = r*r*sin(t)*cos(p)
                f[i,j,k,ith]  = cos(p)*sin(t)*sin(t)/r
                f[i,j,k,iphi] = numpy.sqrt(r)*cos(p)*cos(p)*sin(t)

    return f

# exact curl
def exact(rad, theta, phi):
    nr = len(rad)
    nth = len(theta)
    nphi = len(phi)
    nq = 3

    #                i    j     k   l
    f = numpy.empty((nth, nphi, nr, nq))

    for k in range(nr):
        r = rad[k]
        for j in range(nphi):
            p = phi[j]
            for i in range(nth):
                t = theta[i]

                f[i,j,k,ir]   = r*r*sin(t)*cos(p)
                f[i,j,k,ith]  = cos(p)*sin(t)*sin(t)/r
                f[i,j,k,iphi] = numpy.sqrt(r)*cos(p)*cos(p)*sin(t)

    return f

def test():

    # set up theta, phi, r, quantities
    rlo = 6.2e10
    rhi = 6.96e10
    thi = 40.0*numpy.pi/180.
    tlo = -40.0*numpy.pi/180.
    phi = 40.0*numpy.pi/180.
    plo = -40.0*numpy.pi/180.

    dr = (rhi - rlo)/100
    dth = (thi - tlo)/100
    dphi = (phi - plo)/100

    radius = numpy.arange(rlo, rhi+dr, dr)
    theta  = numpy.arange(tlo, thi+dth, dth)
    phi    = numpy.arange(plo, phi+dphi, dphi)

    # data has dimensions of (nth, nphi, nr, nq=1)
    data = func(radius, theta, phi)

    # avgdata has dimensions of (nr, nq)
    curl_data = curl.vectorcurl(iter, case, test=[data, radius, theta, phi])

    true = exact(radius, phi, theta)

    l2norm = numpy.sqrt(numpy.sum((true-curl_data)**2))

    print "\nL2 Norm: ",l2norm
    print

    pylab.clf()

    pylab.plot(radius, true, label='true', color='r', linestyle='--')
    pylab.plot(radius, avgdata, label='avg', color='b', linestyle='-')
    pylab.legend(loc='lower left')

    pylab.twinx()
    pylab.plot(radius, (avgdata - true)/true, label='% err', 
               color='g', linestyle='--')
    pylab.legend(loc='lower right')

    pylab.show()

    return


if __name__ == "__main__":

   test()
