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
from vectorcurl import *

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

    #                     i    j     k   l
    f     = numpy.empty((nth, nphi, nr, nq))
    exact = numpy.empty((nth, nphi, nr, nq))

    for k in range(nr):
        r = rad[k]
        r2 = r*r
        r1 = numpy.sqrt(r)
        for j in range(nphi):
            p = phi[j]
            for i in range(nth):
                t = theta[i]

                f[i,j,k,ir]   = r2*sin(t)*cos(p)
                f[i,j,k,ith]  = cos(p)*sin(t)*sin(t)/r
                f[i,j,k,iphi] = r1*cos(p)*cos(p)*sin(t)

                # exact curl
                exact[i,j,k,ir] = 2.*cos(p)**2*cos(t)/r1 + sin(p)*sin(t)/r2
                exact[i,j,k,ith] = -r*sin(t) - 1.5*cos(p)**2*sin(t)/r1
                exact[i,j,k,iphi] = -r*cos(t)*cos(p)

    return f, exact

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

    nr   = len(radius)
    nth  = len(theta)
    nphi = len(phi)
    print nr,nth,nphi

    # data has dimensions of (nth, nphi, nr, nq=1)
    print "\nSet Up Function & Exact Curl..."
    data, true = func(radius, theta, phi)

    # avgdata has dimensions of (nr, nq)
    iter = None
    case = None
    print "\nTake Curl..."
    curl_data = vectorcurl(iter, case, test=[data, radius, theta, phi])
    curl_data = curl_data[0] # vectorcurl returns 4 objects, we only 
                             # care about the first one

    it1 = int(0.3*nth)
    it2 = int(0.6*nth)
    ip1 = int(0.3*nphi)
    ip2 = int(0.6*nphi)
    ir1 = int(0.3*nr)
    ir2 = int(0.6*nr)
    # R-Component
    print "\nRadial Component:"
    # varying r
    l2_1 = l2norm(true[it1,ip1,:,ir], curl_data[it1,ip1,:,ir])
    l2_2 = l2norm(true[it2,ip2,:,ir], curl_data[it2,ip2,:,ir])
    print "\tvarying r L2 pt1: ",l2_1
    print "\tvarying r L2 pt2: ",l2_2
    # varying theta
    l2_1 = l2norm(true[:,ip1,ir1,ir], curl_data[:,ip1,ir1,ir])
    l2_2 = l2norm(true[:,ip2,ir2,ir], curl_data[:,ip2,ir2,ir])
    print "\tvarying th L2 pt1: ",l2_1
    print "\tvarying th L2 pt2: ",l2_2
    # varying phi
    l2_1 = l2norm(true[it1,:,ir1,ir], curl_data[it1,:,ir1,ir])
    l2_2 = l2norm(true[it2,:,ir2,ir], curl_data[it2,:,ir2,ir])
    print "\tvarying phi L2 pt1: ",l2_1
    print "\tvarying phi L2 pt2: ",l2_2

    # Theta-Component
    print "\nTheta Component:"
    # varying r
    l2_1 = l2norm(true[it1,ip1,:,ith], curl_data[it1,ip1,:,ith])
    l2_2 = l2norm(true[it2,ip2,:,ith], curl_data[it2,ip2,:,ith])
    print "\tvarying r L2 pt1: ",l2_1
    print "\tvarying r L2 pt2: ",l2_2
    # varying theta
    l2_1 = l2norm(true[:,ip1,ir1,ith], curl_data[:,ip1,ir1,ith])
    l2_2 = l2norm(true[:,ip2,ir2,ith], curl_data[:,ip2,ir2,ith])
    print "\tvarying th L2 pt1: ",l2_1
    print "\tvarying th L2 pt2: ",l2_2
    # varying phi
    l2_1 = l2norm(true[it1,:,ir1,ith], curl_data[it1,:,ir1,ith])
    l2_2 = l2norm(true[it2,:,ir2,ith], curl_data[it2,:,ir2,ith])
    print "\tvarying phi L2 pt1: ",l2_1
    print "\tvarying phi L2 pt2: ",l2_2


    # Phi-Component
    print "\nPhi Component:"
    # varying r
    l2_1 = l2norm(true[it1,ip1,:,iphi], curl_data[it1,ip1,:,iphi])
    l2_2 = l2norm(true[it2,ip2,:,iphi], curl_data[it2,ip2,:,iphi])
    print "\tvarying r L2 pt1: ",l2_1
    print "\tvarying r L2 pt2: ",l2_2
    # varying theta
    l2_1 = l2norm(true[:,ip1,ir1,iphi], curl_data[:,ip1,ir1,iphi])
    l2_2 = l2norm(true[:,ip2,ir2,iphi], curl_data[:,ip2,ir2,iphi])
    print "\tvarying th L2 pt1: ",l2_1
    print "\tvarying th L2 pt2: ",l2_2
    # varying phi
    l2_1 = l2norm(true[it1,:,ir1,iphi], curl_data[it1,:,ir1,iphi])
    l2_2 = l2norm(true[it2,:,ir2,iphi], curl_data[it2,:,ir2,iphi])
    print "\tvarying phi L2 pt1: ",l2_1
    print "\tvarying phi L2 pt2: ",l2_2


    print

    #pylab.clf()

    #pylab.plot(radius, true, label='true', color='r', linestyle='--')
    #pylab.plot(radius, avgdata, label='avg', color='b', linestyle='-')
    #pylab.legend(loc='lower left')

    #pylab.twinx()
    #pylab.plot(radius, (avgdata - true)/true, label='% err', 
               #color='g', linestyle='--')
    #pylab.legend(loc='lower right')

    #pylab.show()

    return

def l2norm(x, y):

    return numpy.sqrt(numpy.sum((x-y)**2))

if __name__ == "__main__":

   test()
