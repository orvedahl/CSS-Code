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
                exact[i,j,k,ith] = -r*sin(p) - 1.5*cos(p)**2*sin(t)/r1
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

    plot = True

    it1 = int(0.3*nth)
    it2 = int(0.6*nth)
    ip1 = int(0.3*nphi)
    ip2 = int(0.6*nphi)
    ir1 = int(0.3*nr)
    ir2 = int(0.6*nr)
    # R-Component
    print "\nRadial Component:"
    # varying r
    x1 = true[it1,ip1,:,ir]; y1 = curl_data[it1,ip1,:,ir]
    x2 = true[it2,ip2,:,ir]; y2 = curl_data[it2,ip2,:,ir]
    l2_1 = l2norm(x1, y1)
    l2_2 = l2norm(x2, y2)
    print "\tvarying r L2 pt1: ",l2_1
    print "\tvarying r L2 pt2: ",l2_2
    show_plot("Radial-Comp, vary radius", x1, x2, y1, y2, radius, plot)
    # varying theta
    x1 = true[:,ip1,ir1,ir]; y1 = curl_data[:,ip1,ir1,ir]
    x2 = true[:,ip2,ir2,ir]; y2 = curl_data[:,ip2,ir2,ir]
    l2_1 = l2norm(x1, y1)
    l2_2 = l2norm(x2, y2)
    print "\tvarying th L2 pt1: ",l2_1
    print "\tvarying th L2 pt2: ",l2_2
    show_plot("Radial-Comp, vary theta", x1, x2, y1, y2, theta, plot)
    # varying phi
    x1 = true[it1,:,ir1,ir]; y1 = curl_data[it1,:,ir1,ir]
    x2 = true[it2,:,ir2,ir]; y2 = curl_data[it2,:,ir2,ir]
    l2_1 = l2norm(x1, y1)
    l2_2 = l2norm(x2, y2)
    print "\tvarying phi L2 pt1: ",l2_1
    print "\tvarying phi L2 pt2: ",l2_2
    show_plot("Radial-Comp, vary phi", x1, x2, y1, y2, phi, plot)


    # Theta-Component
    print "\nTheta Component:"
    # varying r
    x1 = true[it1,ip1,:,ith]; y1 = curl_data[it1,ip1,:,ith]
    x2 = true[it2,ip2,:,ith]; y2 = curl_data[it2,ip2,:,ith]
    l2_1 = l2norm(x1, y1)
    l2_2 = l2norm(x2, y2)
    print "\tvarying r L2 pt1: ",l2_1
    print "\tvarying r L2 pt2: ",l2_2
    show_plot("Theta-Comp, vary radius", x1, x2, y1, y2, radius, plot)
    # varying theta
    x1 = true[:,ip1,ir1,ith]; y1 = curl_data[:,ip1,ir1,ith]
    x2 = true[:,ip2,ir2,ith]; y2 = curl_data[:,ip2,ir2,ith]
    l2_1 = l2norm(x1, y1)
    l2_2 = l2norm(x2, y2)
    print "\tvarying th L2 pt1: ",l2_1
    print "\tvarying th L2 pt2: ",l2_2
    show_plot("Theta-Comp, vary theta", x1, x2, y1, y2, theta, plot)
    # varying phi
    x1 = true[it1,:,ir1,ith]; y1 = curl_data[it1,:,ir1,ith]
    x2 = true[it2,:,ir2,ith]; y2 = curl_data[it2,:,ir2,ith]
    l2_1 = l2norm(x1, y1)
    l2_2 = l2norm(x2, y2)
    print "\tvarying phi L2 pt1: ",l2_1
    print "\tvarying phi L2 pt2: ",l2_2
    show_plot("Theta-Comp, vary phi", x1, x2, y1, y2, phi, plot)


    # Phi-Component
    print "\nPhi Component:"
    # varying r
    x1 = true[it1,ip1,:,iphi]; y1 = curl_data[it1,ip1,:,iphi]
    x2 = true[it2,ip2,:,iphi]; y2 = curl_data[it2,ip2,:,iphi]
    l2_1 = l2norm(x1, y1)
    l2_2 = l2norm(x2, y2)
    print "\tvarying r L2 pt1: ",l2_1
    print "\tvarying r L2 pt2: ",l2_2
    show_plot("Phi-Comp, vary radius", x1, x2, y1, y2, radius, plot)
    # varying theta
    x1 = true[:,ip1,ir1,iphi]; y1 = curl_data[:,ip1,ir1,iphi]
    x2 = true[:,ip2,ir2,iphi]; y2 = curl_data[:,ip2,ir2,iphi]
    l2_1 = l2norm(x1, y1)
    l2_2 = l2norm(x2, y2)
    print "\tvarying th L2 pt1: ",l2_1
    print "\tvarying th L2 pt2: ",l2_2
    show_plot("Phi-Comp, vary theta", x1, x2, y1, y2, theta, plot)
    # varying phi
    x1 = true[it1,:,ir1,iphi]; y1 = curl_data[it1,:,ir1,iphi]
    x2 = true[it2,:,ir2,iphi]; y2 = curl_data[it2,:,ir2,iphi]
    l2_1 = l2norm(x1, y1)
    l2_2 = l2norm(x2, y2)
    print "\tvarying phi L2 pt1: ",l2_1
    print "\tvarying phi L2 pt2: ",l2_2
    show_plot("Phi-Comp, vary phi", x1, x2, y1, y2, phi, plot)

    print

    return

def show_plot(title, true1, true2, data1, data2, x, plot):

    if (not plot): return

    pylab.clf()

    pylab.title(title)

    l1 = "-"
    l2 = "--"
    colt = "r"  # true data
    cold = "b"  # data
    colr = "g"  # residual

    pylab.plot(x, true1, label='true-1', color=colt, linestyle=l1)
    pylab.plot(x, true2, label='true-2', color=colt, linestyle=l2)

    pylab.plot(x, data1, label='data-1', color=cold, linestyle=l1)
    pylab.plot(x, data2, label='data-2', color=cold, linestyle=l2)
    pylab.legend(loc='upper left')

    pylab.twinx()
    pylab.plot(x, (data1-true1)/true1, label='data-1 % err', color=colr, 
               linestyle=l1)
    pylab.plot(x, (data2-true2)/true2, label='data-2 % err', color=colr, 
               linestyle=l2)
    pylab.legend(loc='upper right')

    pylab.show()

    return

def l2norm(x, y):

    l2 = numpy.sqrt(numpy.sum((x - y)**2))

    return l2

if __name__ == "__main__":

   test()
