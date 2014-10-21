#!/usr/bin/env python
#
# python module for testing the shell_avg code
#
# 10-9-2014 Orvedahl R.
#

import os
import numpy
import pylab
from numpy import sin, cos
import relative_import
unit_test_dir = os.path.dirname(os.path.realpath(__file__))
relative_import.append_path(unit_test_dir+"/..") # import modules from ./.. dir
import shell_avg

# function to average
def func(rad, theta, phi):
    nr = len(rad)
    nth = len(theta)
    nphi = len(phi)
    nq = 1

    #                i    j     k   l
    f = numpy.empty((nth, nphi, nr, nq))

    for l in range(nq):
        for k in range(nr):
            r = rad[k]
            for j in range(nphi):
                p = phi[j]
                for i in range(nth):
                    t = theta[i]

                    # version 1
                    #f[i,j,k,l] = r*sin(t) + r*r*cos(p)*cos(t)

                    # version 2
                    f[i,j,k,l] = r*sin(t)*cos(p)

    return f

# exact shell avg
def exact(r, phi, th):
    # get bounds
    phihi = phi[-1]
    philo = phi[0]
    thhi = th[-1]
    thlo = th[0]

    # version 1:
    # numerator
    #a = r*(phihi - philo)
    #a *= (0.5*thhi - 0.25*sin(2.*thhi)) - \
    #     (0.5*thlo - 0.25*sin(2.*thlo))

    #b = r*r*(sin(phihi) - sin(philo))
    #b *= -0.25*(cos(2.*thhi) - cos(2.*thlo))

    # version 2:
    # numerator
    a = r*(sin(phihi) - sin(philo))
    a *= (0.5*(thhi - thlo) - 0.25*(sin(2.*thhi) - sin(2.*thlo)))
    b = 0.

    # denominator
    c = (philo - phihi)*(numpy.cos(thhi) - numpy.cos(thlo))

    return (a + b) / c

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
    avgdata = shell_avg.shell_avg(data, theta, phi, method='simp')

    # only consider first quantity
    avgdata = avgdata[:,0]

    true = exact(radius, phi, theta)

    l2norm = numpy.sqrt(numpy.sum((true-avgdata)**2))

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
