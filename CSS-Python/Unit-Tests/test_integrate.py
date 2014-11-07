#!/usr/bin/env python
#
# python module for testing the integrals
#
# 9-27-2014 Orvedahl R.
#

import os
import relative_import
unit_test_dir = os.path.dirname(os.path.realpath(__file__)) # path of this file
relative_import.append_path(unit_test_dir+"/..")
#from integrate import integrate1D as int1D  # pure python version
import integratef90                          # fortran f2py version
int1D = integratef90.integratef90.integrate1d   # f90: case insensitive
import numpy
from numpy import exp, cos, sin
import pylab
import time

def func(x):
    return numpy.sin(numpy.pi*x)
    #return -2.0*0.1*x*exp(-0.1*x*x)*cos(x) - exp(-0.1*x*x)*sin(x)

def exact(xlo, xhi):
    return (1./numpy.pi)*(-numpy.cos(numpy.pi*xhi) + numpy.cos(numpy.pi*xlo))
    #def f(x):
    #    return exp(-0.1*x*x)*cos(x)
    #return f(xhi) - f(xlo)

def test():

    starttime = time.time()

    xlo = 0.
    xhi = 1.0

    # choose number of slabs to use
    N = [3, 7, 15, 31, 63, 127]

    dx  = []
    err = []

    for n in N:

        delta = (xhi - xlo) / n

        # inclusive: [xlo, xhi] & len(x) = n+1
        #  --> very important to be inclusive! will be 2nd order otherwise
        x = numpy.arange(xlo, xhi+delta, delta)
        f = func(x)

        integral = int1D(f, x, delta, method='simp')

        e = abs(integral - exact(xlo, xhi))
        dx.append(delta)
        err.append(e)
        print n, delta, integral, e

    endtime = time.time()
    print "\nelapsed time: ", endtime-starttime
    print

    dx = numpy.array(dx)
    err = numpy.array(err)

    pylab.plot(dx, err, linestyle="", marker="x", color='r')

    #pylab.plot(dx, err[0]*(dx/dx[0])**1,color='r',linestyle="--",label="1st")
    pylab.plot(dx, err[0]*(dx/dx[0])**2,color='b',linestyle="--",label="2nd")
    pylab.plot(dx, err[0]*(dx/dx[0])**4,color='g',linestyle="--",label="4th")

    pylab.legend(loc='lower right')

    ax = pylab.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')

    pylab.xlabel(r"$\Delta x$")
    pylab.ylabel(r"Abs Err")

    pylab.show()

    return


if __name__ == "__main__":

   test()

