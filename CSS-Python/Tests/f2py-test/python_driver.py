#!/usr/bin/env python
#
# play with f2py and try to incorporate "use module" in the fortran code
#
# R. Orvedahl 11-7-2014

import numpy
import functionf90
import integratef90

def main():

    nx = 127

    xlo = 0.0
    xhi = 1.0

    method = "simp"

    # test the "use stuff" function
    print "\nTEST 'use module'"
    #                 module file - module name  - function
    numerical, exact = functionf90.functionf90.def_func_and_integrate(\
                                                 nx, xlo, xhi, method)

    print nx, numerical, exact, abs(exact-numerical)

    # directly call the integration routine
    dx = (xhi-xlo)/float(nx-1)
    x = numpy.zeros((nx))
    fx = numpy.zeros((nx))
    for i in range(nx):
        x[i] = xlo + float(i)*dx
        fx[i] = numpy.sin(numpy.pi*x[i])

    print "\nTEST direct call"
    #           module file - module name   - function
    numerical = integratef90.integratef90.integrate1d(fx, x, dx, method=method)

    print nx, numerical, exact, abs(exact-numerical)

    return

main()

