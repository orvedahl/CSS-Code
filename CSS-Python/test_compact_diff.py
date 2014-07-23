#!/usr/bin/env python
#
# python module for testing the derivatives
#
# 7-8-2014 Orvedahl R.
#

import derivatives as d
import numpy

def test_diff(nr, ibc, dtype):

    r1 = 0.75
    r2 = 1.
    dr = (r2-r1)/float(nr-1)
    dri = 1./(r2-r1)
    r = dr*numpy.array(range(nr))+r1

    f = numpy.sin(2.*numpy.pi*(r-r1)/(r2-r1))
    df = None
    exact1 = 2.*numpy.pi*numpy.cos(2.*numpy.pi*(r-r1)/(r2-r1))/(r2-r1)

    b1 = r1 - dr*(numpy.array(range(4))+1)
    b1 = b1[::-1] # reverse the array
    b1 = numpy.sin(2.*numpy.pi*(b1-r1)/(r2-r1))
    b2 = r2 + dr*(numpy.array(range(4))+1)
    b2 = numpy.sin(2.*numpy.pi*(b2-r1)/(r2-r1))

    if (dtype < 3):
        if (ibc == 0):
            b1 = 0.
            b2 = 0.

        if (ibc == 2):
            b1 = 0.

        if (ibc == 3):
            b2 = 0.

        if (ibc == 4):
            b1 = 2.*numpy.pi/dri/(r2-r1)
            b2 = 2.*numpy.pi/dri/(r2-r1)
            df = f.copy()
            df[:] = 0.
            df[0] = b1
            df[nr-1] = b2

    if (dtype == 2):
        b2 = r2 + dr*(numpy.array(range(4))+1)

    if (dtype == 1):
        b1 = r1 - dr*(numpy.array(range(4))+1)
        b1 = b1[::-1]

    df1 = d.compact_fd6(dri, f, b1, b2, ibc, dtype, darr=df)

    local_err = df1 - exact1
    global_err = numpy.sqrt(numpy.sum(local_err*local_err))

    print "\nglobal error 1st deriv: ",global_err

    return

