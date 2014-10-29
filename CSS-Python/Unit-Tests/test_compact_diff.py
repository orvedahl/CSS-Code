#!/usr/bin/env python
#
# python module for testing the derivatives
#
# 7-8-2014 Orvedahl R.
#

import os
import sys
import getopt
import relative_import
relative_import.append_path(os.path.dirname(os.path.realpath(__file__))+"/..")
import derivatives as d
import numpy
from numpy import sin, cos, pi
import pylab

def test_diff(nr, ibc, dtype, plot):

    r1 = 0.75
    r2 = 1.
    dr = (r2-r1)/float(nr-1)
    dri = 1./(r2-r1)
    r = dr*numpy.array(range(nr))+r1
    #dri *= len(r)-1

    f = sin(2.*pi*(r-r1)/(r2-r1))
    df = None
    exact1 = 2.*pi*cos(2.*pi*(r-r1)/(r2-r1))/(r2-r1)

    b1 = r1 - dr*(numpy.array(range(4))+1)
    b1 = b1[::-1] # reverse the array
    b1 = sin(2.*pi*(b1-r1)/(r2-r1))
    b2 = r2 + dr*(numpy.array(range(4))+1)
    b2 = sin(2.*pi*(b2-r1)/(r2-r1))

    if (dtype < 3):
        if (ibc == 0):
            b1 = 0.
            b2 = 0.

        if (ibc == 2):
            b1 = 0.

        if (ibc == 3):
            b2 = 0.

        if (ibc == 4):
            b1 = 2.*pi/dri/(r2-r1)
            b2 = 2.*pi/dri/(r2-r1)
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
    # the dr multiplier is to make it grid invariant
    global_err = numpy.sqrt(dr*numpy.sum(local_err*local_err))

    print "\nibc: ", ibc
    print "dtype: ", dtype
    print "nr:",len(r)
    print "dr:",dr
    print "L2 norm 1st deriv: ",global_err
    print "deriv at r[nr/4] : ", df1[nr/4]
    print "exact at r[nr/4] : ", exact1[nr/4]

    if (plot):
        pylab.clf()
        pylab.plot(r, df1, color='b', linestyle="-", label='df')
        pylab.plot(r, exact1, color='b', linestyle='--', label='df exact')
        pylab.ylabel("df")
        pylab.xlabel("radius")
        pylab.title("ibc %d, dtype %d" %(ibc, dtype))
        pylab.legend(loc='lower left')

        pylab.twinx()
        pylab.plot(r, local_err, color='r', label='abs error')
        pylab.ylabel("df/dr - exact")
        pylab.legend(loc='lower right')
        pylab.show()

    return

if __name__ == "__main__":

    try:
       opts, args = getopt.getopt(sys.argv[1:], "n:pi:d:", 
                                  ["nr=","plot", "ibc=", "dtype="])
    except getopt.GetoptError:
       print "\nUnknown command line arg"
       print sys.argv[1:]
       print "\nUsage:\n"
       print "\t./test_compact_diff.py [options]\n"
       print "\t-n <nr>, --nr <nr>        Specify the number of grid points\n"
       print "\t-p, --plot                Plot the differences\n"
       print "\t-i <i>, --ibc <i>         Specify the boundary type\n"
       print "\t-d <d>, --dtype <d>       Specify the domain type\n"
       sys.exit(2)

    nr = 128
    ibc = 0
    dtype = 0
    plot = False

    for opt, arg in opts:
        if opt in ("-n", "--nr"):
            nr = int(arg)
        elif opt in ("-i", "--ibc"):
            ibc = int(arg)
        elif opt in ("-d", "--dtype"):
            dtype = int(arg)
        elif opt in ("-p", "--plot"):
            plot = True

    test_diff(nr, ibc, dtype, plot)

