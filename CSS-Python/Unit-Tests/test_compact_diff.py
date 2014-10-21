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
    # the dr multiplier is to make it grid invariant
    global_err = numpy.sqrt(dr*numpy.sum(local_err*local_err))

    print "\nglobal error 1st deriv: ",global_err
    print "deriv at r[nr/4]", df1[nr/4]
    print "exact at r[nr/4]", exact1[nr/4]
    print "dr:",dr
    print "nr:",len(r)

    return

if __name__ == "__main__":

    try:
       opts, args = getopt.getopt(sys.argv[1:], "n:", ["--nr"])
    except getopt.GetoptError:
       print "\nUnknown command line arg"
       print "\nUsage:\n"
       print "\t./test_compact_diff.py -n <nr>"
       print "\t./test_compact_diff.py --nr <nr>"
       print
       sys.exit(2)

    nr = 128

    for opt, arg in opts:
        if opt in ("-n", "--nr"):
            nr = int(arg)

    test_diff(nr, 0, 0)

