#!/usr/bin/env python
#
# estimate convective turnover time from checkpoint files
#
# Orvedahl R. 8-11-2014
#

import sys
import getopt
import numpy
from read_checkpoint import *

def turnover_time(iter, case):

    # option 1:
    #          tau = 2*L / V_rms_radial
    #            L = size of radial domain
    #            V_rms_radial = RMS of radial velocity

    # option 2:
    #          tau = L*( 1 / V_rms_up + 1 / V_rms_down)
    #            L = size of radial domain
    #            V_rms_up = RMS of radial velocity in upward direction
    #            V_rms_down = RMS of radial velocity in downward direction

    # extract Vr from data, it is stored in checkpoint.u
    data, radius, theta, phi, time, header, ierr = read_checkpoint(\
                                                   iter, case, vel_curl=True)

    if (ierr):
        print "\n---ERROR: could not read checkpoint data---"
        print "\theader file: "+header+"\n"
        sys.exit(2)

    # data is size (nth, nphi, nr, 3) where last indices are w(0), v(1), u(0)
    # with w being theta, v being phi and u being radial velocities
    full_vr = data[:,:,:,2]

    # size of domain
    radial_domain = radius[-1] - radius[0]

    nr = len(radius)
    nphi = len(phi)
    nth = len(theta)

    # sum over theta and phi
    vr = numpy.mean(numpy.mean(full_vr, axis=0), axis=0)

    # calculate the RMS velocity
    vr_rms = rms(vr)

    tau_1 = 2. * radial_domain / vr_rms
    print "\nOption 1:"
    print "\tuse RMS"
    print "\tTau (sec) : "+str(tau_1)
    print "\tTau (days): "+str(tau_1/86400.)

    print "\nOption 2:"
    print "\tuse +/- RMS"
    print "\tTau (sec) : "+str(tau_1)
    print "\tTau (days): "+str(tau_1/86400.)


def rms(x, axis=None):

    y = numpy.sqrt(numpy.mean(x**2, axis=axis))

    return y


def usage():

    print "\nPython script to estimate convective turnover time of CSS data\n"
    print "Usage:\n"
    print "\t    ./turnover_time.py [options]\n"
    print "\t    -i <i>, --iter=<i>       Use iteration <i>, i.e. which"
    print "\t                                checkpoint file to use\n"
    print "\t    -c <c>, --case=<c>       Which run to use\n"
    print "\t    -h, --help               Display help message\n"
    print 


if __name__ == "__main__":

    try:
       opts, args = getopt.getopt(sys.argv[1:], "hc:i:", 
                             ["case=", "iter=" "help"])

    except getopt.GetoptError:
       print "\n---ERROR: Unknown Command Line Option---\n"
       usage()
       sys.exit(2)

    # defaults
    iter = None
    case = ""

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            usage()
            sys.exit(2)
        elif opt in ("-i", "--iter"):
            iter = int(arg)
        elif opt in ("-c", "--case"):
            case = arg

    turnover_time(iter, case)

