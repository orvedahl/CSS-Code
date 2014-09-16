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
import integrate

import pylab

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

    # option 3:
    #                   r2
    #          tau = int   1./ vr * dr ~= (r2-r1)/vr[nr/2]
    #                   r1
    #           do integration using midpoint rule (vr[nr/2])
    #           trap rule has issues with divide by 0

    # extract Vr from data, it is stored in checkpoint.u
    data, radius, theta, phi, time, header, ierr = read_checkpoint(\
                                                   iter, case, vel_curl=True)

    if (ierr):
        print "\n---ERROR: could not read checkpoint data---"
        print "\theader file: "+header+"\n"
        sys.exit(2)

    # data is size (nth, nphi, nr, 3) where last indices are w(0), v(1), u(2)
    # with w being phi, v being theta and u being radial velocities
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

    # rms of 1./vr
    inv_vr = 1./vr
    ind_0 = numpy.where(vr==0.)
    inv_vr[ind_0] = 0.
    vrinv_rms = rms(inv_vr)

    # find up-flows and down-flows
    vr_up   = vr[vr > 0]
    vr_down = vr[vr < 0] 

    if (len(vr_up) > 0):
        vr_rms_u = rms(vr_up)
    else:
        vr_rms_u = 0.
    vr_rms_d = rms(vr_down)

    tau_1 = 2. * radial_domain / vr_rms

    print "\nOption 1:"
    print "\tuse RMS"
    print "\tTau (sec) : "+str(tau_1)
    print "\tTau (days): "+str(tau_1/86400.)

    if (len(vr_up) == 0):
        tau_2 = 2. * radial_domain / vr_rms_d
    elif (len(vr_down) == 0):
        tau_2 = 2. * radial_domain / vr_rms_u
    else:
        tau_2 = radial_domain * ( 1./vr_rms_u + 1./vr_rms_d )

    print "\nOption 2:"
    print "\tuse +/- RMS"
    print "\tTau (sec) : "+str(tau_2)
    print "\tTau (days): "+str(tau_2/86400.)

    # midpoint rule (trap rule has issues with 1/0)
    tau_3 = integrate.integrate_inv_vr(abs(vr), radius, method='simp')
    # tau_3 = radial_domain * abs(1./vr[nr/2])

    print "\nOption 3:"
    print "\tuse integral"
    print "\tTau (sec) : "+str(tau_3)
    print "\tTau (days): "+str(tau_3/86400.)

    tau_4 = 2. * radial_domain * vrinv_rms

    print "\nOption 4:"
    print "\tuse RMS of 1/Vr"
    print "\tTau (sec) : "+str(tau_4)
    print "\tTau (days): "+str(tau_4/86400.)

    frac_u = 0.8
    frac_d = 0.05
    rlow = radius[0] + frac_d*radial_domain
    rhi = radius[0] + frac_u*radial_domain
    rad_inner = numpy.where((rlow < radius) & (radius < rhi))
    vr_restrict = vr[rad_inner]
    vr_restrict_rms = rms(vr_restrict)

    # restricted domain = rhi - rlow
    tau_5 = 2.*(rhi - rlow) / vr_restrict_rms

    print "\nOption 5:"
    print "\tuse RMS of restricted Vr"
    print "\tlower bdry: "+str(frac_d)
    print "\tupper bdry: "+str(frac_u)
    print "\tTau (sec) : "+str(tau_5)
    print "\tTau (days): "+str(tau_5/86400.)

    tau_6 = 2.*radial_domain / numpy.amax(abs(vr))

    print "\nOption 6:"
    print "\tuse abs(max(Vr))"
    print "\tTau (sec) : "+str(tau_6)
    print "\tTau (days): "+str(tau_6/86400.)
    print

    pylab.clf()
    pylab.plot([rlow,rlow],[numpy.amin(vr),numpy.amax(vr)], 'k--')
    pylab.plot([rhi,rhi],[numpy.amin(vr),numpy.amax(vr)], 'k--')
    pylab.plot(radius, vr, label='vr')
    pylab.plot(radius, -numpy.amax(abs(vr))*numpy.ones(nr), label='max vr')
    pylab.plot(radius, -vr_rms*numpy.ones(nr), label='rms')
    pylab.plot(radius, -vr_restrict_rms*numpy.ones(nr), label='rest. rms')
    pylab.xlabel('radius')
    pylab.ylabel('velocity')
    pylab.legend()
    pylab.show()


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

