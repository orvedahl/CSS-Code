#!/usr/bin/env python
#
# plot state variables (rho, s, u, v, w) vs. radius from checkpoint files
#
# Orvedahl R. 8-12-2014
#

import sys
import getopt
import matplotlib
import pylab
from read_checkpoint import *


def plot_vars(iter, case, eps, dpi, show, output, tex):

    if (tex):
        matplotlib.rc('text', usetex=True)
        vrl = r"$v_r$"
        vtl = r"$v_\theta$"
        vpl = r"$v_\phi$"
        rhol = r"$\rho$"
        sl   = r"$s$"
    else:
        matplotlib.rc('text', usetex=False)
        vrl = r"v-r"
        vtl = r"v-th"
        vpl = r"v-phi"
        rhol = r"rho"
        sl   = r"s"

    title1 = r"Velocity Components vs Radius"
    title2 = r"Density and Entropy vs Radius"

    # get data
    data, radius, theta, phi, time, header, ierr = read_checkpoint(iter, case)
    if (ierr):
        print "\n---ERROR: could not read checkpoint data---"
        print "\theader file: "+header+"\n"
        sys.exit(2)

    # extract the data & average over theta and phi
    #    rho = index 0 (density)
    #    w   = index 1 (v-phi)
    #    v   = index 2 (v-theta)
    #    u   = index 3 (v-r)
    #    s   = index 4 (entropy)
    rho     = numpy.mean(numpy.mean(data[:,:,:,0], axis=0), axis=0)
    vphi    = numpy.mean(numpy.mean(data[:,:,:,1], axis=0), axis=0)
    vth     = numpy.mean(numpy.mean(data[:,:,:,2], axis=0), axis=0)
    vr      = numpy.mean(numpy.mean(data[:,:,:,3], axis=0), axis=0)
    entropy = numpy.mean(numpy.mean(data[:,:,:,4], axis=0), axis=0)

    # plot it
    pylab.clf()

    # first the velocities
    pylab.subplot(211)

    pylab.plot(radius, vr,   label=vrl, color='r')
    pylab.plot(radius, vth,  label=vtl, color='b')
    pylab.plot(radius, vphi, label=vpl, color='g')

    pylab.xlabel("Radius (cm)")
    pylab.ylabel("Velocity (cm/s)")
    pylab.title(title1)

    pylab.legend(loc='upper right')
    leg = pylab.gca().get_legend()
    leg.draw_frame(False)

    fact = 0.05
    dx = radius[-1] - radius[0]
    xmax = radius[-1]
    xmin = radius[0]
    pylab.xlim(xmin-fact*dx, xmax+fact*dx)

    # next the density and entropy
    sp = pylab.subplot(212)

    pylab.title(title2)
    pylab.plot(radius, rho, label=rhol, color='r')

    if (tex):
        pylab.ylabel("Density (g cm$^{-3}$)")
    else:
        pylab.ylabel("Density (g/cm^3)")

    pylab.legend(loc='lower left')
    leg = pylab.gca().get_legend()
    leg.draw_frame(False)

    sp2 = pylab.twinx()

    pylab.plot(radius, entropy,  label=sl, color='b')
    if (tex):
        pylab.ylabel("Entropy (erg K$^{-1}$)")
    else:
        pylab.ylabel("Entropy (erg/K)")

    pylab.xlabel("Radius (cm)")
    fact = 0.05
    dx = radius[-1] - radius[0]
    xmax = radius[-1]
    xmin = radius[0]
    pylab.xlim(xmin-fact*dx, xmax+fact*dx)

    pylab.legend(loc='upper right')
    leg = pylab.gca().get_legend()
    leg.draw_frame(False)

    pylab.subplots_adjust(hspace=0.5)

    if (show):
        pylab.show()
    else:
        if (eps):
            fmt = 'eps'
            output += '.eps'
        else:
            fmt = 'png'
            output += '.png'

        pylab.savefig(output, dpi=dpi, format=fmt)
        print "\nSaved Image: "+output

    print "\n---Complete---\n"

    return


def usage():

    print "\nScript to plot state variables as a function of radius"
    print "\nUsage:\n"
    print "\tplot_variables.py [options]\n"
    print "\t-i <i>, --iter=<i>      Use iteration <i>, i.e. which"
    print "\t                           checkpoint file to use\n"
    print "\t-c <c>, --case=<c>      Which run to use\n"
    print "\t-o <o>, --output=<o>    Set output file basename to <o>\n"
    print "\t-d <d>, --dpi=<d>       Set dpi to <d> for PNG images\n"
    print "\t-e, --eps               Generate EPS images\n"
    print "\t-s, --show              Display the plot\n"
    print "\t--no-tex                Disable TeX features\n"
    print "\t-h, --help              Display help message\n"


if __name__=="__main__":

    try:
       opts, args = getopt.getopt(sys.argv[1:], "hi:c:o:d:se", 
                                  ["iter=", "output=", "eps", "dpi=",
                                   "case=", "show", "no-tex", "help"])

    except getopt.GetoptError:
       print "\n---ERROR: Unknown Command Line Option---\n"
       usage()
       sys.exit(2)

    # defaults
    eps = False
    dpi = 100
    output = ""
    iter = None
    case = ""
    show = False
    tex = True

    # parse options
    for opt, arg in opts:

        if opt in ("-h", "--help"):
            usage()
            sys.exit(2)
        elif opt in ("-e", "--eps"):
            eps = True
        elif opt in ("-d", "--dpi"):
            dpi = float(arg)
        elif opt in ("-o", "--output"):
            output = arg
        elif opt in ("-s", "--show"):
            show = True
        elif opt in ("-i", "--iter"):
            iter = int(arg)
        elif opt in ("-c", "--case"):
            case = arg
        elif opt in ("--no-tex"):
            tex = False

    plot_vars(iter, case, eps, dpi, show, output, tex)


