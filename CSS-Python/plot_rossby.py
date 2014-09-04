#!/usr/bin/env python
#
# plot rossby number as function of radius from checkpoint files
#
# Orvedahl R. 8-12-2014
#

import sys
import getopt
import matplotlib
import pylab
import numpy
from rossby_from_checkpoint import *


def plot_rossby(iter, case, omega0, eps, dpi, show, output, tex, method):

    if (tex):
        matplotlib.rc('text', usetex=True)
        title = r"Rossby Number vs Radius, $\Omega_0 = %.3e$" % (omega0)
    else:
        matplotlib.rc('text', usetex=False)
        title = r"Rossby Number vs Radius, Omega_0 = %.3e" % (omega0)

    # get rossby number
    rossby, radius = rossby_from_chk(iter, case, omega0, method=method, 
                                       return_rad=True)

    # plot it
    pylab.clf()

    pylab.plot(radius, rossby, label='Rossby')

    ones = numpy.ones(len(radius))
    pylab.plot(radius, ones, label='Unity', linestyle=":", color='k')

    pylab.xlabel("Radius (cm)")
    pylab.ylabel("Rossby Number")
    pylab.title(title)

    pylab.legend(loc='upper right')
    leg = pylab.gca().get_legend()
    leg.draw_frame(False)

    fact = 0.05
    dx = radius[-1] - radius[0]
    xmax = radius[-1]
    xmin = radius[0]
    pylab.xlim(xmin-fact*dx, xmax+fact*dx)

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

    print "\nScript to plot rossby number as function of radius"
    print "\nUsage:\n"
    print "\tplot_rossby.py [options] <list of files>\n"
    print "\t-i <i>, --iter=<i>      Use iteration <i>, i.e. which"
    print "\t                           checkpoint file to use\n"
    print "\t-c <c>, --case=<c>      Which run to use\n"
    print "\t--omega=<omega0>        Specify rotation rate\n"
    print "\t--method=<m>            Get Rossby using <m>: 'adv', 'ens'\n"
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
                                   "case=", "omega=", "show", "no-tex",
                                   "method=", "help"])

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
    omega0 = 2.66e-6 # solar
    show = False
    tex = True
    method = 'adv'

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
        elif opt in ("--omega"):
            omega0 = float(arg)
        elif opt in ("--no-tex"):
            tex = False
        elif opt in ("--method"):
            method = arg

    plot_rossby(iter, case, omega0, eps, dpi, show, output, tex, method)


