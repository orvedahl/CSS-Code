#!/bin/env python
#
# python module to plot scalar data
#
# 6-12-2014 Orvedahl R.
#

import read_scalar
import defaults
import sys
import getopt
import pylab
from matplotlib import rc

def plot_scalar(css_case, eps=False, dir=None, isi=None, ie=None, log=None, 
            show=True, mkuniq=None, mag=None, imag=None, tex=True):

    if (css_case == None or css_case == ""):
        print "\n---ERROR: Must supply CSS Case---\n"
        sys.exit(2)

    if (tex):
        rc('text', usetex=True)
    else:
        rc('text', usetex=False)

    if (dir == None):
        dir = defaults.dir
    filename = dir + css_case + "/scalar.dat"

    iters, data = read_scalar.read_scalar(filename, mkuniq=mkuniq, 
                                          mag=mag, imag=imag)

    if (isi == None):
        isi = 2
    if (ie == None):
        ie = len(iters)-1

    # the extra "+1" is to make isi:ie inclusive
    times = data[isi:ie+1,0]
    dt    = data[isi:ie+1,1]
    ke    = data[isi:ie+1,2]
    cke   = data[isi:ie+1,3]
    mcke  = data[isi:ie+1,4]
    drke  = data[isi:ie+1,5]
    mach  = data[isi:ie+1,6]
    if (mag != None):
        me   = data[isi:ie+1,7]
        bmax = data[isi:ie+1,8]

        yrng = [min(min(cke),min(me)), max(ke)]
        brng = [min(bmax), max(bmax)]
    else:
        yrng = [min(cke), max(ke)]

    day = 24.0*3600.0

    pylab.clf()
    if (mag != None):
        pylab.subplot(311)
    else:
        pylab.subplot(211)
    pylab.title("Volume Averaged KE")
    pylab.xlabel("time (days)")
    if (tex):
        pylab.ylabel(r"$<KE>$")
    else:
        pylab.ylabel(r"avg KE")
    pylab.plot(times/day, ke,   label="KE",   color="r")
    pylab.plot(times/day, cke,  label="CKE",  color="g")
    pylab.plot(times/day, mcke, label="MCKE", color="b")
    pylab.plot(times/day, drke, label="DRKE", color="k")
    if (mag != None):
        pylab.plot(times/day, me, label="ME", color="m")
        pylab.plot(times/day, me+ke, label="ME+KE", color="c")

    pylab.legend(loc='upper left')
    leg = pylab.gca().get_legend()
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    leg.draw_frame(False)
    pylab.ylim(yrng)

    if (log != None):
        pylab.yscale('log')

    if (mag != None):
        pylab.subplot(312)
    else:
        pylab.subplot(212)
    pylab.title("Maximum Mach \#")
    pylab.xlabel("time (days)")
    if (tex):
        pylab.ylabel(r"$<M>$")
        pylab.title("Maximum Mach \#")
    else:
        pylab.ylabel(r"M")
        pylab.title("Maximum Mach #")
    pylab.plot(times/day, mach, color="r")
    if (log != None):
        pylab.yscale('log')

    if (mag != None):
        pylab.subplot(313)
        pylab.title("Maximum B Field")
        pylab.xlabel("time (days)")
        if (tex):
            pylab.ylabel(r"$B_{\rm{max}}$")
        else:
            pylab.ylabel(r"B max")
        pylab.plot(times/day, mach, color="r")
        pylab.ylim(byrng)
        if (log != None):
            pylab.yscale('log')

    pylab.subplots_adjust(hspace=.5)

    if (show):
        pylab.show()
    else:
        output = "scalar"
        if (eps):
            output += ".eps"
        else:
            output += ".png"
        pylab.savefig(output, pad_inches=0.33)
        print "\tSaved Image: "+output

    return

def usage():

    print "\nPython script to plot scalar quantities from CSS\n"
    print "Usage:\n"
    print "    ./plot_scalar.py --case=<case> [options]\n"
    print "    -c <case>, --case=<case> Required location of data\n"
    print "    --dir=<dir>              Use <dir> as root directory of <case>\n"
    print "    --isi=<isi>              Specify data start index\n"
    print "    --ie=<ie>                Specify data end index (inclusive)\n"
    print "    --imag=<imag>            Skip <imag> lines\n"
    print "    --mag                    Data includes magnetic fields\n"
    print "    --mkuniq                 Sort and remove duplicate records\n"
    print "    --no-tex                 Do not use TeX supported labels\n"
    print "    --log                    Plot y-axis on log scale\n"
    print "    --eps                    Generate an EPS figure\n"
    print "    --save                   Save figure\n"
    print "    -h, --help               Display help message\n"
 
 
if __name__ == "__main__":

    try:
       opts, args = getopt.getopt(sys.argv[1:], "hc:", ["case=", "no-tex", 
                              "dir=", "isi=", "ie=", "imag=", "mag", 
                              "mkuniq", "log", "eps", "save", "help"])

    except getopt.GetoptError:
       print "\n---ERROR: Unknown Command Line Option---\n"
       usage()
       sys.exit(2)

    # defaults
    case = ""
    dir = None
    isi = None
    ie = None
    imag = None
    mag = None
    mkuniq = None
    tex = True
    log = None
    eps = False
    show = True

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            usage()
            sys.exit(2)
        elif opt in ("-c", "--case"):
            case = arg
        elif opt in ("--no-tex"):
            tex = False
        elif opt in ("--dir"):
            dir = arg
        elif opt in ("--isi"):
            isi = int(arg)
        elif opt in ("--ie"):
            ie = int(arg)
        elif opt in ("--imag"):
            imag = int(arg)
        elif opt in ("--mag"):
            mag = True
        elif opt in ("--mkuniq"):
            mkuniq = True
        elif opt in ("--log"):
            log = True
        elif opt in ("--eps"):
            eps = True
        elif opt in ("--save"):
            show = False

    plot_scalar(case, eps=eps, dir=dir, isi=isi, ie=ie, log=log, show=show,
                mkuniq=mkuniq, mag=mag, imag=imag, tex=tex)

