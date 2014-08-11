#!/usr/bin/env python
#
# plot shell slices for each file in a list
#
# Orvedahl R. 8-8-2014
#

import sys
import getopt
import numpy
import pylab

from read_shell_slice import *

def plot_many_shell_slice(variable, base, eps, dpi, vmin, vmax, filelist, 
                          nrad, show_labels, show_cbar, from_file):


    if (len(filelist) < 1):
        print "\n---ERROR: must specify files---\n"
        print usage()
        sys.exit(2)

    # read files from file
    if (from_file):
        file = filelist[0]
        filelist = []
        mf = open(file, 'r')
        for line in mf:
            fil = line.replace("\n", "")
            filelist.append(fil)
        mf.close()

    # set title, labels, ranges, colormaps, etc
    title = "Shell Slice of "+variable
    xlabel = "Phi (deg)"
    ylabel = "Theta (deg)"
    cmap = "jet"

    if (eps):
        fmt = 'eps'
        suffix = '.eps'
    else:
        fmt = 'png'
        suffix = '.png'

    # read first file for quantities
    data, nrecs, theta, phi, radii, quants, quant_names, ierr = \
                       read_shell_slice(filelist[0], rThetaPhiQuantities=True)
    if (ierr):
        print "\n---ERROR: read failed on first file: "+fil
        sys.exit(2)

    if (nrad >= len(radii)):
        print "\n---ERROR: desired radial slice > # radial slices---"
        print "\tDesired index: "+str(nrad)
        print "\tMax allowed  : "+str(len(radii)-1)
        sys.exit(2)

    # standard spherical coords: smallest theta is toward top of plot (ymax)
    tmin = numpy.amax(theta)*180./numpy.pi
    tmax = numpy.amin(theta)*180./numpy.pi
    pmin = numpy.amin(phi)*180./numpy.pi
    pmax = numpy.amax(phi)*180./numpy.pi

    # free the memory
    data = None
    theta = None
    phi = None
    radii = None

    # store quantities
    nq = len(quant_names)
    quantities = {}
    for i in range(nq):
        quantities[i] = quant_names[i]
    quant_names = None
    quants = None

    # get index of desired quantity
    iq = None
    for id, qnt in quantities.iteritems():
        if (qnt == variable): 
            iq = id
            break
    if (iq == None):
        print "\n---ERROR: quantity not found---"
        print "\tDesired value: "+variable
        print "\tPossible quantities:"
        for i in range(nq):
            print "\t\t"+str(quantities[i])
        sys.exit(2)

    mint = 1.e30
    maxt = -1.e30

    get_vmin = False
    get_vmax = False
    if (vmin == None):
        get_vmin = True
    if (vmax == None):
        get_vmax = True

    # loop over files
    i = 0
    for fil in filelist:

        print "\n\tBegin "+fil
        # read file & extract data
        data, nrecs, theta, phi, radii, quants, quant_names, ierr = \
                            read_shell_slice(fil, rThetaPhiQuantities=True)

        if (not ierr):
            # pick chosen radial slice, variable and first record
            data = data[:,:,nrad,iq,0]

            extent = (pmin, pmax, tmin, tmax)

            if (get_vmin):
                vmin = numpy.amin(data)
            if (get_vmax):
                vmax = numpy.amax(data)

            pylab.clf()
            cax = pylab.imshow(data, interpolation='quadric', cmap=cmap,
                               extent=extent, origin='upper', vmin=vmin,
                               vmax=vmax, aspect='auto')
            if (show_cbar):
                cb = pylab.colorbar(cax)
                cb.set_clim(vmin, vmax)

            if (show_labels):
                pylab.title(title)
                pylab.xlabel(xlabel)
                pylab.ylabel(ylabel)
            else:
                pylab.axis('off')

            # save image
            out = "_%d" % (i)
            out = base + out + suffix
            pylab.savefig(out, format=fmt, dpi=dpi)
            print "\tSaved "+out
            i += 1

            # find min vmin and max vmax over all files
            vmint = numpy.amin(data)
            vmaxt = numpy.amax(data)
            if (vmint < mint):
                mint = vmint
            if (vmaxt > maxt):
                maxt = vmaxt

        else:
            print "\n---WARNING: read failed: "+fil
            continue

        # free the memory
        data = None
        theta = None
        phi = None
        radii = None
        quants = None
        quant_names = None

    print "\nAbsolute minimum colorbar: "+str(mint)
    print "Absolute maximum colorbar: "+str(maxt)

    print "\n---Complete---\n"

    return


def usage():

    print "\nScript to plot shell slice images from many different files"
    print "\nUsage:\n"
    print "\tplot_many_shell_slice.py [options] <list of files>\n"
    print "\t-v <v>, --var=<v>       What variable to plot\n"
    print "\t-o <o>, --output=<o>    Set output file basename to <o>\n"
    print "\t-e, --eps               Generate EPS images\n"
    print "\t-d <d>, --dpi=<d>       Set dpi to <d> for PNG images\n"
    print "\t-r <r>, --rad-index=<r> Use radial slice index <r> (def=-1)\n"
    print "\t--vmin=<vmin>           Set min colorbar value for all images\n"
    print "\t--vmax=<vmax>           Set max colorbar value for all images\n"
    print "\t--no-show-labels        Do not display the titles and xy labels\n"
    print "\t--no-show-cbar          Do not display the colorbar\n"
    print "\t--from-file             Read filenames from single file given in"
    print "\t                            <list of files>\n"
    print "\t-h, --help              Display help message\n"


if __name__=="__main__":

    try:
       opts, args = getopt.getopt(sys.argv[1:], "hv:o:ed:r:", 
                                  ["var=", "output=", "eps", "dpi=",
                                   "rad-index=", "vmin=", "vmax=", 
                                   "no-show-labels", "no-show-cbar", 
                                   "from-file", "help"])

    except getopt.GetoptError:
       print "\n---ERROR: Unknown Command Line Option---\n"
       usage()
       sys.exit(2)

    # defaults
    variable = "Vr"
    base = variable
    eps = False
    dpi = 100
    nrad = -1  # default is last radial slice i.e. surface
    vmin = None
    vmax = None
    show_cbar = True
    show_labels = True
    from_file = False

    # parse options
    for opt, arg in opts:

        if opt in ("-h", "--help"):
            usage()
            sys.exit(2)
        elif opt in ("-v", "--var"):
            variable = arg
        elif opt in ("--vmin"):
            vmin = float(arg)
        elif opt in ("--vmax"):
            vmax = float(arg)
        elif opt in ("-e", "--eps"):
            eps = True
        elif opt in ("-d", "--dpi"):
            dpi = float(arg)
        elif opt in ("-r", "--rad-index"):
            nrad = int(arg)
        elif opt in ("-o", "--output"):
            base = arg
        elif opt in ("--no-show-cbar"):
            show_cbar = False
        elif opt in ("--no-show-labels"):
            show_labels = False
        elif opt in ("--from-file"):
            from_file = True

    filelist = args[:]

    plot_many_shell_slice(variable, base, eps, dpi, vmin, vmax, filelist, 
                          nrad, show_labels, show_cbar, from_file)


