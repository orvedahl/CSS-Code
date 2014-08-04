#!/bin/env python
#
# python module to plot 2D r-theta polar data with curved plot boundaries 
#  in the theta direction to show the curvature
#
# Orvedahl R. 7-22-2014
#

import numpy
import pylab
from mpl_toolkits.axes_grid import make_axes_locatable
from matplotlib import rc

# ORIGINAL VERSION USED axes_grid1
# from mpl_toolkits.axes_grid1 import make_axes_locatable

#####################################################################
# plot polar 2D colormap
#####################################################################
def polar_image(data, r, theta, radians=True, r_bcz=[], equator=True, 
           vmin=None, vmax=None, cont=False, cmap=None, cbar=True, 
           aspect=1, add_cont=True, add_data=None, extent=None, cb_title=''):

    # data               : 2D array of size (ntheta, nradius)
    # r                  : radius data
    # theta              : theta data
    # radians            : is theta in radians
    # equator            : plot dotted line at equator
    # r_bcz              : add horizontal line(s) at this position
    # vmin, vmax         : min, max colorbar values
    # cont               : plot contours of data
    # cmap               : matplotlib colormap instance to use
    # cbar               : include a colorbar
    # add_cont           : add contours of different data, add_data
    # add_data           : extra data for contour drawing
    # extent             : (rmin, rmax, theta_min, theta_max)
    # cb_title           : title for colorbar

    if (not radians):
        theta = theta * numpy.pi / 180.

    costheta = numpy.cos(theta)
    sintheta = numpy.sin(theta)

    # set defaults
    if (vmin == None):
        vmin = numpy.amin(data)
    if (vmax == None):
        vmax = numpy.amax(data)

    if (cmap == None):
        cmap = pylab.get_cmap("jet")

    nth = len(costheta)
    nr = len(r)

    # set up grid
    rtmp = numpy.reshape(r, (1,nr))
    cthtmp = numpy.reshape(costheta, (nth,1))
    sthtmp = numpy.reshape(sintheta, (nth,1))
    xr = numpy.dot(cthtmp, rtmp)
    yr = numpy.dot(sthtmp, rtmp)

    # use pcolormesh
    im = pylab.pcolormesh(yr, xr, data, cmap=cmap)

    # plot horizontal lines
    for rbcz in r_bcz:
        pylab.plot(rbcz*sintheta, rbcz*costheta, 'k--')

    # plot equator
    if (equator):
        pylab.plot([r[0],r[-1]], [0,0], 'k--')

    # turn off background cartesian axis & set pcolormesh limits
    # and set aspect ratio
    pylab.axes().set_aspect(aspect)
    pylab.clim(vmin, vmax)
    pylab.axis('off')


    # define the corners of the image
    # left boundaries
    xleftt = r[0]*sintheta[0] # top
    yleftt = r[0]*costheta[0]
    xleftb = r[0]*sintheta[-1] # bottom
    yleftb = r[0]*costheta[-1]
    # right boundaries
    xrightt = r[-1]*sintheta[0] # top
    yrightt = r[-1]*costheta[0]
    xrightb = r[-1]*sintheta[-1] # bottom
    yrightb = r[-1]*costheta[-1]

    if (extent == None):
        xmin = min(xleftt, xleftb)
        xmax = numpy.amax(r[-1]*sintheta)
        ymin = min(yleftb, yrightb)
        ymax = max(yleftt, yrightt)
        pylab.xlim(xmin, xmax)
        pylab.ylim(ymin, ymax)
    else:
        pylab.xlim(extent[0], extent[1])
        if (extent[2] > extent[3]):
            pylab.ylim(extent[3], extent[2])
        else:
            pylab.ylim(extent[2], extent[3])

    # plot boundaries
    pylab.plot(r[0]*sintheta, r[0]*costheta, 'k') # inner
    pylab.plot(r[-1]*sintheta, r[-1]*costheta, 'k') # outer
    pylab.plot([xleftt, xrightt], [yleftt, yrightt], 'k') # top
    pylab.plot([xleftb, xrightb], [yleftb, yrightb], 'k') # bottom

    # plot contours
    if (cont):
        levs = vmin+numpy.linspace(0,1,10)*(vmax-vmin)
        pylab.contour(yr, xr, data, colors='w', levels=levs)

    # plot extra data contours
    if (add_cont):
        mini = numpy.amin(add_data)
        maxi = numpy.amax(add_data)
        levs = mini+numpy.linspace(0,1,10)*(maxi-mini)
        pylab.contour(yr, xr, add_data, colors='k', levels=levs)

    # plot colorbar
    if (cbar):
        cb = pylab.colorbar(im)
        cb.set_clim(vmin, vmax)
        cb.set_label(cb_title)

    return


