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
    #pylab.axis('tight')

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
        fact = 0.02
        pylab.xlim(xmin-fact*(xmax-xmin), xmax+fact*(xmax-xmin))
        pylab.ylim(ymin-fact*(ymax-ymin), ymax+fact*(ymax-ymin))
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

    #===================================================================
    # plot axes --> big ol' cludge
    #===============================
    # coords of axes, the fact + xx_width is so it lies beyond the plot
    fact1r = 0.03 # distance between plot and theta axis
    fact1t = 0.03 # distance between plot and radial axis
    factr = 0.1 # distance between theta axis and theta labels
    factt = 0.09 # distance between radial axis and radial labels
    theta_width = theta[-1] - theta[0]
    radial_width = r[-1] - r[0]

    # coords for radial axis
    theta_raxis = theta[-1] + fact1t*theta_width
    theta_rlabs = theta[-1] + factt*theta_width
    raxisxl = r[0]*numpy.sin(theta_raxis)
    raxisyl = r[0]*numpy.cos(theta_raxis)
    raxisxr = r[-1]*numpy.sin(theta_raxis)
    raxisyr = r[-1]*numpy.cos(theta_raxis)

    # coords for theta axis
    radius_taxis = r[0] - fact1r*radial_width
    radius_tlabs = r[0] - factr*radial_width

    pylab.plot([raxisxl, raxisxr], [raxisyl, raxisyr], 'k')
    pylab.plot(radius_taxis*sintheta, radius_taxis*costheta, 'k')

    # plot tick marks at lower, midpoint, upper values
    th_axis_tick_width = 0.008*radial_width
    r_axis_tick_width = 0.008*theta_width

    # radial tick marks
    sintheta_ax = numpy.sin(numpy.arange(theta_raxis-r_axis_tick_width,
                       theta_raxis+r_axis_tick_width, 0.1*r_axis_tick_width))
    costheta_ax = numpy.cos(numpy.arange(theta_raxis-r_axis_tick_width,
                       theta_raxis+r_axis_tick_width, 0.1*r_axis_tick_width))
    rlo = r[0]
    rhi = r[-1]
    rhalf = rlo + 0.5*(rhi-rlo)
    pylab.plot(rlo*sintheta_ax, rlo*costheta_ax, 'k')
    pylab.plot(rhalf*sintheta_ax, rhalf*costheta_ax, 'k')
    pylab.plot(rhi*sintheta_ax, rhi*costheta_ax, 'k')

    # radial labels
    eps = 0.994
    x1 = eps*rlo*numpy.sin(theta_rlabs)
    y1 = eps*rlo*numpy.cos(theta_rlabs)
    x2 = eps*rhalf*numpy.sin(theta_rlabs)
    y2 = eps*rhalf*numpy.cos(theta_rlabs)
    x3 = eps*rhi*numpy.sin(theta_rlabs)
    y3 = eps*rhi*numpy.cos(theta_rlabs)
    s1 = '%.2e' % (rlo)
    s2 = '%.2e' % (rhalf)
    s3 = '%.2e' % (rhi)
    pylab.annotate(s1, xy=(x1,y1), xycoords='data')
    pylab.annotate(s2, xy=(x2,y2), xycoords='data')
    pylab.annotate(s3, xy=(x3,y3), xycoords='data')

    # theta tick marks
    tlo = theta[0]
    thi = theta[-1]
    thalf = tlo + 0.5*(thi-tlo)
    rlo = radius_taxis-th_axis_tick_width
    rhi = radius_taxis+th_axis_tick_width
    lox = [rlo*numpy.sin(tlo), rhi*numpy.sin(tlo)]
    loy = [rlo*numpy.cos(tlo), rhi*numpy.cos(tlo)]
    hix = [rlo*numpy.sin(thi), rhi*numpy.sin(thi)]
    hiy = [rlo*numpy.cos(thi), rhi*numpy.cos(thi)]
    halfx = [rlo*numpy.sin(thalf), rhi*numpy.sin(thalf)]
    halfy = [rlo*numpy.cos(thalf), rhi*numpy.cos(thalf)]
    pylab.plot(lox, loy, 'k')
    pylab.plot(halfx, halfy, 'k')
    pylab.plot(hix, hiy, 'k')

    # theta labels
    eps = 1.003
    x1 = radius_tlabs*numpy.sin(tlo*eps)
    y1 = radius_tlabs*numpy.cos(tlo*eps)
    x2 = radius_tlabs*numpy.sin(thalf*eps)
    y2 = radius_tlabs*numpy.cos(thalf*eps)
    x3 = radius_tlabs*numpy.sin(thi*eps)
    y3 = radius_tlabs*numpy.cos(thi*eps)
    s1 = '%d' % (tlo*180./numpy.pi)
    s2 = '%d' % (thalf*180./numpy.pi)
    s3 = '%d' % (thi*180./numpy.pi)
    pylab.annotate(s1, xy=(x1,y1), xycoords='data')
    pylab.annotate(s2, xy=(x2,y2), xycoords='data')
    pylab.annotate(s3, xy=(x3,y3), xycoords='data')

    # axes labels
    pylab.annotate("Theta (deg)", xy=(0.15,0.5), xycoords='figure fraction',
                   verticalalignment='center', rotation='vertical')
    pylab.annotate("Radius (cm)", xy=(0.4,0.07), xycoords='figure fraction')
    #===================================================================

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


