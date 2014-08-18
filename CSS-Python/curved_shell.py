#!/bin/env python
#
# python module to plot 2D theta-phi shell-slice data with curved plot 
# boundaries to show the curvature
#
# Orvedahl R. 8-18-2014
#

import numpy
import pylab

# CURRENT pyASH VERSION USES matplotlib.basemap.Basemap

#####################################################################
# plot curved 2D colormap
#####################################################################
def curved_image(fig, data, equator=True, vmin=None, 
                 vmax=None, cmap="jet", cbar=True, aspect=1, extent=None, 
                 cb_title='', proj='mollweide'):

    # fig                : matplotlib figure instance on which to plot data
    # data               : 2D array of size (ntheta, nphi)
    # equator            : plot dotted line at equator
    # vmin, vmax         : min, max colorbar values
    # cmap               : matplotlib colormap instance to use
    # cbar               : include a colorbar
    # extent             : (rmin, rmax, theta_min, theta_max)
    # cb_title           : title for colorbar
    # proj               : projection to use ('mollweide', 'aitoff', 
    #                       'hammer', 'lamber', 'polar', 'rectilinear')

    # set defaults
    if (vmin == None):
        vmin = numpy.amin(data)
    if (vmax == None):
        vmax = numpy.amax(data)

    # set projection
    ax = fig.add_subplot(1,1,1, projection=proj)

    # plot equator
    if (equator):
        pylab.plot([r[0],r[-1]], [0,0], 'k--')

    if (extent == None):
        extent = None
        #fact = 0.02
        #pylab.xlim(xmin-fact*(xmax-xmin), xmax+fact*(xmax-xmin))
        #pylab.ylim(ymin-fact*(ymax-ymin), ymax+fact*(ymax-ymin))
    else:
        pylab.xlim(extent[0], extent[1])
        if (extent[2] > extent[3]):
            pylab.ylim(extent[3], extent[2])
        else:
            pylab.ylim(extent[2], extent[3])

    # use imshow
    im = ax.imshow(data, vmin=vmin, vmax=vmax, extent=extent, cmap=cmap,
                      origin='upper', clip_on=False, aspect=aspect, 
                      interpolation='quadric')

    # plot colorbar
    if (cbar):
        cb = pylab.colorbar(im)
        cb.set_clim(vmin, vmax)
        cb.set_label(cb_title)

    return


