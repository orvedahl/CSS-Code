#!/bin/env python
#
# python module to plot 2D theta-phi shell-slice data with curved plot 
# boundaries to show the curvature
#
# Orvedahl R. 8-18-2014
#

import numpy
import pylab
from mpl_toolkits.basemap import Basemap

#####################################################################
# plot curved 2D colormap
#####################################################################
def curved_image(data, phi, theta, radians=True, grid=True, vmin=None, 
                 vmax=None, cmap="jet", cbar=True, extent=None, 
                 cb_title='', proj='ortho', title='', xlabel='Longitude (deg)',
                 ylabel='Latitude (deg)', lat_0=0, lon_0=0, nlon=5, nlat=5):

    # data               : 2D array of size (ntheta, nphi)
    # phi                : 1D array of phi values
    # theta              : 1D array of theta values
    # lon_0              : central longitude (phi direction)
    # lat_0              : central latitude (theta direction)
    # radians            : are phi, theta, lat_0 & lon_0 in radians or degrees
    # grid               : plot latitude/longitude grid lines
    # vmin, vmax         : min, max colorbar values
    # cmap               : matplotlib colormap instance to use
    # cbar               : include a colorbar
    # extent             : (phi_min, phi_max, theta_min, theta_max)
    # cb_title           : title for colorbar
    # nlon               : number of longitude lines to draw
    # nlat               : number of latitude lines to draw
    # proj               : projection to use ('moll', 'tmerc', 'merc', 'omerc'
    #                        'robin', 'ortho', ...
    # title              : give a title to the plot
    # xlabel             : give a label to the x-axis
    # ylabel             : give a label to the y-axis

    # convert everything to degrees
    if (radians):
        theta = theta * 180. / numpy.pi
        phi = phi * 180. / numpy.pi
        lat_0 = lat_0 * 180. / numpy.pi
        lon_0 = lon_0 * 180. / numpy.pi
        if (extent != None):
            extent[0] = extent[0] * 180. / numpy.pi # longitude min
            extent[1] = extent[1] * 180. / numpy.pi # longitude max
            extent[2] = extent[2] * 180. / numpy.pi # latitude min
            extent[3] = extent[3] * 180. / numpy.pi # latitude max

    nphi = len(phi)
    nth  = len(theta)

    # set defaults
    if (vmin == None):
        vmin = numpy.amin(data)
    if (vmax == None):
        vmax = numpy.amax(data)

    # set plot extent
    if (extent == None):
        lonmin = phi[0]
        lonmax = phi[-1]
        latmin = 90. - theta[-1]
        latmax = 90. - theta[0]

    else:
        lonmin = extent[0]
        lonmax = extent[1]
        latmin = extent[2]
        latmax = extent[3]

    fact = 0.0
    ll_lon = lonmin - fact*(lonmax - lonmin) # lon of lower left corner
    ll_lat = latmin - fact*(latmax - latmin) # lat of lower left corner
    ur_lon = lonmax + fact*(lonmax - lonmin) # lon of upper right corner
    ur_lat = latmax + fact*(latmax - latmin) # lat of upper right corner

    # set up main map
    map = Basemap(projection=proj, lon_0=lon_0, lat_0=lat_0)

    # set boundaries of map
    ll_x, temp = map(ll_lon, lat_0)  # lat_0 because it is the widest
    temp, ll_y = map(lon_0, ll_lat)  # lon_0 because it is the widest
    ur_x, temp = map(ur_lon, lat_0)  # lat_0 because it is the widest
    temp, ur_y = map(lon_0, ur_lat)  # lon_0 because it is the widest

    map.llcrnrx = ll_x
    map.llcrnry = ll_y
    map.urcrnrx = ur_x
    map.urcrnry = ur_y

    # set up grid
    lon_tmp = phi
    lat_tmp = 90. - theta
    x, y = map(*numpy.meshgrid(lon_tmp, lat_tmp))

    # STANDARD pcolormesh:
    # if data is data(x,y) then calling sequence is:
    #       pcolormesh(y, x, data, ...)
    # or
    #       pcolormesh(x, y, data.transpose(), ...)
    # data array has dimensions (nth, nphi)
    #
    # MAP.pcolormesh:
    # if data is data(x,y) then calling sequence is:
    #       pcolormesh(x, y, data, ...)
    im = map.pcolormesh(x, y, data, cmap=cmap, latlon=True, 
                        vmin=vmin, vmax=vmax)
    #map.drawmapboundary()

    # plotting grid lines
    if (grid):
        dlon = (lonmax-lonmin)/(nlon+1)
        dlat = (latmax-latmin)/(nlat+1)
        parallels = numpy.arange(latmin, latmax+dlon, dlon)
        meridians = numpy.arange(lonmin, lonmax+dlat, dlat)
        col = "0.1"
        ls = ":"
        for par in parallels:
            xl, yl = map(lonmin, par)
            xr, yr = map(lonmax, par)

            # line itself
            map.plot([xl,xr], [yl,yr], linestyle=ls, color=col)

            # label
            titl = str(round(par,1))
            x, y = map(lonmin-0.10*(lonmax-lonmin), par)
            pylab.text(x, y, titl, fontsize=7)

        for mer in meridians:
            latl, lonl = latmin, mer
            latr, lonr = latmax, mer

            # line itself
            map.drawgreatcircle(lonl,latl,lonr,latr,linestyle=ls,color=col)

            # label
            titl = str(round(mer,1))
            x, y = map(mer, latmin-0.03*(latmax-latmin))
            pylab.text(x, y, titl, fontsize=7)

        # OLD WAY:
        # labels=[left, right, top, bottom]
        # only label latitude lines that intersect the left
        # only label longitude lines that intersect the bottom
        #map.drawparallels(parallels, latmax=latmax, labels=[1,0,0,0])
        #map.drawmeridians(meridians, latmax=latmax, labels=[0,0,0,1])

    # plot colorbar
    if (cbar):
        cb = pylab.colorbar(im, shrink=0.9, pad=0.05)
        cb.set_clim(vmin, vmax)
        cb.set_label(cb_title)

    pylab.title(title)
    pylab.xlabel(xlabel, labelpad=18)
    pylab.ylabel(ylabel, labelpad=33)

    return


