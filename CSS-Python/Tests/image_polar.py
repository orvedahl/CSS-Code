#!/bin/env python
#
# python module to plot polar data
#
# 7-3-2014 Orvedahl R.
#

import numpy as np
import pylab
from mpl_toolkits.axes_grid import make_axes_locatable
from matplotlib import rc

# ORIGINAL VERSION USED axes_grid1
#from mpl_toolkits.axes_grid1 import make_axes_locatable

rc('text', usetex=True)

#####################################################################
# plot polar colormap
#####################################################################
def image_polar(field,radius,costheta,sintheta,r_bcz=[0.71],mini=-1,maxi=-1,\
                cont=True,mycmap='jet',cbar=True,add_c=False,add_f=0,tit=''):
    """ Basic routine for polar plot
    --------------------------------------------------------------------------
    field              : 2D array of size (ntheta, nradius)
    radius             : radius
    costheta, sintheta : ...
    r_bcz              : add a horiztonal line at this position (given in Rsun)
    mini               : min color table value 
    maxi               : max color table value 
    cont               : overplot the contours of field
    mycmap             : python colormap to use
    cbar               : plot colorbar
    add_c              : add contours of another field 'add_f'
    add_f              : additional field for contour drawing
    tit                : title
    --------------------------------------------------------------------------
    """

    r = radius/6.9599e10
    n_r=len(r)
    n_t=len(costheta)
    rtmp = r.reshape(1,n_r)
    cthtmp = costheta.reshape(n_t,1)
    sthtmp = sintheta.reshape(n_t,1)
    xr = np.dot(cthtmp,rtmp)
    yr = np.dot(sthtmp,rtmp)

    if (mini == -1):
        mini=field.min()
    if (maxi == -1):
        maxi=max(abs(mini),field.max())

    pylab.hold(True)
    im=pylab.pcolormesh(yr,xr,field,cmap=mycmap)

    # plot dotted lines at specified radii
    for rbcz in r_bcz:
        pylab.plot(rbcz*sintheta,rbcz*costheta,'k--')#,[0,1],[0,0],'k--')

    # plot equator
    pylab.plot([r[0],r[-1]], [0,0], 'k--')

    # turn background cartesion axes off
    pylab.axis('off')
    pylab.clim((mini,maxi)) # set colorbar limits for pcolormesh

    # left boundaries
    xleftt = r[0]*sintheta[0]  # top
    yleftt = r[0]*costheta[0]
    xleftb = r[0]*sintheta[-1] # bottom
    yleftb = r[0]*costheta[-1]
    # right boundaries
    xrightt = r[-1]*sintheta[0]  # top
    yrightt = r[-1]*costheta[0]
    xrightb = r[-1]*sintheta[-1] # bottom
    yrightb = r[-1]*costheta[-1]
    
    xmin = min(xleftt, xleftb)
    xmax = r[-1]
    ymax = yrightt
    ymin = yrightb
    pylab.xlim((xmin, xmax))
    pylab.ylim((ymin, ymax))
    #xlim((0,1))
    #ylim((-1,1))
    # plot inner boundary
    pylab.plot(r[0]*sintheta,r[0]*costheta,'k')
    # plot outer boundary
    pylab.plot(r[n_r-1]*sintheta,r[n_r-1]*costheta,'k')
    pylab.scatter(r[0],0,marker='x')
    pylab.scatter(r[-1],0,marker='x')
    pylab.scatter(xleftt, yleftt, marker='x')
    pylab.scatter(xleftb, yleftb, marker='x')
    pylab.scatter(xrightt, yrightt, marker='x')
    pylab.scatter(xrightb, yrightb, marker='x')
    pylab.annotate(tit,xy=(0.05,0.05))
    pylab.hold(False)

    if (cont):
        pylab.hold(True)
        levs=mini+np.linspace(0,1,10)*(maxi-mini)
        pylab.contour(yr,xr,field,colors='w',levels=levs)
        pylab.hold(False)

    if (add_c):
        pylab.hold(True)
        mmini = add_f.min()
        mmaxi  = add_f.max()
        levs=mmini+np.linspace(0,1,10)*(mmaxi-mmini)
        levs = array([1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, \
                      1e-1, 5e-1, 1.0, 2.0])*1e23
        pylab.contour(yr,xr,add_f,colors='k',levels=levs)
        pylab.hold(False)

    pylab.hold(True)
    if (cbar):
        #divider = make_axes_locatable(pylab.gca())
        #cax = divider.append_axes("right", "5%", pad="3%")
        cb = pylab.colorbar(im)#,cax=cax)
        cb.set_clim(mini, maxi) # set limits, same as above
        cb.set_label(tit)
    pylab.hold(False)

