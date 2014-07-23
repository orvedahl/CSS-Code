#!/bin/env python
#
# python module to plot polar data
#
# 7-3-2014 Orvedahl R.
#

import numpy as np
from pylab import *
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
    field              : 2D array 
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
    xr = dot(cthtmp,rtmp)
    yr = dot(sthtmp,rtmp)

    if (mini == -1):
        mini=field.min()
    if (maxi == -1):
        maxi=max(abs(mini),field.max())

    hold(True)
    im=pcolormesh(yr,xr,field,cmap=mycmap)
    for rbcz in r_bcz:
        plot(rbcz*sintheta,rbcz*costheta,'k--',[0,1],[0,0],'k--')
    #plot((4.83/6.9599)*sintheta,(4.83/6.9599)*costheta,'k--',[0,1],[0,0],'k--')
    axis('equal')
    axis('off')
    clim((mini,maxi))
    xlim((0,1))
    ylim((-1,1))
    plot(r[0]*sintheta,r[0]*costheta,'k')
    plot(r[n_r-1]*sintheta,r[n_r-1]*costheta,'k')
    plot([0,0],[-r[n_r-1],r[n_r-1]],'k--')
    plot([0,0],[-r[0],-r[n_r-1]],'k',[0,0],[r[n_r-1],r[0]],'k')
    annotate(tit,xy=(0.05,0.05))
    hold(False)

    if (cont):
        hold(True)
        levs=mini+np.linspace(0,1,10)*(maxi-mini)
        contour(yr,xr,field,colors='w',levels=levs)
        hold(False)

    if (add_c):
        hold(True)
        mmini = add_f.min()
        mmaxi  = add_f.max()
        levs=mmini+np.linspace(0,1,10)*(mmaxi-mmini)
        levs = array([1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, \
                      1e-1, 5e-1, 1.0, 2.0])*1e23
        contour(yr,xr,add_f,colors='k',levels=levs)
        hold(False)

    hold(True)
    if (cbar):
        divider = make_axes_locatable(gca())
        #cax = divider.append_axes("right", "5%", pad="3%")
        colorbar(im)#,cax=cax)
    hold (False)

