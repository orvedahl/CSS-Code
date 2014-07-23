#!/bin/env python
#
# python module to show matplotlib colormaps
#
# 6-24-2014 Orvedahl R.
#

from pylab import *
from numpy import outer
rc('text', usetex=False)
a=outer(arange(0,1,0.01),ones(10))
a=a.transpose()
#figure(figsize=(6, 10))
figure()
subplots_adjust(top=0.98,bottom=0.02,left=0.2,right=0.95, hspace=.2)
maps=[m for m in cm.datad if not m.endswith("_r")]
maps.sort()
l=len(maps)+1

#build rectangle in axes coords
left, width = -.02, .1
bottom, height = .0, 1.
right = left + width
top = bottom + height
for i, m in enumerate(maps):
    subplot(l,1,i+1)
    axis("off")
    imshow(a,aspect='auto',cmap=get_cmap(m),origin="lower")
    ax = gca()
    ax.text(left, 0.5*(bottom+top), m, horizontalalignment="right",
            verticalalignment="center", rotation="horizontal",
            transform=ax.transAxes, fontsize=8)
savefig("colormaps.png",dpi=300,facecolor='gray')
