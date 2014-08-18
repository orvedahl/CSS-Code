
import pylab
import numpy
from curved_shell import *

plot_original = False

phi   = numpy.arange(-40, 40, 0.1)
theta = numpy.arange(50., 130., 0.1)
theta = theta*numpy.pi/180.
phi   = phi*numpy.pi/180.

a = numpy.zeros((len(theta), len(phi)))
for i in range(len(theta)):
    th = theta[i]
    for j in range(len(phi)):
        p = phi[j]
        a[i,j] = numpy.sin(th)*numpy.cos(p)

fig = pylab.figure(figsize=(7,7), dpi=100)

cb_title = "Hi"
aspect = 'auto'
cmap = pylab.get_cmap("jet")

vmin=None
vmax=None

extent=None

if (plot_original):
    axes = fig.add_subplot(111)

    im = axes.imshow(a, aspect=aspect, cmap=cmap, origin="upper", 
                      extent=extent, interpolation='quadric', 
                      vmin=vmin, vmax=vmax)
    cb = fig.colorbar(im)
    cb.set_clim(vmin, vmac)
    cb.set_label=cb_title

    pylab.title("Hello")
    pylab.show()

else:

    curved_image(fig, a, equator=False, vmin=vmin, vmax=vmax,
                 cmap=cmap, cbar=True, cb_title=cb_title, 
                 aspect=aspect, extent=extent, proj='aitoff')
    pylab.title("Hello")

    pylab.show()

