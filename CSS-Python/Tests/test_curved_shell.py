
import sys
import getopt
import pylab
import numpy
from curved_shell import *

err = False
proj = 'ortho'
small = False
plot_original = False
try:
    opts, args = getopt.getopt(sys.argv[1:], "p:os", ["proj","orig","small"])
except getopt.GetoptError:
    print "\n---ERROR: Unknown Command Line Option---\n"
    err = True

if (not err):
    for opt, arg in opts:
        if opt in ("-p", "--proj"):
           proj = arg
        elif opt in ("-s", "--small"):
           small = True
        elif opt in ("-o", "--orig"):
           plot_original = True

if small:
   phi_d   = numpy.arange(0., 20.+0.1, 0.1)
   theta_d = numpy.arange(80., 100.+0.1, 0.1)
else:
   phi_d   = numpy.arange(-40, 40.+0.5, 0.5)
   theta_d = numpy.arange(50., 130.+0.5, 0.5)

theta = theta_d*numpy.pi/180.
phi   = phi_d*numpy.pi/180.

a = numpy.zeros((len(theta), len(phi)))
for i in range(len(theta)):
    th = theta[i]
    for j in range(len(phi)):
        p = phi[j]

        # should be positive in top-right & bottom-left corners
        # and negative in top-left & bottom-right corners and zero elsewhere
        # (this assumes theta includes 90 deg and phi includes 0 deg)
        #a[i,j] = numpy.cos(th)*numpy.sin(p)

        # should be bulls-eye with peak of 1.0 centered on th=90, phi=0
        # should also be elongated with largest edge in theta direction
        # (this assumes theta includes 90 deg and phi includes 0 deg)
        a[i,j] = numpy.sin(th)*numpy.cos(2.*p)

fig = pylab.figure(figsize=(7,7), dpi=100)

cb_title = "Hi"
aspect = 'auto'
cmap = pylab.get_cmap("jet")

vmin=None
vmax=None

extent=None

if (plot_original):
    axes = fig.add_subplot(111)

    extent=(numpy.min(phi_d), numpy.max(phi_d), 
            numpy.max(theta_d), numpy.min(theta_d))

    im = axes.imshow(a, aspect=aspect, cmap=cmap, origin="upper", 
                      extent=extent, interpolation='quadric', 
                      vmin=vmin, vmax=vmax)
    cb = fig.colorbar(im)
    cb.set_clim(vmin, vmax)
    cb.set_label=cb_title

    pylab.title("Hello")
    pylab.xlabel("phi")
    pylab.ylabel("theta")
    pylab.show()

else:
    nphi=len(phi)
    nth = len(theta)
    extent = [phi[nphi/4], phi[3*nphi/4], 
              0.5*numpy.pi-theta[nth/4], 0.5*numpy.pi - theta[3*nth/4]]
    extent = None
    lat0 = 0.5*numpy.pi - theta[len(theta)/2]
    lon0 = phi[len(phi)/2]
    curved_image(a, phi, theta, radians=True, grid=True, 
                 vmin=vmin, vmax=vmax,
                 cmap=cmap, cbar=True, cb_title=cb_title, 
                 extent=extent, proj=proj, lat_0=lat0, lon_0=lon0, 
                 title='Hello', xlabel='long', ylabel='lat')

    pylab.show()

