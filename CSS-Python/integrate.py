#
# module to perform integrations
#
# R. Orvedahl 9-16-2014

import sys
import numpy

##########################################################################
# interface to 1D integration routines
##########################################################################
def integrate1D(func, x, dx, method="simp"):

    if (method == "simp"):
        integral = simpson_integrate(func, dx)

    elif (method == "trap"):
        integral = trap_integrate(func, x)

    else:
        print "\n---ERROR: unknown integration method: "+method
        print
        sys.exit(2)

    return integral


##########################################################################
# Integrate 1./v(r) over r
##########################################################################
def integrate_inv_vr(vr, rad, method='simp'):

    intgrnd = 1./vr

    # deal with possible 1./0.
    ind_0 = numpy.where(vr == 0.)
    intgrnd[ind_0] = 0.

    # simpson's method with even/odd number of slabs
    if (method == "simp"):
        dr = rad[1] - rad[0]    # assumes a uniform grid
        integral = simpson_integrate(intgrnd, dr)

    # trapezoid rule
    elif (method == "trap"):
        integral = trap_integrate(intgrnd, rad)

    return integral


##########################################################################
# Trapezoidal Rule Integration (1D)
##########################################################################
def trap_integrate(func, x):

    integral = 0.
    N = len(x) - 1

    n = 0
    while (n < N):
        #           <---width---> <---height in middle---->
        integral += (x[n+1]-x[n])*0.5*(func[n] + func[n+1])
        n += 1

    return integral


##########################################################################
# Simpson Rule Integration (1D)
##########################################################################
def simpson_integrate(func, dx):

    nslabs = len(func) - 1

    integral = 0.

    if nslabs%2==0:
       odd = False
    else:
       odd = True

    n = 0
    while (n <= nslabs - 2):
        # simpson integration over even number of slabs
        integral += (1./3.)*dx*(func[n] + 4.*func[n+1] + func[n+2])
        n += 2

    if (odd):
        # include last slab if it exists
        integral += dx/12.*(\
                      -func[nslabs-2] + 8.*func[nslabs-1] + 5.*func[nslabs])

    return integral


