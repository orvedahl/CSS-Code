#!/usr/bin/env python
#
# Convert from (latitude, longitude) ---> Projection (x, y) coords
#
# See:
#	http://en.wikipedia.org/wiki/Mollweide_projection
#
#	http://en.wikipedia.org/wiki/Kavrayskiy_VII_projection
#
#       http://en.wikipedia.org/wiki/Orthographic_projection_in_cartography
#
# MOLLWEIDE:
# lam = longitude
# lam0 = central longitude (prime meridian)
# phi = latitude
# theta = intermediate angle
# x = cartesian x coordinate in Mollweide projection
# y = cartesian y coordinate in Mollweide projection
#
#      x = 2*sqrt(2)/pi * (lam - lam0) * cos(theta)
#      y = sqrt(2) * sin(theta)
#
# wheretheta is defined as:
#      2*theta + sin(2*theta) = pi*sin(phi)
#
# there is an inverse, given (x,y) return (lat, lon)
#      lam = lam0 + pi*x / (2*sqrt(2)*cos(theta))
#      sin(theta) = y/sqrt(2)
#      sin(phi) = (2*theta + sin(2*theta)) / pi
#
# KAVRAYSKIY VII:
# lam = longitude (radians)
# phi = latitude (radians)
# x = cartesian x coordinate in Kav 7 projection
# y = cartesian y coordinate in Kav 7 projection
#
#      x = 1.5*lam*sqrt(1/3 - (phi/pi)**2)
#      y = phi
#
# ORTHOGRAPHIC:
# lam = longitude
# phi = latitude
# R = radius of sphere
# lam0 = center longitude
# phi0 = center latitude
# x = cartesian x coordinate in Ortho projection
# y = cartesian y coordinate in Ortho projection
#
#      x = R * cos(phi) * sin(lam-lam0)
#      y = R * [cos(phi0)*sin(phi) - sin(phi0)*cos(phi)*cos(lam-lam0)]
#
#
# Orvedahl R. 8-20-2014
#

import sys
import numpy

def lat_lon_to_xy(lat, lon, proj='ortho', lon0=0, tol=1.e-10, radians=True):

    if (not radians):
        lat = lat*numpy.pi/180.
        lon = lon*numpy.pi/180.
        lon0 = lon0*numpy.pi/180.

    if (proj == 'kav7'):

        nlat = len(lat)
        nlon = len(lon)
        x = []
        for i in range(nlat):
            phi = lat[i]
            for j in range(nlon):
                lam = lon[j]
                x.append(1.5*lam*numpy.sqrt(1./3. - (phi/numpy.pi)**2))
        x = numpy.array(x)
        y = lat

    elif (proj == 'moll'):
        # find theta first
        theta = theta_newton(lat, tol)

        sth = numpy.sin(theta)
        cth = numpy.cos(theta)

        # coefficients
        Ay = numpy.sqrt(2.)
        Ax = 2. * Ay / numpy.pi

        # (x,y) coordinates
        x = Ax * (lon - lon0) * cth
        y = Ay * sth

    else:
        print "\n---ERROR: Unknown projection type: "+proj
        print "\n"
        sys.exit(2)

    return x, y


# Newton find to get theta, based on given latitude
def theta_newton(phi, tol):

    theta = numpy.zeros((len(phi)))

    i = 0
    for p in phi:

        # Newton find is analytic
        # theta_0 = phi
        #                        2.*theta_n + sin(2.*theta_n) - pi*sin(phi)
        # theta_n+1 = theta_n -  ------------------------------------------
        #                                 2. + 2.*cos(2.*theta_n)

        # initial guess for newton find
        th0 = p

        sinp = numpy.sin(p)

        num = 2*th0 + numpy.sin(2.*th0) - numpy.pi*sinp
        den = 2. + 2.*numpy.cos(2.*th0)
        dth = -num/den

        th = th0 + dth

        # handle the poles differently
        if (abs(p) == 0.5*numpy.pi):
            th = p
            dth = 0.1*tol

        # newton find
        while (abs(dth) > tol):

            num = 2*th + numpy.sin(2.*th) - numpy.pi*sinp
            den = 2. + 2.*numpy.cos(2.*th)
            dth = -num/den

            th += dth

        # store theta
        theta[i] = th
        i += 1

    return theta

