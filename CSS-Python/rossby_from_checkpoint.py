#!/usr/bin/env python
#
# calculate Rossby number from checkpoint files
#
# Orvedahl R. 8-7-2014
#

import sys
import numpy
import vectorcurl
from read_checkpoint import *
import derivatives
import shell_avg

def rossby_from_chk(iter, case, omega0, method='adv', return_rad=False):

    if (method == 'adv'):
        # Use the advective term:
        #      Ro =  U / (2*L*Omega) where U = magnitude of vel

        # get velocity, dont need rho or s --> vel_curl=True
        data, radius, theta, phi, time, header, ierr = read_checkpoint(\
                                                  iter, case, vel_curl=True)
        if (ierr):
            print "\n---ERROR: could not read checkpoint data---"
            print "\theader file: "+header+"\n"
            sys.exit(2)

        # average velocities over theta and phi
        # data is of size (nth, nphi, nr, 3) where last indices are
        #    0, w --> u_phi
        #    1, v --> u_theta
        #    2, u --> u_r

        # average velocity components over theta,phi
        avgdata = shell_avg.shell_avg_rms(data, theta, phi)
        vr   = avgdata[:, 2]
        vth  = avgdata[:, 1]
        vphi = avgdata[:, 0]

        # magnitude of velocity
        U = numpy.sqrt(vr*vr + vth*vth + vphi*vphi)

        nr = len(radius)
        rossby = numpy.zeros((nr))

        radial_domain = radius[-1] - radius[0]

        # get Omega(r)
        omega = numpy.zeros((nr))
        omega[:] = omega0 + vphi[:]/radius[:]

        # calculate Rossby number
        rossby[:] = U / (2. * radial_domain * omega)

    elif (method == 'scale'):
        # Use Pressure Scale Height
        #      Ro = Vrms / (2. * Omega * H)
        #      1/H = drho/dr / rho

        # get data
        # data is of size (nth, nphi, nr, 5) where last indices are
        #    0, rho --> density
        #    1, w --> u_phi
        #    2, v --> u_theta
        #    3, u --> u_r
        #    4, s --> entropy
        data, radius, theta, phi, time, header, ierr = read_checkpoint(\
                                                  iter, case, vel_curl=False)
        if (ierr):
            print "\n---ERROR: could not read checkpoint data---"
            print "\theader file: "+header+"\n"
            sys.exit(2)

        nr = len(radius)
        rossby = numpy.zeros((nr))

        # average velocity components over theta,phi
        avgdata = shell_avg.shell_avg_rms(data, theta, phi)
        vr   = avgdata[:, 3]
        vth  = avgdata[:, 2]
        vphi = avgdata[:, 1]
        rho  = avgdata[:, 0]

        # magnitude of velocity
        U = numpy.sqrt(vr*vr + vth*vth + vphi*vphi)

        # get Omega(r)
        omega = numpy.zeros((nr))
        omega[:] = omega0 + vphi[:]/radius[:]

        # get scale height: 1/H = -1/rho * drho/dr
        dri = 1./(radius[-1] - radius[0])
        drhodr = derivatives.compact_fd6(dri, rho, 0., 0., 1, dtype=0)
        H = -rho / drhodr

        # calculate Rossby number
        rossby[:] = U / (2. * H * omega)


    elif (method == 'ens'):
        # Use the Enstrophy
        #      Ro = |curl(v)| / Omega

        # get curl
        curl, radius, theta, phi = vectorcurl.vectorcurl(iter, case)

        nr = len(radius)
        rossby = numpy.zeros((nr))

        # average curl over theta, phi
        avgcurl = shell_avg.shell_avg_rms(curl, theta, phi)
        ens = numpy.sqrt(avgcurl[:,0]*avgcurl[:,0] +
                         avgcurl[:,1]*avgcurl[:,2] +
                         avgcurl[:,2]*avgcurl[:,2])
        # get magnitude of curl --> sqrt(enstrophy)
        #ens = numpy.sqrt(curl[:,:,:,0]*curl[:,:,:,0] +
        #                 curl[:,:,:,1]*curl[:,:,:,1] +
        #                 curl[:,:,:,2]*curl[:,:,:,2])

        #for ir in range(nr):
        #    rossby[ir] = numpy.average(ens[:,:,ir])

        # divide by rotation rate of the frame
        rossby[:] = ens[:]/omega0
        #rossby[:] = rossby[:]/omega0

    else:
        print "\n---ERROR: unknown method in rossby_from_chk\n"
        sys.exit(2)

    if (return_rad):
        return rossby, radius
    else:
        return rossby


