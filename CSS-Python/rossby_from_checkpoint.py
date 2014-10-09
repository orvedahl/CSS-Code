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
import pylab

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
        avgdata = shell_avg.shell_avg(data, theta, phi)
        vr   = avgdata[:, 2]
        vth  = avgdata[:, 1]
        vphi = avgdata[:, 0]

        # magnitude of convective velocity (so Vr)
        U = numpy.sqrt(vr*vr) # + vth*vth + vphi*vphi)

        nr = len(radius)
        rossby = numpy.zeros((nr))

        radial_domain = radius[-1] - radius[0]

        # get Omega(r)
        omega = numpy.zeros((nr))
        omega[:] = omega0 + vphi[:]/radius[:]

        # calculate Rossby number, use Omega_0 instead of Omega(r)
        rossby[:] = U / (2. * radial_domain * omega0)

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
        avgdata = shell_avg.shell_avg(data, theta, phi)
        vr   = avgdata[:, 3]
        vth  = avgdata[:, 2]
        vphi = avgdata[:, 1]
        rho  = avgdata[:, 0]

        # magnitude of convective velocity (so Vr)
        U = numpy.sqrt(vr*vr) # + vth*vth + vphi*vphi)

        # get Omega(r)
        omega = numpy.zeros((nr))
        omega[:] = omega0 + vphi[:]/radius[:]

        # get scale height: 1/H = -1/rho * drho/dr
        dri = 1./(radius[-1] - radius[0])
        drhodr = derivatives.compact_fd6(dri, rho, 0., 0., 1, dtype=0)
        H = -rho / drhodr

        # calculate Rossby number, use Omega_0 instead of Omega(r)
        rossby[:] = U / (2. * H * omega0)


    elif (method == 'ens'):
        # Use the Enstrophy
        #      Ro = |curl(v)| / Omega

        # get curl with dimensions (nth, nphi, nr, 3)
        curl, radius, theta, phi = vectorcurl.vectorcurl(iter, case)

        nr = len(radius)
        rossby = numpy.zeros((nr))

        # average ens over theta, phi where ens = |curl|
        ens = numpy.sqrt(curl[:,:,:,0]*curl[:,:,:,0] +
                         curl[:,:,:,1]*curl[:,:,:,1] +
                         curl[:,:,:,2]*curl[:,:,:,2])
        ens = shell_avg.shell_avg(ens, theta, phi, numq=1)

        # ens has dimensions (nr, nq=1) so reshape it to 1D array
        ens = numpy.reshape(ens, (nr))

        #for ir in range(nr):
        #    rossby[ir] = numpy.average(ens[:,:,ir])

        # divide by rotation rate of the frame
        rossby[:] = ens[:]/(2.*omega0)
        #rossby[:] = rossby[:]/omega0

    elif (method == 'turn'):
        # Use the local convective turnover time
        #      Ro = P_rot / tau_c
        #      P_rot = 2*pi/Omega0 -- rotation period
        #      tau_c = alpha*Hp/v -- local convective turnover time
        #      -1/Hp = d(log p)/dr where p is pressure (pressure scale height)
        #      alpha = between 1 and 2, comes from MLT
        #      v = convective velocity

        alpha = 1.4 # standard value from Mixing Length Theory

        # get Cv and gamma
        gamma = 1.525     # from input file
        Cp = 3.5e8        # from input file
        Cv = Cp/gamma     # by definition

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
        avgdata = shell_avg.shell_avg(data, theta, phi)
        vr   = avgdata[:, 3]
        vth  = avgdata[:, 2]
        vphi = avgdata[:, 1]
        rho  = avgdata[:, 0]
        s    = avgdata[:, 4]

        # magnitude of convective velocity (so Vr)
        U = numpy.sqrt(vr*vr) # + vth*vth + vphi*vphi)

        # get Omega(r)
        omega = numpy.zeros((nr))
        omega[:] = omega0 + vphi[:]/radius[:]

        # get scale height: 1/H = -1/p * dp/dr
        dri = 1./(radius[1] - radius[0])
        logp = gamma*numpy.log(rho) + s/Cv
        dlogpdr = derivatives.compact_fd6(dri, logp, 0., 0., 1, dtype=0)
        H = - 1.0 / dlogpdr

        # get convective turnover time
        tau_c = numpy.zeros((nr))
        tau_c = alpha*H/U

        pylab.clf()
        print 'logp',numpy.mean(logp)
        print 'H',numpy.mean(H)
        print 'U',numpy.mean(U)
        print 't',numpy.mean(tau_c)
        print 'O',numpy.mean(omega), omega0
        print 'Ro = 6/t/O'
        pylab.plot(radius, vr, color='r', label='vr')
        pylab.plot(radius, vth, color='g', label='vt')
        pylab.plot(radius, vphi, color='b', label='vp')
        pylab.plot(radius, dlogpdr, color='r', label='dlogpdr')
        pylab.plot(radius, H, color='g', label='H')
        pylab.plot(radius, U, color='b', label='U')
        pylab.plot(radius, tau_c, color='k', label='tau')
        pylab.yscale('log')

        # calculate Rossby number, using Omega_0 instead of Omega(r)
        rossby[:] = 2.*numpy.pi / (tau_c * omega0)
        pylab.plot(radius, rossby, color='k', label='Ro')
        pylab.legend()
        pylab.show()

    else:
        print "\n---ERROR: unknown method in rossby_from_chk\n"
        sys.exit(2)

    print "\n\tAverage Rossby: "+str(numpy.mean(rossby))
    print   "\tMedian Rossby : "+str(numpy.median(rossby))+"\n"
    if (return_rad):
        return rossby, radius
    else:
        return rossby


