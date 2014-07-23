#!/usr/bin/env python
#
# python module for calculating derivatives
#
# 7-8-2014 Orvedahl R.
#

import linalg
import numpy

#######################################################################
# Compact Finite Difference Derivative (6th Order)
#
# dri = 1/(x2-x1) where x2 and x1 are the endpoints of the domain
# arr = array to be differenced
# b1 = boundary value on x1 (specify if ibc = 3 or 4)
# b2 = boundary value on x2 (specify if ibc = 2 or 4)
#
# ibc = boundary condition type:
#  0 --> darr/dr = 0
#  1 --> Dirichlet on x1 & x2
#  2 --> Neumann on x1 & Dirichlet on x2
#  3 --> Vice Versa
#  4 --> Neumann on x1 & x2
#
# dtype = domain type (typically use dtype=0)
#  0   --> full domain
#  1,2 --> boundary
#  3   --> internal
#
#######################################################################
def compact_fd6(dri, arr, b1, b2, ibc, dtype=0, darr=None):

    n1 = len(arr)
    if (darr == None):
        darr = numpy.zeros((n1))

    if (ibc == 5):
        dli = float(n1)
    else:
        dli = float(n1-1)

    c1 = dli*7./9.
    c2 = dli/36.

    if ((dtype == 1) or (dtype == 3)):
        darr[0] = c1*(arr[1]-b1[3]) + c2*(arr[2]-b1[2])
        darr[1] = c1*(arr[2]-arr[0]) + c2*(arr[3]-b1[3])
    elif ((dtype == 2) or (dtype == 3)):
        darr[n1-2] = c1*(arr[n1-1]-arr[n1-3]) + c2*(b2[0]-arr[n1-4])
        darr[n1-1] = c1*(b2[0]-arr[n1-2]) + c2*(b2[1]-arr[n1-3])

    for i in range(2,n1-3+1):
        darr[i] = c1*(arr[i+1]-arr[i-1]) + c2*(arr[i+2]-arr[i-2])

    alpha = 1./3.
    uppr = numpy.zeros((n1))
    lowr = numpy.zeros((n1))
    diag = numpy.zeros((n1))
    diag[:] = 1.0
    uppr[:] = alpha
    lowr[:] = alpha

    if (dtype < 3):

        if ((dtype == 1) or (dtype == 0)):
            c1 = -43./96*dli
            c2 = -5./6.*dli
            c3 = 9./8.*dli
            c4 = 1./6.*dli
            c5 = -1./96.*dli
            darr[n1-2] = -(c1*arr[n1-1] + c2*arr[n1-2] + c3*arr[n1-3] + \
                           c4*arr[n1-4] + c5*arr[n1-5])

            # fifth order coeff.
            c1 = -10./3.*dli
            c2 = -3.*dli
            c3 = 6.*dli
            c4 = 1./3.*dli
            if (ibc == 0):
                darr[n1-1] = 0.
            if ((ibc == 1) or (ibc == 2)):
                darr[n1-1] = -(c1*arr[n1-1] + c2*arr[n1-2] + c3*arr[n1-3] + \
                               c4*arr[n1-4])

            if (dtype == 1):
                darr[0] += -alpha*(\
                       -b1[0]+9.*b1[1]-45.*b1[2]+45.*arr[0]-9.*arr[1]+arr[2])
                darr[0] *= dli/60.

            # set up matrix
            # here is the 6th order interior value
            # here are the pentadiagonal and 5th order values for the boundary
            alpha2 = 3./4.
            gamma2 = 1./8.
            alpha1 = 6.
            beta1 = 3.

            # precondition the matrix to make it tridiagonal
            const = 1./(alpha2 - beta1*gamma2)
            up1 = const*(alpha1*alpha2-beta1)

            if ((ibc == 1) or (ibc == 2)):
                darr[n1-1] = const*(alpha2*darr[n1-1] - beta1*darr[n1-2])

            lowr[n1-2] = alpha2
            uppr[n1-2] = gamma2

            if (ibc >= 0):
                if ((ibc != 1) and (ibc != 2)):
                    lowr[n1-1] = 0.
                else:
                    lowr[n1-1] = up1

        if ((dtype == 2) or (dtype == 0)):
            c1 = -43./96*dli
            c2 = -5./6.*dli
            c3 = 9./8.*dli
            c4 = 1./6.*dli
            c5 = -1./96.*dli
            darr[1] = c1*arr[0]+c2*arr[1]+c3*arr[2]+c4*arr[3]+c5*arr[4]

            # 5th order coeff
            c1 = -10./3.*dli
            c2 = -3.*dli
            c3 = 6.*dli
            c4 = 1./3.*dli
            if (ibc == 0):
                darr[0] = 0.

            # if ibc =2,4 then dxdy(:,1) must be set by the calling routine
            if ((ibc == 1) or (ibc == 3)):
                darr[0] = c1*arr[0]+c2*arr[1]+c3*arr[2]+c4*arr[3]

            if (dtype == 2):
                darr[n1-1] += -alpha*(\
                                      -arr[n1-3]+9.*arr[n1-2]-45.*arr[n1-1]+\
                                      45.*b2[1]-9.*b2[2]+b2[3])
                darr[n1-1] *= dli/60.

            # set up the matrix
            # here is the 6th order interior value
            # here are the pentadiagonal and 5th order values for the boundary
            alpha2 = 3./4.
            gamma2 = 1./8.
            alpha1 = 6.
            beta1 = 3.

            # precondition the matrix to make it tridiagonal
            const = 1./(alpha2-beta1*gamma2)
            up1 = const*(alpha1*alpha2-beta1)
            if (ibc % 2 != 0):
                darr[0] = const*(alpha2*darr[0] - beta1*darr[1])
            if (ibc % 2 == 0):
                uppr[0] = 0.
            else:
                # 5th order boundary conditions
                uppr[0] = up1

            uppr[1] = alpha2
            lowr[1] = gamma2

        lowr[0] = 0.
        uppr[n1-1] = 0.

        darr = linalg.tri_solver(lowr, diag, uppr, darr, n1, its=0)

    else:
        lowr[0] = 0.
        uppr[n1-1] = 0.

        # dtype = 3 (internal boundaries only)
        darr[0] += -alpha*(\
                      -b1[0]-9.*b1[1]-45.*b1[2]+45.*arr[0]-9.*arr[1]+arr[2])
        darr[0] *= dli/60.

        darr[n1-1] += -alpha*(\
               -arr[n1-3]+9.*arr[n1-2]-45.*arr[n1-1]+45.*b2[1]-9.*b2[2]+b2[3])
        darr[n1-1] *= dli/60.

        darr = linalg.tri_solver(lowr, diag, uppr, darr, n1)

    return darr*dri


