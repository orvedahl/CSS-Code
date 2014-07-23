#!/usr/bin/env python
#
# python module for linear algebra
#
# 7-8-2014 Orvedahl R.
#

import numpy

def tri_solver(lowr, diag, uppr, rhs, n, its=0):

    if (its == 3):

        bp = numpy.zeros(n)
        cp = numpy.zeros(n-1)
        dp = numpy.zeros(n)
        result = numpy.zeros(n)

        bp[0] = diag[0]
        cp[0] = uppr[0]
        dp[0] = rhs[0]

        for i in range(1,n-2+1):
            cp[i] = uppr[i]*bp[i-1]
            bp[i] = diag[i]*bp[i-1]-cp[i-1]*lowr[i]
            dp[i] = rhs[i]*diag[i-1]-dp[i-1]*lowr[i]

        bp[n-1] = diag[n-1]*bp[n-2]-cp[n-2]*lowr[n-1]
        dp[n-1] = rhs[n-1]*diag[n-2]-dp[n-2]*lowr[n-1]

        result[n-1] = dp[n-1]/bp[n-1]
        for i in range(n-2, 2-1, -1):
            result[i] = (dp[i]-cp[i]*result[i+1])/bp[i]

        # zero pivot
        result[1] = 5.5*rhs[2] - 5.5*result[2] - result[3]
        result[0] = rhs[0] = 10.*result[1]

    elif (its == 2):

        bp = numpy.zeros(n)
        cp = numpy.zeros(n)
        dp = numpy.zeros(n)
        result = numpy.zeros(n)

        dp[0] = rhs[0]/diag[0]
        bp[0] = 1.0

        bet = diag[1]
        dp[1] = rhs[1]/bet
        cp[1] = uppr[0]/bet
        bp[1] = lowr[1]/bet
        for i in range(2, n-1+1):
            cp[i] = uppr[i-1]/bet
            bet = diag[i]-lowr[i]*cp[i]
            dp[i] = (rhs[i]-lowr[i]*dp[i-1])/bet
            bp[i] = -bp[i-1]*lowr[i]/bet

        for i in range(n-2,0-1,-1):
            dp[i] = dp[i]-cp[i+1]*dp[i+1]
            bp[i] = bp[i]-cp[i+1]*bp[i+1]

        result[0] = dp[0]/bp[0]
        for i in range(1,n-1+1):
            result[i] = dp[i] - result[0]*bp[i]

    elif (its == 1):

        bp = numpy.zeros(n)
        cp = numpy.zeros(n)
        dp = numpy.zeros(n)
        result = numpy.zeros(n)

        dp[n-1] = rhs[n-1]/diag[n-1]
        bp[n-1] = 1.0
        bet = diag[n-2]
        dp[n-2] = rhs[n-2]/bet
        cp[n-2] = uppr[n-1]/bet
        bp[n-2] = lowr[n-2]/bet
        for i in range(n-3, 0-1, -1):
            cp[i] = uppr[i+1]/bet
            bet = diag[i]-lowr[i]*cp[i]
            dp[i] = (rhs[i]-lowr[i]*dp[i+1])/bet
            bp[i] = -bp[i+1]*lowr[i]/bet

        for i in range(1,n-1+1):
            dp[i] = dp[i] - cp[i-1]*dp[i-1]
            bp[i] = bp[i] - cp[i-1]*bp[i-1]

        result[n-1] = dp[n-1]/bp[n-1]
        for i in range(n-2,0-1,-1):
            result[i] = dp[i]-result[n-1]*bp[i]

    else:

        cp = numpy.zeros(n-1)
        dp = numpy.zeros(n)
        result = numpy.zeros(n)

        cp[0] = uppr[0]/diag[0]
        dp[0] = rhs[0]/diag[0]

        for i in range(1,n-2+1):
            cp[i] = uppr[i]/(diag[i]-cp[i-1]*lowr[i])
            dp[i] = (rhs[i]-dp[i-1]*lowr[i])/(diag[i]-cp[i-1]*lowr[i])

        dp[n-1] = (rhs[n-1]-dp[n-2]*lowr[n-1])/(diag[n-1]-cp[n-2]*lowr[n-1])

        result[n-1] = dp[n-1]
        for i in range(n-2, 0-1, -1):
            result[i] = dp[i]-cp[i]*result[i+1]

    return result


