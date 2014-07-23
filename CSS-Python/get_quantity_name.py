#!/usr/bin/env python
#
# return quantity names
#
# Orvedahl R. 6-16-2014
#

def get_quantity_name(quantities):

    # define names as dictionary
    names = {}

    names[1 ] = 'Rho'
    names[2 ] = 'Temperature'
    names[3 ] = 'Entropy'
    names[4 ] = 'Pressure'
    names[5 ] = 'Vr'
    names[6 ] = 'Vtheta'
    names[7 ] = 'Vphi'
    names[8 ] = 'KE'
    names[9 ] = 'ds/dr'
    names[10] = 'Br'
    names[11] = 'Btheta'
    names[12] = 'Bphi'
    names[13] = 'ME'
    names[14] = 'Poloidal Mag'
    names[15] = 'F_ks'
    names[16] = 'F_en'
    names[17] = 'F_ke'
    names[18] = 'F_ac'
    names[19] = 'F_sl'

    # if quantities is a single integer, it cannot be iterated over
    if (type(quantities) == int):
        return f(quantities, names)

    # use list to store names
    q_names = []
    for q in quantities:
        q_names.append(f(q, names))

    return q_names


def f(q, names):

    if (q not in range(1,len(names)+1)):
        name = "Unknown"
    else:
        name = names[q]

    return name

