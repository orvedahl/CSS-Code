#!/usr/bin/env python
#
# read unformatted data into python that was generated in fortran as:
#	open(unit, file, form='unformatted', access='direct', recl=irec)
#
# 6-12-2014 Orvedahl R.

import FortranFile
import numpy

def main():

    filename = "unformatted.dat"

    #f = FortranFile.FortranFile(filename)#, endian='=')
    #dims = f.readInts()
    #data = f.readReals(prec='f')
    #data = f.readRecord()
    #f.close()
    #print type(dims)
    #print dims

    data = read_record(filename)

    nx = data[0]
    ny = data[1]+2
    nz = data[2]

    nx = int(nx)
    ny = int(ny)
    nz = int(nz)
    nrec = len(data) / (nx*ny*nz)
    #print nx,ny,nz,len(data),nx*ny*nz,nrec
    data = numpy.reshape(data, (nx,ny,nz,nrec), order='F')

    nth  = data[0,0,0,0]
    nr   = data[1,0,0,0]
    nq   = data[2,0,0,0]
    r1   = data[3,0,0,0]
    r2   = data[4,0,0,0]

    print "\nExpect:"
    print 384, 190, 200
    print 6.2e10, 6.96e10
    
    print "Actual:"
    print nth, nr, nq
    print r1, r2

    print "\nExpect:"
    print 384, 192, 200

    print "Actual:"
    print numpy.shape(data)
    print

    for r in range(numpy.shape(data)[-1]):

        print "---Record: "+str(r+1)
        print "a(3, 7, 2)   = " + str(data[3-1,7-1,2-1,r])
        print "a(13,72,25)  = " + str(data[13-1,72-1,25-1,r])
        print "a(1,17,12)   = " + str(data[1-1,17-1,12-1,r])
        print "a(33,1,77)   = " + str(data[33-1,1-1,77-1,r])
        print "a(39,77,1)   = " + str(data[39-1,77-1,1-1,r])
        print "a(7, 11, 5)  = " + str(data[7-1,11-1,5-1,r])
        print "a(11, 7, 5)  = " + str(data[11-1,7-1,5-1,r])
        print "a(11, 5, 7)  = " + str(data[11-1,5-1,7-1,r])
        print "a(7,  5, 11) = " + str(data[7-1,5-1,11-1,r])
        print "-----------------------------------------"

    print

    return

def read_record(filename):

    f = open(filename, 'rb')

    #f.seek(8*record*nx*ny*nz)

    #field = numpy.fromfile(f, dtype='float64', count=nx*ny*nz)
    field = numpy.fromfile(f, dtype='float64')

    #field = numpy.reshape(field, (nx,ny,nz))

    f.close()

    return field

main()

