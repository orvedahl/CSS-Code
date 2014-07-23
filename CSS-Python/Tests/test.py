#!/usr/bin/env python

import read_scalar
import numpy

def test1():
    f = 'test.dat'
    nlines = read_scalar.getData(f, nlines=True)

    print "\nRead data from: "+f
    print "There are ",nlines," lines of data total"

    data = read_scalar.getData(f, skip=4)

    print "Skipping 4 lines:"
    print data

    return

def test2():
    x = numpy.zeros((6,3))

    for i in range(6):
        x[i,1] = i
        x[i,2] = i

    iters = numpy.array([1,3,2,3,2,4])
    x[0,0]=1
    x[1,0]=3
    x[2,0]=2
    x[3,0]=3
    x[4,0]=2
    x[5,0]=4

    print "Original"
    print iters
    print x
    print "Now Sort and Remove Duplicates"
    iters, x = read_scalar.remove_duplicates(iters,x)
    print iters
    print x

    return

def test3a():

    x = numpy.arange(2,8)
    print "before call 1:",x
    test3b()
    print "after call 1:",x

    x = numpy.arange(2,8)
    print "before call 2:",x
    test3b(array=x)
    print "after call 2:",x

def test3b(array=None):

    if (array == None):
        array = numpy.arange(2,8)
        array[0] = 6
        return
    else:
        array[0] = 0

    return

