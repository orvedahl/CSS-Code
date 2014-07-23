#!/usr/bin/env python
#
# python module for reading scalar data
#
# 6-11-2014 Orvedahl R.
#

import numpy

def read_scalar(filename, mkuniq=None, mag=None, imag=None):

    # get number of lines in file
    niters = getData(filename, nlines=True)

    if (imag != None):
        niters = niters - imag

    print "Reading niters=",niters

    ltmp = 0
    if (mag != None):
        values = numpy.zeros((niters, 9), dtype=numpy.float64)
        if (imag == None):
            # read data
            data = getData(filename)

            # extract iterations and data values
            iters = data[:,0]                
            values[:,:] = data[:,1:]

        else:
            print "Skipping ", imag, " lines."
            # read data
            data = getData(filename, skip=imag)

            # extract iterations and data values
            iters = data[:,0]                
            values[:,:] = data[:,1:]
            
    else:
        values = numpy.zeros((niters, 7), dtype=numpy.float64)
        # read data
        data = getData(filename)

        # extract iterations and data values
        iters = data[:,0]                
        values[:,:] = data[:,1:]

    if (mkuniq != None):

        # remove duplicate records and sort iterations
        iters, values = remove_duplicates(iters, values)

        niters = len(iters)

        # write out fixed file
        mf = open(filename, 'w')
        if (mag != None):
            format = "%10d"+9*"%15.7E" + "\n"
            for i in range(niters):
                mf.write(format % (iters[i], \
                                    values[i,0],\
                                    values[i,1],\
                                    values[i,2],\
                                    values[i,3],\
                                    values[i,4],\
                                    values[i,5],\
                                    values[i,6],\
                                    values[i,7],\
                                    values[i,8]))

        else:
            format = "%10d"+7*"%15.7E" + "\n"
            for i in range(niters):
                mf.write(format % (iters[i], \
                                    values[i,0],\
                                    values[i,1],\
                                    values[i,2],\
                                    values[i,3],\
                                    values[i,4],\
                                    values[i,5],\
                                    values[i,6]))
        mf.close() 

    return iters, values


#####################################################################
# remove duplicate entries and sort array
#####################################################################
def remove_duplicates(x1, x2):

    # x1[k]   = iteration number
    # x2[k,j] = jth data point for iteration k

    # sort both arrays based on the iteration array
    ind_sort = numpy.argsort(x1)
    x1 = x1[ind_sort]
    x2 = x2[ind_sort,:]

    # find number of unique elements (numpy.unique only operates on 1D data)
    num_unique = len(numpy.unique(x1))
    ny = len(x2[0,:])

    # allocate new storage
    iters = numpy.zeros((num_unique))
    values = numpy.zeros((num_unique, ny))

    # store unique elements (first element is assumed to be unique)
    iters[0] = x1[0]
    values[0,:] = x2[0,:]
    j = 1
    for k in range(1,len(x1)):

        if (x1[k] == x1[k-1]):
            continue

        # store unique elements
        iters[j] = x1[k]
        values[j,:] = x2[k,:]
        j += 1

    return iters, values


#####################################################################
# get numerical data from file
#####################################################################
def getData(filename, separator=None, skip=0, nlines=False):

    # if separator is None, all whitespace will be stripped
    # otherwise separator will be used to split the line and white space
    # needs to be taken care of

    f = open(filename, "r")

    if (skip > 0):
        f = f.readlines()[skip:]

    # get number of data points
    index = 0
    for line in f:

        # skip empty lines and lines starting with "#"
        if (not(line.lstrip().startswith("#") or line.lstrip() == "")):

            # only do this for the first value line
            if index == 0:

               # get how many columns of data to read
               fields = line.split(separator)
               size = len(fields)

            # update number of valid lines
            index = index + 1

    if (skip == 0):
        f.close()

    # return number of lines
    if (nlines):
        return index

    # allocate space for data
    data = numpy.zeros((index,size), numpy.float64)

    index = 0
    f = open(filename, "r")

    if (skip > 0):
        f = f.readlines()[skip:]

    # read data into arrays
    for line in f:

        if (not(line.lstrip().startswith("#") or line.lstrip() == "")):

            # remove newline characters and return characters
            line = line.replace("\n", "")
            line = line.replace("\r", "")

            # split the line using "separator" to delimit different columns
            fields = line.split(separator)

            # remove all white space from each entry
            k = 0
            for fil in fields:
                fil = fil.replace(" ", "")
                fil = fil.replace("\t", "")
                fil = fil.replace("D", "E")
                fil = float(fil)
                fields[k] = fil
                k += 1

            # fill array
            data[index,:] = fields

            index = index + 1

    if (skip == 0):
        f.close()

    return data


