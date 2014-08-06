#!/usr/bin/env python
#
# read Checkpoint files
#
# Orvedahl R. 8-5-2014
#


import sys
import string
import numpy
import defaults
import unformatted_read

#####################################################################
# read Checkpoint files
#####################################################################
def read_checkpoint(iter, case):

    ierr = False

    #--------------------------------------------------------------
    # read header file
    #--------------------------------------------------------------
    # write iteration to string and strip whitespace
    stemp = "%07d" % (iter)
    iterstr = string.strip(stemp)

    # construct filename
    read_dir = defaults.dir + case + "/Checkpoints/" + iterstr + "/"
    fname = string.strip(read_dir + "header")

    # read both lines of header
    temp1, temp2 = read_header(fname)

    # error trap
    if (len(temp1) != 5 or len(temp2) != 8):
        print "\nERROR: header has incorrect number of variables"
        print "\theader file: "+fname
        print
        sys.exit(2)

    # extract the numbers from each line
    nr   = int(temp1[0])
    nth  = int(temp1[1])
    nph  = int(temp1[2])
    nq   = int(temp1[3])
    iter = int(temp1[4])

    print temp1

    time = numpy.zeros(2)
    time[0] = temp2[0] # dt
    time[1] = temp2[1] # t
    r1      = temp2[2]
    r2      = temp2[3]
    th1     = temp2[4]
    th2     = temp2[5]
    ph1     = temp2[6]
    ph2     = temp2[7]

    print temp2

    # construct each checkpoint.x filename
    prefix =read_dir+"checkpoint."
    names = ["rho", "w", "v", "u", "s", "bt", "bp", "br"]
    files = []
    for n in names:
        files.append(prefix+n)

    # loop over each filename and add the data
    data = numpy.empty((nth, nph, nr, nq))
    for iv in range(nq):

        t_data, err = unformatted_read.read_f90(files[iv])
        if (err):
            ierr = True
            break

        # reshape the data
        t_data = numpy.reshape(t_data, (nth, nph, nr), order='F')

        # store data
        data[:,:,:,iv] = t_data[:,:,:]

    if (not ierr):
        # get radius, theta and phi data
        rad = (r2-r1)*numpy.array(range(nr))/float(nr-1) + r1
        theta = (th2-th1)*numpy.array(range(nth))/float(nth-1) + th1
        phi = (ph2-ph1)*numpy.array(range(nph))/float(nph-1) + ph1

        return data, rad, theta, phi, time, ierr

    else:

        return None, None, None, None, None, ierr


#####################################################################
# read header file
#####################################################################
def read_header(fname):

    f = open(fname, 'r')

    nline = 0
    for line in f:

        # skip comments and empty lines
        if (not(line.lstrip().startswith("#") or line.lstrip() == "")):

            nline += 1

            # remove newline and return characters
            line = line.replace("\n", "")
            line = line.replace("\t", "")

            # split line based on white space
            numbers = line.split()

            # remove all white space, change scientific notation D to E
            k = 0
            for num in numbers:
                num = num.replace(" ", "")
                num = num.replace("\t", "")
                num = num.replace("D", "E")
                num = float(num)
                numbers[k] = num
                k += 1

            # header should only have 2 lines
            if (nline == 1):
                line1 = numbers
            elif (nline == 2):
                line2 = numbers

    f.close()

    return line1, line2


