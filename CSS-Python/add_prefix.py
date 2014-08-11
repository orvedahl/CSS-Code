#!/usr/bin/env python
#
# add a prefix to every line in a file and return the new file
#
# Orvedahl R. 8-11-2014
#

import sys
import getopt

def add_prefix(prefix, suffix, file):

    f = open(file, 'r')

    outfile = file+".new"
    mf = open(outfile, 'w')

    for line in f:

        # strip trailing "\n"
        s = line[:-1]

        # add prefix, suffix and "\n"
        mf.write(prefix+s+suffix+"\n")

    f.close()
    mf.close()

    print "\nRead from file: "+file
    print "\tPrefix: -->"+prefix+"<--"
    print "\tSuffix: -->"+suffix+"<--"
    print "\nNew file: "+outfile+"\n"

    return


def usage():

    print "\nPython script to add a prefix to every line of a file\n"
    print "Usage:\n"
    print "\t    ./add_prefix.py --prefix=<prefix> --file=<file>\n"
    print "\t    -f <f>, --file=<f>       Specify the file to adjust\n"
    print "\t    -p <p>, --prefix=<p>     Add <p> before every line in file\n"
    print "\t    -s <s>, --suffix=<s>     Add <s> after every line in file\n"
    print "\t    -h, --help               Display help message\n"
    print 


if __name__ == "__main__":

    try:
       opts, args = getopt.getopt(sys.argv[1:], "hf:p:s:", 
                             ["file=", "prefix=", "suffix=", "help"])

    except getopt.GetoptError:
       print "\n---ERROR: Unknown Command Line Option---\n"
       usage()
       sys.exit(2)

    # defaults
    file = None
    prefix = ""
    suffix = ""

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            usage()
            sys.exit(2)
        elif opt in ("-f", "--file"):
            file = arg
        elif opt in ("-s", "--suffix"):
            suffix = arg
        elif opt in ("-p", "--prefix"):
            prefix = arg

    add_prefix(prefix, suffix, file)

