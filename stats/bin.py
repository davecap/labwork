#!/usr/bin/python
# -*- coding: utf-8 -*-

import optparse
from numpy import *

def main():
    usage = """
    usage: %prog <input file>

    Input file has the following format:
x1 y1
x2 y2
.. ..
xN yN
    """
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) == 0:
        parser.error("No input file specified")


    f = open('processed_diheds')
    bins = zeros((360,360))
    for line in f:
        sl = line.rstrip().split(' ')
        bins[sl[0], int(float(sl[1]))] += 1
    print "#X Y count log(count)"
    for i, row in enumerate(bins):
        for j, col in enumerate(row):
            print "%d %d %d %s" % (i, j, col, log(col))

if __name__ == '__main__':
    main()
