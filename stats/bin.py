#!/usr/bin/python
# -*- coding: utf-8 -*-

import optparse
from numpy import *

def main():
    usage = """
    usage: %prog [options] <input file>

    Input file has the following format:
#0 1
x1 y1
x2 y2
.. ..
xN yN
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-x", dest="col_x", default=0, help="X Column [default: %default]")
    parser.add_option("-y", dest="col_y", default=1, help="Y Column [default: %default]")

    options, args = parser.parse_args()
    if len(args) == 0:
        parser.error("No input file specified")

    f = open(args[0])

    #bins = [[] for i in range(360)]
    #col_x = int(options.col_x)
    #col_y = int(options.col_y)

    bins = zeros((360,360))
    max = zeros(360)
    col_x = int(options.col_x)
    col_y = int(options.col_y)

    for line in f:
        if line.find('#') == 0:
            continue

        sl = line.rstrip().split(' ')
        bin_x = int(float(sl[col_x]))
        bin_y = int(float(sl[col_y]))

        bins[bin_x, bin_y] += 1
        if bins[bin_x, bin_y] > max[bin_x]:
            max[bin_x] = bin_y

    print "#X Y count %"
    for i, row in enumerate(bins):
        #print row
        for j, col in enumerate(row):
            print "%d %d %d %f" % (i, j, col, col/max[i]*100)

if __name__ == '__main__':
    main()
