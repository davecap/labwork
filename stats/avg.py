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
    
    bins = [[] for i in range(360)]
    col_x = int(options.col_x)
    col_y = int(options.col_y)

# for each line, fill append the y_values to the correct bin
    for line in f:
        if line.find('#') == 0:
            continue
        sl = line.rstrip().split(' ')
        bin_x = int(float(sl[col_x]))
        y_val = float(sl[col_y])
        bins[bin_x].append(y_val)

    print "#X Y SD SDm SEm"
    for x_val, y_array in enumerate(bins):
        if len(y_array) < 2:
            print "%d - - - -" % (x_val)
        else:
            avg = mean(y_array)
            count = len(y_array)
# calculate the sample SD and SD
            sum_sq = 0
            for y_val in y_array:
                sum_sq += (y_val - avg) * (y_val - avg)
# Bessel's correction (N - 1)
            sample_sd = sqrt(sum_sq/(count - 1))
            sd = sqrt(sum_sq/count)

            sd_mean = sd/sqrt(count)
            se_mean = sample_sd/sqrt(count)
            print "%d %f %f %f %f" % (x_val, avg, sd, sd_mean, se_mean)

if __name__ == '__main__':
    main()
