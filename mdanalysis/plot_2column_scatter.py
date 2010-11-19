#!/usr/bin/python
# -*- coding: utf-8 -*-

import optparse
import numpy
import numpy.linalg
import scipy.stats
from math import pi

from analysis import *
import tables

import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt
import matplotlib.font_manager

def rad2deg(rad):
    deg = rad*180./pi
    if deg < 0:
        deg += 360.
    return deg

def main():
    usage = """
        usage: %prog [options] <H5 File 1> <H5 File 2> <H5 File N>
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-x", dest="x_column", default=None, help="X Column REQUIRED")
    parser.add_option("-y", dest="y_column", default=None, help="Y Column REQUIRED")
    
    options, args = parser.parse_args()
    
    if not args:
        parser.error("No input file(s) specified")
    
    if not options.x_column or not options.y_column:        
        print "Two paths are required! showing all possible paths..."
        h5f = tables.openFile(args[0], mode="r")
        for g in h5f.root._v_children.keys():
            for t in h5f.root._v_children[g]._v_children.keys():
                for c in h5f.root._v_children[g]._v_children[t].description._v_names:
                    print "/%s/%s/%s" % (g, t, c)
        h5f.close()
        parser.error("X and Y paths are required.")
    
    column1 = options.x_column.split('/')[-1]
    ps1 = options.x_column.split('/')
    ps1.pop()
    table1 = '/'.join(ps1)

    column2 = options.y_column.split('/')[-1]
    ps2 = options.y_column.split('/')
    ps2.pop()
    table2 = '/'.join(ps2)
    
    data1 = numpy.array([])
    data2 = numpy.array([])
    
    vfunc = numpy.vectorize(rad2deg)
    
    while args:
        h5_file = args.pop()
        print "Extracting data from %s" % h5_file
        h5f = tables.openFile(h5_file, mode="r")
        data1 = numpy.append(data1, vfunc(h5f.getNode(table1).read(field=column1)))
        data2 = numpy.append(data2, vfunc(h5f.getNode(table2).read(field=column2)))
        h5f.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax.scatter(data1, data2, c=range(len(data1)), cmap=plt.cm.Blues, alpha=0.75)
    ax.scatter(data1, data2)
    
    # ax.set_xlim(0, 360)
    # ax.set_ylim(0, 360)
    ax.set_xlabel(column1)
    ax.set_ylabel(column2)
    ax.set_title(column1 + ' vs ' + column2)
    ax.grid(True)
    plt.savefig(column1 + ' vs ' + column2 + '.eps')
    
    # bins = numpy.arange(-180, 180)
    # freq, bins = numpy.histogram(data_chi1, bins=bins, range=None, normed=False, weights=None, new=None)
    # print freq
    # 
    # freq, bins = numpy.histogram(data_chi2, bins=bins, range=None, normed=False, weights=None, new=None)
    # print freq
    
    # pylab.plot(.5*(bins[1:]+bins[:-1]), n)
    # pylab.show()    

if __name__ == '__main__':
    main()