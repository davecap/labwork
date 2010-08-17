#!/usr/bin/python
# -*- coding: utf-8 -*-

import optparse
import numpy
import numpy.linalg
import scipy.stats
from math import pi

from analysis import *
import tables
import matplotlib.pyplot as plt

def main():
    usage = """
        usage: %prog [options] <H5 File> <path 1> <path 2>
    """
    parser = optparse.OptionParser(usage)
    
    options, args = parser.parse_args()
    
    try:
        h5_file = args[0]
    except:
        parser.error("No input file specified")
    
    try:
        path1 = args[1]
        path2 = args[2]
    except:
        print "Two paths are required! showing all possible paths..."
        h5f = tables.openFile(h5_file, mode="r")
        for g in h5f.root._v_children.keys():
            for t in h5f.root._v_children[g]._v_children.keys():
                for c in h5f.root._v_children[g]._v_children[t].description._v_names:
                    print "/%s/%s/%s" % (g, t, c)
        h5f.close()
        parser.error("No path(s) specified")
    
    column1 = path1.split('/')[-1]
    ps1 = path1.split('/')
    ps1.pop()
    table1 = '/'.join(ps1)
    
    column2 = path2.split('/')[-1]
    ps2 = path2.split('/')
    ps2.pop()
    table2 = '/'.join(ps2)
    
    h5f = tables.openFile(h5_file, mode="r")
    data1 = h5f.getNode(table1).read(field=column1)
    data2 = h5f.getNode(table2).read(field=column2)
    
    h5f.close()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(data1, data2, c=range(len(data1)), cmap=plt.cm.Blues, alpha=0.75)
    # ax.set_xlim(0, 360)
    # ax.set_ylim(0, 360)
    ax.set_xlabel(column1)
    ax.set_ylabel(column2)
    ax.set_title(column1 + ' vs ' + column2 + ' (' + h5_file + ')')
    ax.grid(True)
    plt.show()
    
    # bins = numpy.arange(-180, 180)
    # freq, bins = numpy.histogram(data_chi1, bins=bins, range=None, normed=False, weights=None, new=None)
    # print freq
    # 
    # freq, bins = numpy.histogram(data_chi2, bins=bins, range=None, normed=False, weights=None, new=None)
    # print freq
    
    # pylab.plot(.5*(bins[1:]+bins[:-1]), n)
    # pylab.show()
    
    h5f.close()
    

if __name__ == '__main__':
    main()