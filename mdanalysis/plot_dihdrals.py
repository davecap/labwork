#!/usr/bin/python
# -*- coding: utf-8 -*-

import optparse
import numpy
import numpy.linalg
import scipy.stats
from MDAnalysis import *
from MDAnalysis import collection, SelectionError
import MDAnalysis.core.rms_fitting
from MDAnalysis.core.AtomGroup import Residue, AtomGroup
from math import pi

from analysis import *

import tables

def main():
    usage = """
        usage: %prog [options] <H5 File>
    """
    parser = optparse.OptionParser(usage)
    
    options, args = parser.parse_args()
    
    if len(args) < 1:
        parser.error("No input file specified")
    
    h5_file = args[0]
    
    import matplotlib.pyplot as plt
    
    h5f = tables.openFile(h5_file, mode="r")
    tbl = h5f.getNode('/protein/dihedrals')
    data_chi1 = tbl.read(field='PEPA_139_CHI1')
    data_chi2 = tbl.read(field='PEPA_139_CHI2')

    # data in radians
    data_chi1_deg = []
    data_chi2_deg = []
    
    #rad2deg = (lambda x: x*1)
    rad2deg = (lambda x: x*180./pi)
    
    for rad in zip(data_chi1, data_chi2):
        deg = rad2deg(rad[0])
        if deg < 0:
            deg += 360
        data_chi1_deg.append(deg)
        
        deg = rad2deg(rad[1])
        if deg < 0:
            deg += 360
        data_chi2_deg.append(deg)

    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(data_chi1_deg, data_chi2_deg, c=range(len(data_chi1_deg)), cmap=plt.cm.Blues, alpha=0.75)
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 360)
    ax.set_xlabel(r'PEPA_139_CHI1')
    ax.set_ylabel(r'PEPA_139_CHI2')
    ax.set_title('139 Chi1 vs Chi2 ' + '(' + h5_file + ')')
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