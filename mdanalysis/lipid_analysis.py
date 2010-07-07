#!/usr/bin/python
# -*- coding: utf-8 -*-

import optparse
import numpy
from MDAnalysis import *
import numpy.linalg
from MDAnalysis import collection
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def main():
    usage = """
    """

    parser = optparse.OptionParser(usage)
    # parser.add_option("-x", dest="col_x", default=0, help="X Column [default: %default]")
    # parser.add_option("-y", dest="col_y", default=1, help="Y Column [default: %default]")
    options, args = parser.parse_args()
    if len(args) == 0:
        parser.error("No input file specified")
    
    u = Universe(args[0], args[1])
    
    data_table = { }
    
    nlipids = len(u.selectAtoms('segid DPPC and name P'))/2
    
    top_p = u.selectAtoms('segid DPPC and resid 1:36 and name P')
    bot_p = u.selectAtoms('segid DPPC and resid 37:72 and name P')
        
    # AREA PER LIPID and BILAYER WIDTH
    data_table['Area Per Lipid'] = []
    data_table['Bilayer Width'] = []
    
    for ts in u.trajectory:
        apl = (ts.dimensions[0]*ts.dimensions[1]) / nlipids
        if apl > 0:
            data_table['Area Per Lipid'].append(apl)
        top_z = (numpy.mean([ p.pos[2] for p in top_p]), scipy.stats.sem([ p.pos[2] for p in top_p]))
        bot_z = (numpy.mean([ p.pos[2] for p in bot_p]), scipy.stats.sem([ p.pos[2] for p in bot_p]))
        data_table['Bilayer Width'].append(top_z[0]-bot_z[0])
    
    for k in data_table.keys():
        print "Mean %s: %f +-%f" % (k, numpy.mean(data_table[k]), scipy.stats.sem(data_table[k]))
    
    # DEUTERIUM ORDER PARAMETERS
    tail_carbons = numpy.arange(2,17)
    skip = 5
    order_param = numpy.zeros(len(tail_carbons))
    for i, carbon in enumerate(tail_carbons):
        selection = "resname DPPC and ( name C2%d or name H%dR or name H%dS or name C3%d or name H%dX or name H%dY )" % ((carbon,)*6)
        group = u.selectAtoms(selection)
        data = u.dcd.timeseries(group, skip=skip)
        # There are two deuteriums/carbon atom position in each acyl chain
        cd = numpy.concatenate((data[1::3]-data[0::3], data[2::3]-data[0::3]), axis=0)
        del data
        cd_r = numpy.sqrt(numpy.sum(numpy.power(cd,2), axis=-1))
        # Dot product with the z axis
        cos_theta = cd[...,2]/cd_r
        S_cd = -0.5*(3.*numpy.square(cos_theta)-1)
        # Depending on how you want to treat each deuterium, you can skip the next step
        S_cd.shape = (S_cd.shape[0], S_cd.shape[1]/4, -1) #4 deuterium order parameters/lipid carbon position
        order_param[i] = numpy.average(S_cd)
        # print "%d %f" % (i+2, order_param[i])
    
    pp = PdfPages('multipage.pdf')
    plt.figure(1)
    plt.plot(data_table['Area Per Lipid'])
    plt.plot([64]*800)
    plt.title('Area Per Lipid')
    plt.xlabel('frame')
    plt.ylabel(r'Area per lipid $\AA$')
    pp.savefig()

    plt.figure(2)
    plt.title('Bilayer Width')
    plt.plot(data_table['Bilayer Width'])
    plt.plot([38.3]*800)
    plt.xlabel('frame')
    plt.ylabel(r'Bilayer width $\AA$')
    pp.savefig()
    
    plt.figure(3)
    plt.title('Deuterium Order Parameters')
    plt.plot(numpy.arange(1,16), order_param, 'bo')
    plt.ylabel(r'$S_{CD}$')
    plt.xlabel('tail carbon number')
    pp.savefig()
    
    pp.close()
    
    
    

if __name__ == '__main__':
    main()
