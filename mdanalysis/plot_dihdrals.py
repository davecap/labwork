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
        
    # analysis.add_to_sequence('/protein/rmsd/backbone', RMSD(selection='backbone'))
    # analysis.add_timeseries('/protein/dihedrals/PEPA_139_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 N", "atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG")), pp=(lambda x: x*180./pi))
    # analysis.add_timeseries('/protein/dihedrals/PEPA_139_CHI2', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG", "atom PEPA 139 OD1")), pp=(lambda x: x*180./pi))
    # analysis.add_timeseries('/protein/dihedrals/PEPA_132_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 132 N", "atom PEPA 132 CA", "atom PEPA 132 CB", "atom PEPA 132 CG")), pp=(lambda x: x*180./pi))
    # analysis.add_timeseries('/protein/dihedrals/PEPA_132_CHI2', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 132 CA", "atom PEPA 132 CB", "atom PEPA 132 CG", "atom PEPA 132 OD1")), pp=(lambda x: x*180./pi))
    # analysis.add_timeseries('/protein/dihedrals/PEPA_286_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 286 N", "atom PEPA 286 CA", "atom PEPA 286 CB", "atom PEPA 286 CG")), pp=(lambda x: x*180./pi))
    # 
    # analysis.add_timeseries('/protein/distances/CUA1_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom J 1 CU", "atom J 2 CU")))
    # analysis.add_timeseries('/protein/distances/HDH102NE2_HDH102FE', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 102 NE2", "atom PEPA 102 FE")))
    # analysis.add_timeseries('/protein/distances/HDH419NE2_HDH419FE', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 419 NE2", "atom PEPA 419 FE")))
    # analysis.add_timeseries('/protein/distances/CUB_HDH419FE', Timeseries.Distance("r", trj.selectAtoms("atom I 3 CU", "atom PEPA 419 FE")))
    # analysis.add_timeseries('/protein/distances/HSD411NE2_MG', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 411 NE2", "atom I 1 MG")))
    # analysis.add_timeseries('/protein/distances/ASP412OD2_MG', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 412 OD2", "atom I 1 MG")))
    # analysis.add_timeseries('/protein/distances/HSD333NE2_CUB', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 333 NE2", "atom I 3 CU")))
    # analysis.add_timeseries('/protein/distances/HSE284ND1_CUB', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 284 ND1", "atom I 3 CU")))
    # analysis.add_timeseries('/protein/distances/HSD334NE2_CUB', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 334 NE2", "atom I 3 CU")))
    # analysis.add_timeseries('/protein/distances/CYS252SG_CUA1', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 252 SG", "atom J 1 CU")))
    # analysis.add_timeseries('/protein/distances/CYS252SG_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 252 SG", "atom J 2 CU")))
    # analysis.add_timeseries('/protein/distances/CYS256SG_CUA1', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 256 SG", "atom J 1 CU")))
    # analysis.add_timeseries('/protein/distances/CYS256SG_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 256 SG", "atom J 2 CU")))
    # analysis.add_timeseries('/protein/distances/HSE217ND1_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 217 ND1", "atom J 2 CU")))
    # analysis.add_timeseries('/protein/distances/HSE260ND1_CUA1', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 260 ND1", "atom J 1 CU")))
    # analysis.add_timeseries('/protein/distances/MET263SD_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 263 SD", "atom J 2 CU")))
    # analysis.add_timeseries('/protein/distances/GLU254O_CUA1', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 254 O", "atom J 1 CU")))
    # analysis.add_timeseries('/protein/distances/E254OE1_MG', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 254 OE1", "atom I 1 MG")))
    # 
    
    import matplotlib.pyplot as plt
    
    h5f = tables.openFile(h5_file, mode="r")
    tbl = h5f.getNode('/protein/dihedrals')
    data_chi1 = tbl.read(field='PEPA_139_CHI1')
    data_chi2 = tbl.read(field='PEPA_139_CHI2')

    # data in radians
    data_chi1_deg = []
    data_chi2_deg = []
    
    rad2deg = (lambda x: x*1)
    # rad2deg = (lambda x: x*180./pi)
    
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