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

from rmsd import rmsd_trj

def main():
    usage = """
        usage: %prog [options] <PDB file> <PSF file> <DCD file>
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-s", dest="selection", default="name CA", help="Atom selection [default: %default]")    
    
    options, args = parser.parse_args()
    
    if len(args) < 3:
        parser.error("No input files specified")
    
    psf_file = args[0]
    pdb_file = args[1]
    dcd_file = args[2]
    
    ref = Universe(psf_file, pdb_file)
    trj = Universe(psf_file, dcd_file)
    
    #
    # RMSD of Calphas 
    #
       
    #rmsds_calpha = rmsd_trj(trj, ref, select="name CA")
    
    #
    # RMSD per residue
    #
    
    # TODO
    
    #
    # DIHEDRALS
    #
    
    # dihedral of residue PEPA 139
    
    #  selection of the atoms involved for the psi for resid '%d' %res
    chi1_sel = trj.selectAtoms("atom %d N"%res, "atom %d CA"%res, "atom %d C"%res, "atom %d N" % (res+1))

    # definition collecting the timeseries of a dihedral
    a = collection.addTimeseries(Timeseries.Dihedral(phi_sel))
    b = collection.addTimeseries(Timeseries.Dihedral(psi_sel))

    # collection of the timeseries data for every 10 steps in the traj
    data_phi = trj.dcd.correl(a, skip=10)*180./pi
    data_psi = trj.dcd.correl(b, skip=10)*180./pi

    # finding the avg and stdev for each residue
    avg_phi = mean(data_phi)
    stdev_phi = std(data_phi)
    avg_psi = mean(data_psi)
    stdev_psi = std(data_psi)
    phi.append([res,avg_phi,stdev_phi])
    psi.append([res,avg_psi,stdev_psi])

    # making an array for phi and psi data
    phi = numpy.array(phi)
    psi = numpy.array(psi)

    # plotting and saving the dihe for each resid
    res = range(2, numresidues-1)
    a = errorbar(phi[:,0], phi[:,1], phi[:,2], fmt='ro', label="phi")
    b = errorbar(psi[:,0], psi[:,1], psi[:,2], fmt='bs', label="psi")
    legend((a[0], b[0]), ("phi", "psi"))
    savefig("backbone_dihedrals.png")
    
    # pp = PdfPages('multipage.pdf')
    # plt.figure(1)
    # plt.plot(data_table['Area Per Lipid'])
    # plt.plot([64]*800)
    # plt.title('Area Per Lipid')
    # plt.xlabel('frame')
    # plt.ylabel(r'Area per lipid $\AA$')
    # pp.savefig()
    # 
    # plt.figure(2)
    # plt.title('Bilayer Width')
    # plt.plot(data_table['Bilayer Width'])
    # plt.plot([38.3]*800)
    # plt.xlabel('frame')
    # plt.ylabel(r'Bilayer width $\AA$')
    # pp.savefig()
    # 
    # plt.figure(3)
    # plt.title('Deuterium Order Parameters')
    # plt.plot(numpy.arange(1,16), order_param, 'bo')
    # plt.ylabel(r'$S_{CD}$')
    # plt.xlabel('tail carbon number')
    # pp.savefig()
    # 
    # pp.close()
    
    
    

if __name__ == '__main__':
    main()
