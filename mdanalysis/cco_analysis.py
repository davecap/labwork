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

from analysis import *
from rmsd import RMSD

#
# JSON data store
#

# 'NAMD':
#   'dt' -> 0.00
#   'dcdtime' -> 0
#   'firsttimestep' -> 0
# 'PDB'
#   'file' -> "asdf.pdb"
# 'PSF'
#   'file' 
# 'TRJ'
#   'files' -> ['asdf.dcd',]
#   'frames' -> (0,1,2,3,4)
#   'times' -> (10,100,...)
# 'RMSD'
#   '<selector i>': [<rmsd values>]
# 'DIHEDRALS'
#   '<selector i>': [<dihedral values>]
# 'DISTANCES'
#   ...
#

def main():
    usage = """
        usage: %prog [options] <PSF file> <PDB file> <DCD file>
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-s", dest="selection", default="name CA", help="Atom selection [default: %default]")    
    parser.add_option("-t", dest="first_timestep", default="0", help="Starting timestep [default: %default]")    
    parser.add_option("-f", dest="dcdtime", default="1", help="DCD output frequency [default: %default]")    
    parser.add_option("-d", dest="dt", default="0.002", help="Integration Timestep [default: %default]")    
    
    options, args = parser.parse_args()
    
    if len(args) < 3:
        parser.error("No input files specified")
    
    psf_file = args[0]
    pdb_file = args[1]
    dcd_file = args[2]
    dcdtime = int(options.dcdtime)
    dt = float(options.dt)
    first_timestep = int(options.first_timestep)
    
    if not psf_file.lower().endswith('.psf'):
        print "Warning: PSF filename does not end in .psf: %s" % psf_file
    if not pdb_file.lower().endswith('.pdb'):
        print "Warning: PDB filename does not end in .pdb: %s" % pdb_file
    if not dcd_file.lower().endswith('.dcd'):
        print "Warning: DCD filename does not end in .dcd: %s" % dcd_file
    
    print "Loading reference system: %s, %s" % (psf_file, pdb_file)
    ref = Universe(psf_file, pdb_file)
    print "Loading trajectory: %s" % (dcd_file)
    trj = Universe(psf_file, dcd_file)
    
    num_frames = trj.trajectory.numframes
    # we have first_timestep, dcdtime and dt
    # frame i is first_timestep + dcdtime*i for actual timestep
    # frame i is (first_timestep + dcdtime*i)*dt for time
    frame_range = numpy.array(range(1,num_frames+1))*dcdtime+first_timestep
    time_range = frame_range*dt
    
    analyses = [RMSD]
    analyses = [ a() for a in analyses ]
    
    for a in analyses:
        a.prepare(ref, trj)
        
    frames = trj.trajectory
    for ts in frames:
        for a in analyses:
            a.process(ts)
    
    print a.results()
    
    #
    # RMSD of protein and per residue 
    #
    
    # rmsds[0] is for backbone rmsd per frame
    # rmsds[1...N] is for individual residues corresponding to residues[i+1]
    # rmsds = rmsd.rmsd_trj(trj, ref, select=zip(ref_residues, trj_residues))
    
    #
    # DIHEDRALS
    #
    
    # dihedral of residue PEPA 139
    
    #  selection of the atoms involved for the psi for resid '%d' %res
    #chi1_sel = trj.selectAtoms("atom %d N"%res, "atom %d CA"%res, "atom %d C"%res, "atom %d N" % (res+1))

    # definition collecting the timeseries of a dihedral
    # a = collection.addTimeseries(Timeseries.Dihedral(phi_sel))
    #     b = collection.addTimeseries(Timeseries.Dihedral(psi_sel))
    # 
    #     # collection of the timeseries data for every 10 steps in the traj
    #     data_phi = trj.dcd.correl(a, skip=10)*180./pi
    #     data_psi = trj.dcd.correl(b, skip=10)*180./pi
    # 
    #     # finding the avg and stdev for each residue
    #     avg_phi = mean(data_phi)
    #     stdev_phi = std(data_phi)
    #     avg_psi = mean(data_psi)
    #     stdev_psi = std(data_psi)
    #     phi.append([res,avg_phi,stdev_phi])
    #     psi.append([res,avg_psi,stdev_psi])
    # 
    #     # making an array for phi and psi data
    #     phi = numpy.array(phi)
    #     psi = numpy.array(psi)
    # 
    #     # plotting and saving the dihe for each resid
    #     res = range(2, numresidues-1)
    #     a = errorbar(phi[:,0], phi[:,1], phi[:,2], fmt='ro', label="phi")
    #     b = errorbar(psi[:,0], psi[:,1], psi[:,2], fmt='bs', label="psi")
    #     legend((a[0], b[0]), ("phi", "psi"))
    #     savefig("backbone_dihedrals.png")
    
    

if __name__ == '__main__':
    main()
