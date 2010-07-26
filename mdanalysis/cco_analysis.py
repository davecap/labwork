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
    
    # analyses = [ RMSD() ]
    # for a in analyses:
    #     a.prepare(ref, trj)
    # frames = trj.trajectory
    # for ts in frames:
    #     for a in analyses:
    #         a.process(ts)
    ## rmsds[0] is for backbone rmsd per frame
    ## rmsds[1...N] is for individual residues corresponding to residues[i+1]
    ## rmsds = rmsd.rmsd_trj(trj, ref, select=zip(ref_residues, trj_residues))
    # rmsds = analyses[0].results()
    
    timeseries = (
                    # (name, selection, post_processor, )
                    (   'DIHE_PEPA_139', 
                        Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 N", "atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG")),
                        (lambda x: x*180./pi),
                    ),
                    (   'DIHE_PEPA_132',
                        Timeseries.Dihedral(trj.selectAtoms("atom PEPA 132 N", "atom PEPA 132 CA", "atom PEPA 132 CB", "atom PEPA 132 CG")),
                        (lambda x: x*180./pi),
                    ),
                    (   'DIHE_PEPA_286',
                        Timeseries.Dihedral(trj.selectAtoms("atom PEPA 286 N", "atom PEPA 286 CA", "atom PEPA 286 CB", "atom PEPA 286 CG")),
                        (lambda x: x*180./pi),
                    ),
                    (   'DIST_CUA1_CUA2',
                        Timeseries.Distance("r", trj.selectAtoms("atom J 1 CU", "atom J 2 CU")),
                        None,
                    ),
                )
    
    for ts in timeseries:
        collection.addTimeseries(ts[1])
        
    #data = universe.dcd.correl(collection, stop=5)
    collection.compute(trj.dcd)
    
    for i, ts in enumerate(timeseries):
        print ts[0]
        if ts[2]:
            print ts[2](collection[i][0])
        else:
            print collection[i][0]
            
    

if __name__ == '__main__':
    main()
