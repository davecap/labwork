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

from tables import *


# Datastore Attributes
#   TODO: parse NAMD config file on each analysis

#   dt
#   dcdtime
#   pdb
#   psf

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
    
    analysis = Analysis('some_file.h5')
    
    analysis.add_timeseries('/protein/dihedrals/PEPA_139', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 N", "atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG")), pp=(lambda x: x*180./pi))
    analysis.add_timeseries('/protein/dihedrals/PEPA_132', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 132 N", "atom PEPA 132 CA", "atom PEPA 132 CB", "atom PEPA 132 CG")), pp=(lambda x: x*180./pi))
    analysis.add_timeseries('/protein/dihedrals/PEPA_286', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 286 N", "atom PEPA 286 CA", "atom PEPA 286 CB", "atom PEPA 286 CG")), pp=(lambda x: x*180./pi))
    analysis.add_timeseries('/protein/distances/CUA1_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom J 1 CU", "atom J 2 CU")))
    analysis.add_to_sequence('/protein/rmsd/backbone', RMSD(ref, trj, selection='backbone'))
    
    analysis.setup()
    analysis.run()
    analysis.save()
    analysis.close()
    
    # tbl_metadata = analysis.table('metadata','trajectory')
    # tbl_metadata.add_metadata(...)
    
    # tbl = analysis.table('protein','dihedrals')
    # # here we will add timeseries analyses
    # #   all timeseries types are Int32Col()
    # tbl.add_timeseries('PEPA_139', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 N", "atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG")), pp=(lambda x: x*180./pi))
    # tbl.add_timeseries('PEPA_132', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 132 N", "atom PEPA 132 CA", "atom PEPA 132 CB", "atom PEPA 132 CG")), pp=(lambda x: x*180./pi))
    # tbl.add_timeseries('PEPA_286', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 286 N", "atom PEPA 286 CA", "atom PEPA 286 CB", "atom PEPA 286 CG")), pp=(lambda x: x*180./pi))
    # 
    # tbl2 = analysis.table('protein', 'distances')
    #     tbl2.add_timeseries('CUA1_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom J 1 CU", "atom J 2 CU")))
    #     
    tbl3 = analysis.table('protein', 'rmsd')
    # here we add a sequential analysis
    #   all sequence types are Int32Col()
    tbl3.add_to_sequence('backbone', RMSD(ref, trj, selection='backbone'))
    
    # run all the analyses
    analysis.run()
    # save the data to the database
    analysis.save()
    #analysis.summary()
    analysis.close()
        
    for a in analyses:
        a.prepare(ref, trj)
    frames = trj.trajectory
    for ts in frames:
        for a in analyses:
            a.process(ts)
    for a in analyses:
        a.write()
        
    timeseries = (
                    # (name, selection, post_processor, )
                    (   'PEPA_139',
                        'dihedrals',
                        Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 N", "atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG")),
                        (lambda x: x*180./pi),
                    ),
                    (   'PEPA_132',
                        'dihedrals',
                        Timeseries.Dihedral(trj.selectAtoms("atom PEPA 132 N", "atom PEPA 132 CA", "atom PEPA 132 CB", "atom PEPA 132 CG")),
                        (lambda x: x*180./pi),
                    ),
                    (   'PEPA_286',
                        'dihedrals',
                        Timeseries.Dihedral(trj.selectAtoms("atom PEPA 286 N", "atom PEPA 286 CA", "atom PEPA 286 CB", "atom PEPA 286 CG")),
                        (lambda x: x*180./pi),
                    ),
                    (   'CUA1_CUA2',
                        'distances',
                        Timeseries.Distance("r", trj.selectAtoms("atom J 1 CU", "atom J 2 CU")),
                        None,
                    ),
                )
    
    # open the datastore
    for ts in timeseries:
        # validate datastore against all analyses to be done
        # create dataset if necessary
        collection.addTimeseries(ts[1])
        
    #data = universe.dcd.correl(collection, stop=5)
    collection.compute(trj.dcd)
    
    for i, ts in enumerate(timeseries):
        # write to the datastore
        print ts[0]
        if ts[2]:
            print ts[2](collection[i][0])
        else:
            print collection[i][0]
            
    

if __name__ == '__main__':
    main()
