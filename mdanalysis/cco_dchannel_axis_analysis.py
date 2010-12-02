#!/usr/bin/python
# -*- coding: utf-8 -*-

# import required modules
import optparse # command-line parsing
import numpy # numerical operations
from tables import * # pytables (database)
from MDAnalysis import * # analysis of NAMD trajectories

# my analysis modules (in the same directory as this file)
from analysis import *
from rmsd import RMSD
from hbonds import HydrogenBondAnalysis
from nearby import NearbyCountAnalysis
from cylindersearch import CylinderSearch

def main():
    usage = """
        usage: %prog [options] <PSF file> <PDB file> <DCD file 1> <DCD file 2> <DCD file n>
    """
    
    # initialize the parser with our custom usage string (above)
    parser = optparse.OptionParser(usage)
    # add all the command-line options
    parser.add_option("-o", dest="h5_filename", default="analysis.h5", help="Output analysis file [default: %default]")    
    parser.add_option("-t", dest="first_timestep", default=0, type="int", help="Starting timestep [default: %default]")    
    parser.add_option("-f", dest="dcdtime", default=1, type="int", help="DCD output frequency [default: %default]")    
    parser.add_option("-d", dest="dt", default=0.002, type="float", help="Integration Timestep [default: %default]")    
    
    # parse the command line
    (options, args) = parser.parse_args()
    
    # check to make sure all arguments are there
    if len(args) < 3:
        parser.error("No input files specified")
    
    psf_file = args[0]
    pdb_file = args[1]
    
    # do some error checking for required file paths
    if not os.path.exists(psf_file):
        print "Error: PSF file doesn't exist: %s" % psf_file
        return
    elif not psf_file.lower().endswith('.psf'):
        print "Warning: PSF filename does not end in .psf: %s" % psf_file
    
    if not os.path.exists(pdb_file):
        print "Error: PDB file doesn't exist: %s" % pdb_file
        return
    elif not pdb_file.lower().endswith('.pdb'):
        print "Warning: PDB filename does not end in .pdb: %s" % pdb_file
    
    # load the PSF and PDB into the analysis system
    print "Loading reference system: %s, %s" % (psf_file, pdb_file)
    ref = Universe(psf_file, pdb_file, permissive=True)
    
    for dcd_file in args[2:]:    
        if not os.path.exists(dcd_file):
            print "Error: DCD file doesn't exist: %s" % dcd_file
            return
        elif not dcd_file.lower().endswith('.dcd'):
            print "Warning: DCD filename does not end in .dcd: %s" % dcd_file
    
        print "Loading trajectory: %s" % (dcd_file)
        trj = Universe(psf_file, dcd_file)
    
        # calculate frames and timesteps
        num_frames = trj.trajectory.numframes
        # frame i is first_timestep + dcdtime*i for timestep
        # frame i is (first_timestep + dcdtime*i)*dt for time
        frame_range = numpy.array(range(1,num_frames+1))*options.dcdtime+options.first_timestep
        time_range = frame_range*options.dt
    
        # add the metadata for this analysis to the database
        analysis = Analysis(options.h5_filename, readonly=False)
        analysis.add_metadata('/metadata/trajectory', { 'psf': psf_file.split('/')[-1], 'pdb': pdb_file.split('/')[-1], 'dcd': dcd_file.split('/')[-1], 'frames': num_frames, 'firsttimestep': options.first_timestep, 'dt': options.dt })
    
        # RMSD
        # this takes too long so I'm disabling it for now
        #analysis.add_to_sequence('/protein/rmsd/backbone_PEPA', RMSD(selection='backbone and segname PEPA'))
    
        # coordinates of POT and CA of 132 and 139 and 286
        analysis.add_timeseries('/protein/coords/PEPA139_z', Timeseries.Atom("z", trj.selectAtoms("segid PEPA and resid 139")))
        analysis.add_timeseries('/protein/coords/PEPA139CA_x', Timeseries.Atom("x", trj.selectAtoms("atom PEPA 139 CA")))
        analysis.add_timeseries('/protein/coords/PEPA139CA_y', Timeseries.Atom("y", trj.selectAtoms("atom PEPA 139 CA")))
        analysis.add_timeseries('/protein/coords/PEPA139CA_z', Timeseries.Atom("z", trj.selectAtoms("atom PEPA 139 CA")))
        
        analysis.add_timeseries('/protein/coords/PEPA132_z', Timeseries.Atom("z", trj.selectAtoms("segid PEPA and resid 132")))
        analysis.add_timeseries('/protein/coords/PEPA132CA_x', Timeseries.Atom("x", trj.selectAtoms("atom PEPA 132 CA")))
        analysis.add_timeseries('/protein/coords/PEPA132CA_y', Timeseries.Atom("y", trj.selectAtoms("atom PEPA 132 CA")))
        analysis.add_timeseries('/protein/coords/PEPA132CA_z', Timeseries.Atom("z", trj.selectAtoms("atom PEPA 132 CA")))
        
        analysis.add_timeseries('/protein/coords/PEPA207_z', Timeseries.Atom("z", trj.selectAtoms("segid PEPA and resid 207")))
        analysis.add_timeseries('/protein/coords/PEPA207CA_x', Timeseries.Atom("x", trj.selectAtoms("atom PEPA 207 CA")))
        analysis.add_timeseries('/protein/coords/PEPA207CA_y', Timeseries.Atom("y", trj.selectAtoms("atom PEPA 207 CA")))
        analysis.add_timeseries('/protein/coords/PEPA207CA_z', Timeseries.Atom("z", trj.selectAtoms("atom PEPA 207 CA")))
        
        analysis.add_timeseries('/protein/coords/PEPA286_z', Timeseries.Atom("z", trj.selectAtoms("segid PEPA and resid 286")))
        analysis.add_timeseries('/protein/coords/PEPA286CA_x', Timeseries.Atom("x", trj.selectAtoms("atom PEPA 286 CA")))
        analysis.add_timeseries('/protein/coords/PEPA286CA_y', Timeseries.Atom("y", trj.selectAtoms("atom PEPA 286 CA")))
        analysis.add_timeseries('/protein/coords/PEPA286CA_z', Timeseries.Atom("z", trj.selectAtoms("atom PEPA 286 CA")))
        
        analysis.add_timeseries('/protein/distances/PEPA139CA_PEPA132CA', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 139 CA", "atom PEPA 132 CA")))
        analysis.add_timeseries('/protein/distances/PEPA139CA_PEPA286CA', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 139 CA", "atom PEPA 286 CA")))
        analysis.add_timeseries('/protein/distances/PEPA132CA_PEPA286CA', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 132 CA", "atom PEPA 286 CA")))

        analysis.add_timeseries('/protein/coords/PEPC7CA_x', Timeseries.Atom("x", trj.selectAtoms("atom PEPC 7 CA")))
        analysis.add_timeseries('/protein/coords/PEPC7CA_y', Timeseries.Atom("y", trj.selectAtoms("atom PEPC 7 CA")))
        analysis.add_timeseries('/protein/coords/PEPC7CA_z', Timeseries.Atom("z", trj.selectAtoms("atom PEPC 7 CA")))
        
        analysis.add_timeseries('/protein/coords/PEPC10CA_x', Timeseries.Atom("x", trj.selectAtoms("atom PEPC 10 CA")))
        analysis.add_timeseries('/protein/coords/PEPC10CA_y', Timeseries.Atom("y", trj.selectAtoms("atom PEPC 10 CA")))
        analysis.add_timeseries('/protein/coords/PEPC10CA_z', Timeseries.Atom("z", trj.selectAtoms("atom PEPC 10 CA")))
        
        analysis.add_timeseries('/protein/coords/PEPA549CA_x', Timeseries.Atom("x", trj.selectAtoms("atom PEPA 549 CA")))
        analysis.add_timeseries('/protein/coords/PEPA549CA_y', Timeseries.Atom("y", trj.selectAtoms("atom PEPA 549 CA")))
        analysis.add_timeseries('/protein/coords/PEPA549CA_z', Timeseries.Atom("z", trj.selectAtoms("atom PEPA 549 CA")))

        analysis.run(trj=trj, ref=ref)
        analysis.save()
        analysis.close()    
    

if __name__ == '__main__':
    main()
