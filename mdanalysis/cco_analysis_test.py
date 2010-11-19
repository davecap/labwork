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
        usage: %prog [options] <PSF file> <PDB file> <DCD file>
    """
    
    # initialize the parser with our custom usage string (above)
    parser = optparse.OptionParser(usage)
    # add all the command-line options
    parser.add_option("-o", dest="h5_filename", default="analysis.h5", help="Output analysis file [default: %default]")    
    parser.add_option("-s", dest="selection", default="name CA", help="Atom selection [default: %default]")    
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
    dcd_file = args[2]
    
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
    
    if not os.path.exists(dcd_file):
        print "Error: DCD file doesn't exist: %s" % dcd_file
        return
    elif not dcd_file.lower().endswith('.dcd'):
        print "Warning: DCD filename does not end in .dcd: %s" % dcd_file    

    # # check to make sure we haven't analyzed this DCD yet
    # if os.path.exists(options.h5_filename):
    #     # open the HDF5 database
    #     h5f = tables.openFile(options.h5_filename, mode="r")
    #     # get the metadata table
    #     tbl = h5f.getNode('/metadata/trajectory')
    #     # get the list of analyzed dcd files
    #     analyzed_dcd_files = tbl.read(field='dcd')
    #     h5f.close()
    #     if dcd_file.split('/')[-1] in analyzed_dcd_files:
    #        print "DCD file %s already analyzed... exiting!" % dcd_file
    #        return -1
    
    # load the PSF and PDB into the analysis system
    print "Loading reference system: %s, %s" % (psf_file, pdb_file)
    ref = Universe(psf_file, pdb_file, permissive=True)
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
    
    # Nov 18, 2010, cylinder search for water and potassium
    r132 = 'segid PEPA and resid 132 and ( name CA or name CB or name N )'
    r139 = 'segid PEPA and resid 139 and ( name CA or name CB or name N )'
    r286 = 'segid PEPA and resid 286 and ( name CA or name CB or name N )'
    r312 = 'segid PEPA and resid 312 and ( name CA or name CB or name N )'
    r481 = 'segid PEPA and resid 481 and ( name CA or name CB or name N )'
    rHA3FE = 'atom PEPA 419 FE'
    
    # D-channel
    analysis.add_to_sequence('/protein/cylinder_132_286/TIP3_POT', CylinderSearch(r132, r286, 'resname TIP3 or resname POT', extension=5.0, radius=10.0), array=True)
    analysis.add_to_sequence('/protein/cylinder_132_139/TIP3_POT', CylinderSearch(r132, r139, 'resname TIP3 or resname POT', extension=1.0, radius=5.0), array=True)
    analysis.add_to_sequence('/protein/cylinder_139_286/TIP3_POT', CylinderSearch(r139, r286, 'resname TIP3 or resname POT', extension=1.0, radius=5.0), array=True)
    
    # GLY312 (entrance to K-channel) to FE of Heme A3 (K-channel)
    analysis.add_to_sequence('/protein/cylinder_312_HA3FE/TIP3_POT', CylinderSearch(r312, rHA3FE, 'resname TIP3 or resname POT', extension=1.0, radius=5.0), array=True)
    
    # Active site (286 to 481)
    analysis.add_to_sequence('/protein/cylinder_286_481/TIP3', CylinderSearch(r286, r481, 'resname TIP3', extension=1.0, radius=5.0), array=True)
    
    
    analysis.run(trj=trj, ref=ref)
    analysis.save()
    analysis.close()    
    

if __name__ == '__main__':
    main()
