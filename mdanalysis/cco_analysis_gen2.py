#!/usr/bin/python
# -*- coding: utf-8 -*-

# import required modules
import optparse # command-line parsing
import numpy # numerical operations
from tables import * # pytables (database)
from MDAnalysis import * # analysis of NAMD trajectories

# my analysis modules (in the same directory as this file)
from analysis import *
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
    
    # RMSD
    # broken: check new MDAnalysis RMSD code
    #analysis.add_to_sequence('/protein/rmsd/backbone', RMSD(selection='backbone'))
    
    # Nov. 8, 2010: saw that there is a potassium ion coordinated by N139D and D132
    # Count potassium ions near 132 and 139
    analysis.add_to_sequence('/protein/ions/PEPA_139_POT', NearbyCountAnalysis('segid PEPA and resid 139', 'resname POT', cutoff=4.0))
    analysis.add_to_sequence('/protein/ions/PEPA_132_POT', NearbyCountAnalysis('segid PEPA and resid 132', 'resname POT', cutoff=4.0))
    analysis.add_to_sequence('/protein/ions/PEPA_286_POT', NearbyCountAnalysis('segid PEPA and resid 286', 'resname POT', cutoff=4.0))
    analysis.add_to_sequence('/protein/ions/PEPA_207_POT', NearbyCountAnalysis('segid PEPA and resid 207', 'resname POT', cutoff=4.0))
        
    # HBONDS: 139, 207, 286 <-> water/protein
    analysis.add_to_sequence('/protein/hbonds/PEPA_132_WATER', HydrogenBondAnalysis('segid PEPA and resid 132', 'resname TIP3'))
    analysis.add_to_sequence('/protein/hbonds/PEPA_132_PROTEIN', HydrogenBondAnalysis('segid PEPA and resid 132', 'protein'))
    analysis.add_to_sequence('/protein/hbonds/PEPA_139_WATER', HydrogenBondAnalysis('segid PEPA and resid 139', 'resname TIP3'))
    analysis.add_to_sequence('/protein/hbonds/PEPA_139_PROTEIN', HydrogenBondAnalysis('segid PEPA and resid 139', 'protein'))
    analysis.add_to_sequence('/protein/hbonds/PEPA_207_WATER', HydrogenBondAnalysis('segid PEPA and resid 207', 'resname TIP3'))
    analysis.add_to_sequence('/protein/hbonds/PEPA_207_PROTEIN', HydrogenBondAnalysis('segid PEPA and resid 207', 'protein'))
    analysis.add_to_sequence('/protein/hbonds/PEPA_286_WATER', HydrogenBondAnalysis('segid PEPA and resid 286', 'resname TIP3'))
    analysis.add_to_sequence('/protein/hbonds/PEPA_286_PROTEIN', HydrogenBondAnalysis('segid PEPA and resid 286', 'protein'))
    
    # HEME A and A3 <-> water/protein
    analysis.add_to_sequence('/protein/hbonds/HEMEA3_PROA_WATER', HydrogenBondAnalysis('atom HEM3 1 O1A or atom HEM3 1 O2A', 'resname TIP3', selection1_type='acceptor'))
    analysis.add_to_sequence('/protein/hbonds/HEMEA3_PROA_PROTEIN', HydrogenBondAnalysis('atom HEM3 1 O1A or atom PEHEM3PA 1 O2A', 'protein', selection1_type='acceptor'))
    analysis.add_to_sequence('/protein/hbonds/HEMEA3_PROD_WATER', HydrogenBondAnalysis('atom HEM3 1 O1D or atom HEM3 1 O2D', 'resname TIP3', selection1_type='acceptor'))
    analysis.add_to_sequence('/protein/hbonds/HEMEA3_PROD_PROTEIN', HydrogenBondAnalysis('atom HEM3 1 O1D or atom HEM3 1 O2D', 'protein', selection1_type='acceptor'))
    analysis.add_to_sequence('/protein/hbonds/HEMEA_PROA_WATER', HydrogenBondAnalysis('atom HEMA 1 O1A or atom HEMA 1 O2A', 'resname TIP3', selection1_type='acceptor'))
    analysis.add_to_sequence('/protein/hbonds/HEMEA_PROA_PROTEIN', HydrogenBondAnalysis('atom HEMA 1 O1A or atom HEMA 1 O2A', 'protein', selection1_type='acceptor'))
    analysis.add_to_sequence('/protein/hbonds/HEMEA_PROD_WATER', HydrogenBondAnalysis('atom HEMA 1 O1D or atom HEMA 1 O2D', 'resname TIP3', selection1_type='acceptor'))
    analysis.add_to_sequence('/protein/hbonds/HEMEA_PROD_PROTEIN', HydrogenBondAnalysis('atom HEMA 1 O1D or atom HEMA 1 O2D', 'protein', selection1_type='acceptor'))
    
    # dihedrals of 132, 139, 286
    analysis.add_timeseries('/protein/dihedrals/PEPA_139_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 N", "atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG")))
    analysis.add_timeseries('/protein/dihedrals/PEPA_139_CHI2', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG", "atom PEPA 139 OD1")))
    analysis.add_timeseries('/protein/dihedrals/PEPA_132_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 132 N", "atom PEPA 132 CA", "atom PEPA 132 CB", "atom PEPA 132 CG")))
    analysis.add_timeseries('/protein/dihedrals/PEPA_132_CHI2', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 132 CA", "atom PEPA 132 CB", "atom PEPA 132 CG", "atom PEPA 132 OD1")))
    analysis.add_timeseries('/protein/dihedrals/PEPA_286_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 286 N", "atom PEPA 286 CA", "atom PEPA 286 CB", "atom PEPA 286 CG")))
    analysis.add_timeseries('/protein/dihedrals/PEPA_286_CHI2', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 286 CA", "atom PEPA 286 CB", "atom PEPA 286 CG", "atom PEPA 286 CD")))
    
    # dihedrals of proprionates from Heme A and Heme A3
    analysis.add_timeseries('/protein/dihedrals/HEMEA3_PROA_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom HEM3 1 C1A", "atom HEM3 1 C2A", "atom HEM3 1 CAA", "atom HEM3 1 CBA")))
    analysis.add_timeseries('/protein/dihedrals/HEMEA3_PROA_CHI2', Timeseries.Dihedral(trj.selectAtoms("atom HEM3 1 C2A", "atom HEM3 1 CAA", "atom HEM3 1 CBA", "atom HEM3 1 CGA")))
    analysis.add_timeseries('/protein/dihedrals/HEMEA3_PROD_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom HEM3 1 C4D", "atom HEM3 1 C3D", "atom HEM3 1 CAD", "atom HEM3 1 CBD")))
    analysis.add_timeseries('/protein/dihedrals/HEMEA3_PROD_CHI2', Timeseries.Dihedral(trj.selectAtoms("atom HEM3 1 C3D", "atom HEM3 1 CAD", "atom HEM3 1 CBD", "atom HEM3 1 CGD")))
    analysis.add_timeseries('/protein/dihedrals/HEMEA_PROA_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom HEMA 1 C1A", "atom HEMA 1 C2A", "atom HEMA 1 CAA", "atom HEMA 1 CBA")))
    analysis.add_timeseries('/protein/dihedrals/HEMEA_PROA_CHI2', Timeseries.Dihedral(trj.selectAtoms("atom HEMA 1 C2A", "atom HEMA 1 CAA", "atom HEMA 1 CBA", "atom HEMA 1 CGA")))
    analysis.add_timeseries('/protein/dihedrals/HEMEA_PROD_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom HEMA 1 C4D", "atom HEMA 1 C3D", "atom HEMA 1 CAD", "atom HEMA 1 CBD")))
    analysis.add_timeseries('/protein/dihedrals/HEMEA_PROD_CHI2', Timeseries.Dihedral(trj.selectAtoms("atom HEMA 1 C3D", "atom HEMA 1 CAD", "atom HEMA 1 CBD", "atom HEMA 1 CGD")))
    
    # distances between atoms that are distance restrained
    analysis.add_timeseries('/protein/distances/CUA1_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom CUA 1 CU", "atom CUA 2 CU")))
    
    analysis.add_timeseries('/protein/distances/HSD102NE2_HEMAFE', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 102 NE2", "atom HEMA 1 FE")))
    analysis.add_timeseries('/protein/distances/HSD421NE2_HEMAFE', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 421 NE2", "atom HEMA 1 FE")))
    analysis.add_timeseries('/protein/distances/HSD419NE2_HEM3FE', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 419 NE2", "atom HEM3 1 FE")))
    analysis.add_timeseries('/protein/distances/CUB_HEM3FE', Timeseries.Distance("r", trj.selectAtoms("atom CUB 1 CU", "atom HEM3 1 FE")))
    analysis.add_timeseries('/protein/distances/HSD411NE2_MG', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 411 NE2", "atom MG 1 MG")))
    analysis.add_timeseries('/protein/distances/ASP412OD2_MG', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 412 OD2", "atom MG 1 MG")))
    analysis.add_timeseries('/protein/distances/HSD333NE2_CUB', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 333 NE2", "atom CUB 1 CU")))
    analysis.add_timeseries('/protein/distances/HSE284ND1_CUB', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 284 ND1", "atom CUB 1 CU")))
    analysis.add_timeseries('/protein/distances/HSD334NE2_CUB', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 334 NE2", "atom CUB 1 CU")))
    analysis.add_timeseries('/protein/distances/CYS252SG_CUA1', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 252 SG", "atom CUA 1 CU")))
    analysis.add_timeseries('/protein/distances/CYS252SG_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 252 SG", "atom CUA 2 CU")))
    analysis.add_timeseries('/protein/distances/CYS256SG_CUA1', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 256 SG", "atom CUA 1 CU")))
    analysis.add_timeseries('/protein/distances/CYS256SG_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 256 SG", "atom CUA 2 CU")))
    analysis.add_timeseries('/protein/distances/HSE217ND1_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 217 ND1", "atom CUA 2 CU")))
    analysis.add_timeseries('/protein/distances/HSE260ND1_CUA1', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 260 ND1", "atom CUA 1 CU")))
    analysis.add_timeseries('/protein/distances/MET263SD_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 263 SD", "atom CUA 2 CU")))
    analysis.add_timeseries('/protein/distances/GLU254O_CUA1', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 254 O", "atom CUA 1 CU")))
    analysis.add_timeseries('/protein/distances/E254OE1_MG', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 254 OE1", "atom MG 1 MG")))
    
    # Nov 18, 2010, cylinder search for water and potassium
    r132 = 'segid PEPA and resid 132 and ( name CA or name CB or name N )'
    #r139 = 'segid PEPA and resid 139 and ( name CA or name CB or name N )'
    r139 = 'segid PEPA and resid 139 and ( name CA or name N )'
    r286 = 'segid PEPA and resid 286 and ( name CA or name CB or name N )'
    r312 = 'segid PEPA and resid 312 and ( name CA or name CB or name N )'
    r481 = 'segid PEPA and resid 481 and ( name CA or name CB or name N )'
    rHA3FE = 'atom HEM3 1 FE'

    # D-channel
    analysis.add_to_sequence('/arrays/cylinder_132_286', CylinderSearch(r132, r286, 'resname TIP3 or resname POT', extension=5.0, radius=10.0), array=True)
    analysis.add_to_sequence('/arrays/cylinder_132_139', CylinderSearch(r132, r139, 'resname TIP3 or resname POT', extension=1.0, radius=5.0), array=True)
    analysis.add_to_sequence('/arrays/cylinder_139_286', CylinderSearch(r139, r286, 'resname TIP3 or resname POT', extension=1.0, radius=5.0), array=True)

    # GLY312 (entrance to K-channel) to FE of Heme A3 (K-channel)
    analysis.add_to_sequence('/arrays/cylinder_312_HA3FE', CylinderSearch(r312, rHA3FE, 'resname TIP3 or resname POT', extension=1.0, radius=5.0), array=True)

    # Active site (286 to 481)
    analysis.add_to_sequence('/arrays/cylinder_286_481', CylinderSearch(r286, r481, 'resname TIP3 or resname POT', extension=1.0, radius=5.0), array=True)

    analysis.run(trj=trj, ref=ref)
    analysis.save()
    analysis.close()    
    

if __name__ == '__main__':
    main()
