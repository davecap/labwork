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
    parser.add_option("-o", dest="h5_filename", default="analysis.h5", help="Output analysis file [default: %default]")    
    
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
    
    # TODO: check to see if we have already analyzed the given dcd file
    #   run md5 in subprocess...
    
    analysis = Analysis(options.h5_filename, readonly=False)
    
    analysis.add_metadata('/metadata/trajectory', { 'psf': psf_file, 'pdb': pdb_file, 'dcd': dcd_file, 'frames': num_frames, 'firsttimestep': first_timestep, 'dt': dt })
    
    # RMSD
    # this takes too long so I'm disabling it for now
    #analysis.add_to_sequence('/protein/rmsd/backbone', RMSD(selection='backbone'))
    
    # dihedrals of 132, 139, 286
    analysis.add_timeseries('/protein/dihedrals/PEPA_139_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 N", "atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG")))
    analysis.add_timeseries('/protein/dihedrals/PEPA_139_CHI2', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 139 CA", "atom PEPA 139 CB", "atom PEPA 139 CG", "atom PEPA 139 OD1")))
    analysis.add_timeseries('/protein/dihedrals/PEPA_132_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 132 N", "atom PEPA 132 CA", "atom PEPA 132 CB", "atom PEPA 132 CG")))
    analysis.add_timeseries('/protein/dihedrals/PEPA_132_CHI2', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 132 CA", "atom PEPA 132 CB", "atom PEPA 132 CG", "atom PEPA 132 OD1")))
    analysis.add_timeseries('/protein/dihedrals/PEPA_286_CHI1', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 286 N", "atom PEPA 286 CA", "atom PEPA 286 CB", "atom PEPA 286 CG")))
    analysis.add_timeseries('/protein/dihedrals/PEPA_286_CHI2', Timeseries.Dihedral(trj.selectAtoms("atom PEPA 286 N", "atom PEPA 286 CA", "atom PEPA 286 CB", "atom PEPA 286 CG")))
    
    # distances between atoms that are distance restrained
    analysis.add_timeseries('/protein/distances/CUA1_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom J 1 CU", "atom J 2 CU")))
    analysis.add_timeseries('/protein/distances/HDH102NE2_HDH102FE', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 102 NE2", "atom PEPA 102 FE")))
    analysis.add_timeseries('/protein/distances/HDH419NE2_HDH419FE', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 419 NE2", "atom PEPA 419 FE")))
    analysis.add_timeseries('/protein/distances/CUB_HDH419FE', Timeseries.Distance("r", trj.selectAtoms("atom I 3 CU", "atom PEPA 419 FE")))
    analysis.add_timeseries('/protein/distances/HSD411NE2_MG', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 411 NE2", "atom I 1 MG")))
    analysis.add_timeseries('/protein/distances/ASP412OD2_MG', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 412 OD2", "atom I 1 MG")))
    analysis.add_timeseries('/protein/distances/HSD333NE2_CUB', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 333 NE2", "atom I 3 CU")))
    analysis.add_timeseries('/protein/distances/HSE284ND1_CUB', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 284 ND1", "atom I 3 CU")))
    analysis.add_timeseries('/protein/distances/HSD334NE2_CUB', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 334 NE2", "atom I 3 CU")))
    analysis.add_timeseries('/protein/distances/CYS252SG_CUA1', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 252 SG", "atom J 1 CU")))
    analysis.add_timeseries('/protein/distances/CYS252SG_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 252 SG", "atom J 2 CU")))
    analysis.add_timeseries('/protein/distances/CYS256SG_CUA1', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 256 SG", "atom J 1 CU")))
    analysis.add_timeseries('/protein/distances/CYS256SG_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 256 SG", "atom J 2 CU")))
    analysis.add_timeseries('/protein/distances/HSE217ND1_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 217 ND1", "atom J 2 CU")))
    analysis.add_timeseries('/protein/distances/HSE260ND1_CUA1', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 260 ND1", "atom J 1 CU")))
    analysis.add_timeseries('/protein/distances/MET263SD_CUA2', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 263 SD", "atom J 2 CU")))
    analysis.add_timeseries('/protein/distances/GLU254O_CUA1', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 254 O", "atom J 1 CU")))
    analysis.add_timeseries('/protein/distances/E254OE1_MG', Timeseries.Distance("r", trj.selectAtoms("atom PEPB 254 OE1", "atom I 1 MG")))
    
    # distances between proprionates and some residues
    analysis.add_timeseries('/protein/distances/HDH419O1A_HSD411HD1', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 419 O1A", "atom PEPA 411 HD1")))
    analysis.add_timeseries('/protein/distances/HDH419O2A_HSD411HD1', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 419 O2A", "atom PEPA 411 HD1")))
    analysis.add_timeseries('/protein/distances/HDH419O1A_ASPP407HD2', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 419 O1A", "atom PEPA 407 HD2")))
    analysis.add_timeseries('/protein/distances/HDH419O1A_ASPP407OD1', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 419 O1A", "atom PEPA 407 OD1")))
    analysis.add_timeseries('/protein/distances/HDH419O2A_ASPP407HD2', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 419 O2A", "atom PEPA 407 HD2")))
    analysis.add_timeseries('/protein/distances/HDH419O2A_ASPP407OD1', Timeseries.Distance("r", trj.selectAtoms("atom PEPA 419 O2A", "atom PEPA 407 OD1")))
    
    analysis.run(trj=trj, ref=ref)
    analysis.save()
    analysis.close()    
    

if __name__ == '__main__':
    main()