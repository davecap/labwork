#!/usr/bin/env python

from MDAnalysis import *
from MDAnalysis import collection, SelectionError
import MDAnalysis.core.rms_fitting
from MDAnalysis.core.AtomGroup import Residue, AtomGroup

import numpy.linalg
import scipy.stats
import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
import optparse
import numpy

import sys

def rmsd(a,b):
    """Returns RMSD between two coordinate sets a and b."""
    return numpy.sqrt(numpy.sum(numpy.power(a-b,2))/a.shape[0])

def rmsd_trj(traj,ref,select=['backbone'],):
    """Calculate the RMSD of a selection over a trajectory.

      rmsd_trj(traj, ref, 'backbone or name CB or name OT*')

    :Arguments:
      traj
         trajectory, :class:`MDAnalysis.Universe` object
      ref
         reference coordinates; :class:`MDAnalysis.Universe` object
         (uses the current time step of the object)
      select    
         any valid selection string for
         MDAnalysis.AtomGroup.selectAtom that produces identical
         selections in traj and ref or dictionary {'reference':sel1,
         'target':sel2}.  The fasta2select() function returns such a
         dictionary based on a ClustalW or STAMP sequence alignment.

    Both reference and trajectory must be MDAnalysis.Universe
    instances.
    """
    if len(select) == 0:
        raise SelectionError("No atoms selected!")
    if select[0] is not 'backbone':
        select.insert(0,'backbone')
    
    # setup new dcd
    frames = traj.trajectory
    ref_atoms = []
    traj_atoms = []
    masses = []
    ref_com = []
    ref_coordinates = []
    traj_coordinates = []
    rmsds = []
    
    print "Initializing RMSD selects..."
    for s in select:
        rmsds.append([])
        # if it is a string then we have to search for it
        if type(s) is str:
            ref_atoms.append(ref.selectAtoms(s))
            traj_atoms.append(traj.selectAtoms(s))
        elif type(s) is tuple:
            #otherwise it is a Residue/AtomGroup object
            ref_atoms.append(s[0])
            traj_atoms.append(s[1])
        else:
            print "Invalid selection for RMSD: %s" % (type(s))
            continue
            
        if len(ref_atoms[-1]) != len(traj_atoms[-1]):
            raise SelectionError("Reference and trajectory atom selections do not contain "+
                                 "the same number of atoms: N_ref=%d, N_traj=%d" % \
                                 (len(ref_atoms[-1]), len(traj_atoms[-1])))

        # TODO: could check that ref and traj selection have same masses (see rmsfit_alignment.py)
        masses.append(ref_atoms[-1].masses())

        # reference centre of mass system
        # (compatibility with pre 1.0 numpy: explicitly cast coords to float32)
        ref_com.append(ref_atoms[-1].centerOfMass().astype(numpy.float32))
        ref_coordinates.append(ref_atoms[-1].coordinates() - ref_com[-1])
    
        # allocate the array for selection atom coords
        # (allocating and re-filling is MUCH more efficient that re-allocating for 
        # every frame!)
        traj_coordinates.append(traj_atoms[-1].coordinates().copy())
    print "Done."
    
    # R: rotation matrix that aligns r-r_com, x~-x~com   
    #    (x~: selected coordinates, x: all coordinates)
    # Final transformed traj coordinates: x' = (x-x~_com)*R + ref_com
    for ts in frames:
        print "RMSD Fitting Frame %5d/%d [%5.1f%%]" % (ts.frame,frames.numframes, 100.0*ts.frame/frames.numframes)
        
        # only align once per frame (backbone)
        # index 0 is the backbone selection
        # shift coordinates for rotation fitting
        # selection is updated with the time frame
        x_com = traj_atoms[0].centerOfMass().astype(numpy.float32)
        # re-fill the coordinates array
        traj_coordinates[0][:] = traj_atoms[0].coordinates() - x_com
        R = numpy.matrix(MDAnalysis.core.rms_fitting.rms_rotation_matrix(traj_coordinates[0],ref_coordinates[0],masses[0]),dtype=numpy.float32)
        # Transform each atom in the trajectory (use inplace ops to avoid copying arrays)
        # Note that this destroys the original coordinate information in the current timestep.
        ts._pos   -= x_com
        ts._pos[:] = ts._pos * R # R acts to the left (!) & is broadcasted N times (numpy magic!)
        ts._pos   += ref_com[0]
        
        for i in range(len(ref_atoms)):            
            # rmsd_old = rmsd(ref_atoms[i].coordinates(),traj_coordinates[i])
            rmsd_new = rmsd(ref_atoms[i].coordinates(),traj_atoms[i].coordinates())
            #print "RMSD for selection: %s => %5.2fA" % (select[i], rmsd_new)
            rmsds[i].append(rmsd_new)
            #print "Fitted frame %5d/%d  [%5.1f%%]  %5.2fA --> %5.2fA  |translation|=%.2fA\r" % (ts.frame,frames.numframes,100.0*ts.frame/frames.numframes, rmsd_old, rmsd_new, rmsd(x_com,ref_com))
            #if ts.frame % 10 == 0:
            #    echo("Fitted frame %5d/%d  [%5.1f%%]\r" % \
            #        (ts.frame,frames.numframes,100.0*ts.frame/frames.numframes))
    return numpy.array(rmsds)

if __name__ == '__main__':
    usage = """
        usage: %prog [options] <PSF file> <DCD file>
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-p", dest="pdb_file", default=None, help="Reference PDB file [default: use first frame of trj]")    
    parser.add_option("-s", dest="selection", default="name CA", help="Atom selection [default: %default]")    
    
    options, args = parser.parse_args()
    
    if len(args) < 2:
        parser.error("No input files specified")
    
    psf_file = args[0]
    dcd_file = args[1]

    if options.pdb_file:
        print "Using %s as reference frame." % pdb_file
        ref = Universe(psf_file, pdb_file)
    else:
        print "Using first frame as reference"
        ref = Universe(psf_file, dcd_file)
        ref.trajectory[0]
    
    print "Initializing trajectory: %s" % dcd_file
    trj = Universe(psf_file, dcd_file)
    
    rmsds = rmsd_trj(trj, ref, select=options.selection)
    
    plt.figure(1)
    plt.plot(rmsds)
    plt.title('RMSD for trajectory: %s' % (dcd_file))
    plt.xlabel('frame')
    plt.ylabel(r'RMSD $\AA$')
    plt.show()
