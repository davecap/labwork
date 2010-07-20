#!/usr/bin/env python

from MDAnalysis import *
from MDAnalysis import collection, SelectionError
import MDAnalysis.core.rms_fitting

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

def rmsd_trj(traj,ref,select='backbone',):
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

    # setup new dcd
    frames = traj.trajectory

    ref_atoms = ref.selectAtoms(select)
    traj_atoms = traj.selectAtoms(select)
    if len(ref_atoms) != len(traj_atoms):
        raise SelectionError("Reference and trajectory atom selections do not contain "+
                             "the same number of atoms: N_ref=%d, N_traj=%d" % \
                             (len(ref_atoms), len(traj_atoms)))

    # could check that ref and traj selection have same masses (see rmsfit_alignment.py)
    masses = ref_atoms.masses()

    # reference centre of mass system
    # (compatibility with pre 1.0 numpy: explicitly cast coords to float32)
    ref_com = ref_atoms.centerOfMass().astype(numpy.float32)
    ref_coordinates = ref_atoms.coordinates() - ref_com

    # allocate the array for selection atom coords
    # (allocating and re-filling is MUCH more efficient that re-allocating for 
    # every frame!)
    traj_coordinates = traj_atoms.coordinates().copy()

    # collect the time series
    rmsds = []

    # R: rotation matrix that aligns r-r_com, x~-x~com   
    #    (x~: selected coordinates, x: all coordinates)
    # Final transformed traj coordinates: x' = (x-x~_com)*R + ref_com
    for ts in frames:
        # shift coordinates for rotation fitting
        # selection is updated with the time frame
        x_com = traj_atoms.centerOfMass().astype(numpy.float32)
        # re-fill the coordinates array
        traj_coordinates[:] = traj_atoms.coordinates() - x_com
        R = numpy.matrix(MDAnalysis.core.rms_fitting.rms_rotation_matrix(
                traj_coordinates,ref_coordinates,masses),dtype=numpy.float32)
        # Transform each atom in the trajectory (use inplace ops to avoid copying arrays)
        # Note that this destroys the original coordinate information in the current timestep.
        ts._pos   -= x_com
        ts._pos[:] = ts._pos * R # R acts to the left (!) & is broadcasted N times (numpy magic!)
        ts._pos   += ref_com

        rmsd_old = rmsd(ref_atoms.coordinates(),traj_coordinates)
        rmsd_new = rmsd(ref_atoms.coordinates(),traj_atoms.coordinates())

        rmsds.append(rmsd_new)

        #print "Fitted frame %5d/%d  [%5.1f%%]  %5.2fA --> %5.2fA  |translation|=%.2fA\r" % (ts.frame,frames.numframes,100.0*ts.frame/frames.numframes, rmsd_old, rmsd_new, rmsd(x_com,ref_com))
        #if ts.frame % 10 == 0:
        #    echo("Fitted frame %5d/%d  [%5.1f%%]\r" % \
        #        (ts.frame,frames.numframes,100.0*ts.frame/frames.numframes))
    print "\n"
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
