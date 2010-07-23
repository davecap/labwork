# RMSD analysis class
from analysis import Analysis

import numpy
import numpy.linalg
from MDAnalysis import *
from MDAnalysis.core.AtomGroup import Residue, AtomGroup
import MDAnalysis.core.rms_fitting

def get_residues_for_atoms(atoms):
    # Get all protein segments and residues
    residues = []
    last_res = -1
    last_seg = "NOSEG"
    for a in atoms.atoms:
        if last_res == a.resid and a.segment.name == last_seg:
            continue
        residues.append(a.residue)
        last_res = a.resid
        last_seg = a.segment.name
    return residues

class RMSD(Analysis):
    ref_atoms = []
    traj_atoms = []
    masses = []
    ref_com = []
    ref_coordinates = []
    traj_coordinates = []
    rmsds = []
    select = None
    
    def rmsd(self, a,b):
        """Returns RMSD between two coordinate sets a and b."""
        return numpy.sqrt(numpy.sum(numpy.power(a-b,2))/a.shape[0])
    
    def __init__(self):
        print "__init__ RMSD called"
    
    def prepare(self, ref, trj):
        print "Preparing RMSD ref and trj"
        print "Initializing RMSD selects..."
        print "Getting all residues in ref..."
        ref_residues = get_residues_for_atoms(ref.selectAtoms('protein'))
        print "Found %d protein residues in ref." % len(ref_residues)

        print "Getting all residues in trj..."
        trj_residues = get_residues_for_atoms(trj.selectAtoms('protein'))
        print "Found %d protein residues in ref." % len(trj_residues)
        self.select = zip(ref_residues, trj_residues)
        
        if len(self.select) == 0:
            raise SelectionError("No atoms selected!")
        if self.select[0] is not 'backbone':
            self.select.insert(0,'backbone')
        
        for s in self.select:
            self.rmsds.append([])
            # if it is a string then we have to search for it
            if type(s) is str:
                self.ref_atoms.append(ref.selectAtoms(s))
                self.traj_atoms.append(trj.selectAtoms(s))
            elif type(s) is tuple:
                #otherwise it is a Residue/AtomGroup object
                self.ref_atoms.append(s[0])
                self.traj_atoms.append(s[1])
            else:
                print "Invalid selection for RMSD: %s" % (type(s))
                continue

            if len(self.ref_atoms[-1]) != len(self.traj_atoms[-1]):
                raise SelectionError("Reference and trajectory atom selections do not contain "+
                                     "the same number of atoms: N_ref=%d, N_traj=%d" % \
                                     (len(self.ref_atoms[-1]), len(self.traj_atoms[-1])))

            # TODO: could check that ref and traj selection have same masses (see rmsfit_alignment.py)
            self.masses.append(self.ref_atoms[-1].masses())

            # reference centre of mass system
            # (compatibility with pre 1.0 numpy: explicitly cast coords to float32)
            self.ref_com.append(self.ref_atoms[-1].centerOfMass().astype(numpy.float32))
            self.ref_coordinates.append(self.ref_atoms[-1].coordinates() - self.ref_com[-1])

            # allocate the array for selection atom coords
            # (allocating and re-filling is MUCH more efficient that re-allocating for 
            # every frame!)
            self.traj_coordinates.append(self.traj_atoms[-1].coordinates().copy())
        print "Done RMSD prepare."
        

    def process(self, ts):
        print "Processing frame RMSD"
        print "RMSD Fitting Frame %5d" % (ts.frame)
        
        # only align once per frame (backbone)
        # index 0 is the backbone selection
        # shift coordinates for rotation fitting
        # selection is updated with the time frame
        x_com = self.traj_atoms[0].centerOfMass().astype(numpy.float32)
        # re-fill the coordinates array
        self.traj_coordinates[0][:] = self.traj_atoms[0].coordinates() - x_com
        R = numpy.matrix(MDAnalysis.core.rms_fitting.rms_rotation_matrix(self.traj_coordinates[0],self.ref_coordinates[0],self.masses[0]),dtype=numpy.float32)
        # Transform each atom in the trajectory (use inplace ops to avoid copying arrays)
        # Note that this destroys the original coordinate information in the current timestep.
        ts._pos   -= x_com
        ts._pos[:] = ts._pos * R # R acts to the left (!) & is broadcasted N times (numpy magic!)
        ts._pos   += self.ref_com[0]
        
        for i in range(len(self.ref_atoms)):            
            rmsd_new = self.rmsd(self.ref_atoms[i].coordinates(),self.traj_atoms[i].coordinates())
            self.rmsds[i].append(rmsd_new)

    def results(self):
        print "Returning RMSD results"
        return numpy.array(self.rmsds)