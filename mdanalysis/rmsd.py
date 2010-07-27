# RMSD analysis class
from analysis import Analysis

import numpy
import numpy.linalg
from MDAnalysis import *
from MDAnalysis.core.AtomGroup import Residue, AtomGroup
import MDAnalysis.core.rms_fitting
    
class FrameData(object):
    atoms = None
    coordinates = None
    masses = None
    com = None
    rmsd = None
    
    def __init__(self, atom_group, allocate_only=False):
        self.atoms = atom_group
        if allocate_only:
            self.coordinates = atom_group.coordinates().copy()
        else:
            self.masses = atom_group.masses()
            self.com = atom_group.centerOfMass().astype(numpy.float32)
            self.coordinates = atom_group.coordinates() - self.com

class RMSD(Analysis):
    table_name = 'RMSD'
    _selection = None
    _name = None
    _rmsds = []
    
    def _rmsd(self, a,b):
        """Returns RMSD between two coordinate sets a and b."""
        return numpy.sqrt(numpy.sum(numpy.power(a-b,2))/a.shape[0])
    
    def __init__(self, selection, name):
        self._selection = selection
        self._name = name
    
    def prepare(self, ref, trj):
        ref_atoms = ref.selectAtoms(self._selection)
        trj_atoms = trj.selectAtoms(self._selection)
        self.fit_ref = FrameData(ref.selectAtoms('backbone'))
        self.fit_trj = FrameData(trj.selectAtoms('backbone'))
        self.rmsd_ref = FrameData(ref_atoms)
        self.rmsd_trj = FrameData(trj_atoms, allocate_only=True)
        # print "Done RMSD prepare."
        
    def process(self, ts):
        # print "RMSD Fitting Frame %5d" % (ts.frame) 
        x_com = self.fit_trj.atoms.centerOfMass().astype(numpy.float32)
        self.fit_trj.coordinates[:] = self.fit_trj.atoms.coordinates() - x_com
        R = numpy.matrix(MDAnalysis.core.rms_fitting.rms_rotation_matrix(self.fit_trj.coordinates, self.fit_ref.coordinates, self.fit_ref.masses),dtype=numpy.float32)
        ts._pos   -= x_com
        ts._pos[:] = ts._pos * R
        ts._pos   += self.fit_ref.com
        self._rmsds.append(self._rmsd(self.rmsd_ref.atoms.coordinates(),self.rmsd_trj.atoms.coordinates()))
        # print self._rmsds[-1]

    def results(self):
        return numpy.array(self._rmsds)