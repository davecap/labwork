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
    table_name = 'RMSD'
    description = 'RMSDs of backbone and all protein residues'
    _data = {}
    
    def _rmsd(self, a,b):
        """Returns RMSD between two coordinate sets a and b."""
        return numpy.sqrt(numpy.sum(numpy.power(a-b,2))/a.shape[0])
    
    def __init__(self):
        pass
    
    def prepare(self, ref, trj):
        print "Preparing RMSD ref and trj"
        ref_residues = get_residues_for_atoms(ref.selectAtoms('protein'))
        print "Found %d protein residues in ref." % len(ref_residues)

        print "Getting all residues in trj..."
        trj_residues = get_residues_for_atoms(trj.selectAtoms('protein'))
        print "Found %d protein residues in ref." % len(trj_residues)
        
        if len(ref_residues) != len(trj_residues):
            raise SelectionError("Number of residues is difference between reference (%d) and trajectory (%d)." % (len(ref_residues), len(trj_residues)))
        
        ref_trj_residues = zip(ref_residues, trj_residues)
        ref_trj_residues.insert(0, (ref.selectAtoms('backbone'), trj.selectAtoms('backbone')))
        
        for i, r in enumerate(ref_trj_residues):
            ref_atoms = r[0]
            trj_atoms = r[1]
            
            if i == 0:
                k = 'backbone'
            else:
                k = '%s_%d_%s' % (ref_atoms.segment.name, ref_atoms.id, ref_atoms.name)
            ref_com = ref_atoms.centerOfMass().astype(numpy.float32)
            self._data[k] = {   'rmsds': [],
                                'ref': {    'atoms': ref_atoms, 
                                            'masses': ref_atoms.masses(), 
                                            'com': ref_com, 
                                            'coordinates': ref_atoms.coordinates() - ref_com,
                                        },
                                'trj': {    'atoms': trj_atoms,
                                            'coordinates': trj_atoms.coordinates().copy(),
                                        },
                            }
        print "Done RMSD prepare."
        

    def process(self, ts):
        print "RMSD Fitting Frame %5d" % (ts.frame) 
        
        x_com = self._data['backbone']['trj']['atoms'].centerOfMass().astype(numpy.float32)
        self._data['backbone']['trj']['coordinates'][:] = self._data['backbone']['trj']['atoms'].coordinates() - x_com
        R = numpy.matrix(MDAnalysis.core.rms_fitting.rms_rotation_matrix(self._data['backbone']['trj']['coordinates'], self._data['backbone']['ref']['coordinates'], self._data['backbone']['ref']['masses']),dtype=numpy.float32)
        ts._pos   -= x_com
        ts._pos[:] = ts._pos * R
        ts._pos   += self._data['backbone']['ref']['com']
        
        for k, v in self._data.items():
            rmsd = self._rmsd(v['ref']['atoms'].coordinates(),v['trj']['atoms'].coordinates())
            v['rmsds'].append(rmsd)

    def results(self):
        return [ (k, numpy.array(v['rmsds'])) for k,v in self._data.items() ]