# Hydrogen Bonding Analysis
"""
Hydrogen Bond analysis --- :mod:`MDAnalysis.analysis.hbonds`
===================================================================

:Author: David Caplan
:Year: 2010
:Copyright: GNU Public License v3


This is modeled after the VMD HBONDS plugin (http://www.ks.uiuc.edu/Research/vmd/plugins/hbonds/)

Given a Universe (simulation trajectory with 1 or more frames) measure all hydrogen bonds for each frame between selections 1 and 2.

Options:
  - update_selections (True): update selections at each frame?
  - selection_1_type ('both'): selection 1 is the: donor, acceptor, both
  - donor-acceptor distance (A): 3.0
  - Angle cutoff (degrees): 120.0

Returns hydrogen bond data per frame:
    results[ [ <donor index>, <acceptor index>, <donor string>, <acceptor string>, <distance>, <angle> ], [frame 1], [frame 2] ... ]


Example
-------

TODO

  import MDAnalysis
  import hbonds
  
  u = MDAnalysis.Universe(PSF, PDB, permissive=True)
  h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis('protein', 'resname TIP3', distance=3.0, angle=120.0)
  results = h.run(u)

Classes
-------

.. autoclass:: HydrogenBondAnalysis
   :members:

"""

from __future__ import with_statement

import os
import warnings
import bz2
import numpy

import MDAnalysis
from MDAnalysis.core.distances import self_distance_array
from MDAnalysis.core.AtomGroup import AtomGroup
import MDAnalysis.KDTree.NeighborSearch as NS

import logging
logger = logging.getLogger('hbonds')

class HydrogenBondAnalysis(object):
    """Perform a hydrogen bond analysis

    The analysis of the trajectory is performed with the
    :meth:`HydrogenBondAnalysis.run` method. The result is stored in
    :attr:`HydrogenBondAnalysis.timeseries`. It is an array of one element per frame.
    Each frame (element) is an array of D-H-A tuples (donor, hydrogen, acceptor) for
    the hydrogen bonds found.
    
    The criteria for hydrogen bonds can be set by the DA_distance and HA_distance and angle.
    The DA_distance is the cutoff between donors and acceptors. 
    The HA_distance is the cutoff between hydrogens and acceptors.
    The angle is the D-H-A angle cutoff.
    All three must be satisfied to count as an h-bond.

    Donors: NH of the main chain, water-H1/H2, ARG NE, ASN ND2, HIS NE2, SER OG, TYR OH, ARG NH1, CYS SG, HIS ND1, THR OG1, ARG NH2, GLN NE2, LYS NZ, TRP NE1
    Acceptors: CO main chain, water-OH2, GLN OE1, MET SD, ASP OD1, GLU OE1, SER OG, ASP OD2, GLU OE2, THR OG1,CYH SG, HIS ND1, TYR OH.

    """
    
    donors = ('OD1', 'OD2', 'OW', 'OE1', 'OE2', 'NH', 'OH2', 'NE', 'ND2', 'NE2', 'OG', 'OH', 'NH1', 'SG', 'ND1', 'OG1', 'NH2', 'NE2', 'NZ', 'NE1', )
    acceptors = ('O1A', 'O2A', 'O1D', 'O2D', 'OW', 'OE1', 'OE2', 'CO', 'OH2', 'OD1', 'OE1', 'SD', 'OE1', 'OG', 'OD2', 'OE2', 'OG1', 'SG', 'ND1', 'OH', )
    
    def __init__(self, selection1='protein', selection2='all', selection1_type='both', update_selections=False, DA_distance=3.3, HA_distance=2.5, angle=140.0):
        """Calculate hydrogen bonds between two selections.

        :Arguments:
          *selection1*
            Selection string for first selection
          *selection2*
            Selection string for second selection
          *selection1_type*
            Selection 1 can be 'donor', 'acceptor' or 'both'
          *update_selections*
            Update selections at each frame?
          *DA_distance*
            Distance cutoff for Donor <-> Acceptor
          *HA_distance*
            Distance cutoff for Hydrogen <-> Acceptor
          *angle*
            Angle cutoff for hydrogen bonds
            
        The timeseries accessible as the attribute :attr:`HydrogenBondAnalysis.timeseries`.
        """

        self.selection1 = selection1
        self.selection2 = selection2
        self.selection1_type = selection1_type
        self.update_selections = update_selections
        self.DA_distance = DA_distance
        self.HA_distance = HA_distance
        self.angle = angle
        
        if not (self.selection1 and self.selection2):
            raise Exception('HydrogenBondAnalysis: invalid selections')
        elif self.selection1_type not in ('both', 'donor', 'acceptor'):
            raise Exception('HydrogenBondAnalysis: Invalid selection type %s' % self.selection1_type)

    def run(self, trj):
        """ Analyze trajectory and produce timeseries. """
        self.timeseries = []
        self.prepare(trj=trj)
        for ts in self.u.trajectory:
            logger.debug("Analyzing frame %d" % ts.frame)
            self.process(ts.frame)
        return self.timeseries

    def prepare(self, ref=None, trj=None):
        """ Prepare the trajectory (trj is a Universe object). No reference object is needed. """
        self.u = trj
        self.u.trajectory.rewind()
        self._update_selections()
        self.timeseries = []  # final result
        self._hydrogens = {}

    def process(self, frame):
        """ Process a single trajectory frame """
        frame_results = []
        # calculate donor/acceptor pairs
        pairs = {} # NOTE: atom pairs are stored as tuple keys to avoid duplicates
        for d in self._s1_donor_atoms:
            # search for S2 acceptors within <distance> of donor d
            self._s2_acceptor_atoms_ns = NS.AtomNeighborSearch(self._s2_acceptor_atoms)
            for a in self._s2_acceptor_atoms_ns.search(d.pos, self.DA_distance):
                if d != a:
                    pairs.update({(d,a): True})
        for a in self._s1_acceptor_atoms:
            # search for S2 donors within <distance> of acceptor a
            self._s2_donor_atoms_ns = NS.AtomNeighborSearch(self._s2_donor_atoms)
            for d in self._s2_donor_atoms_ns.search(a.pos, self.DA_distance):
                if d != a:
                    pairs.update({(d,a): True})

        logger.debug("Got %d D/A pairs" % len(pairs.keys()))
        for (d,a) in pairs.keys():
            # get hydrogen atoms for the donor
            hydrogens = self._get_bonded_hydrogens(d)
            for h in hydrogens:
                # measure the D-H-A angle and H-A distance
                if self._calc_angle(d,h,a) >= self.angle and self._calc_eucl_distance(h,a) <= self.HA_distance:
                    frame_results.append((d,h,a))
        self.timeseries.append(frame_results)

    def results(self):
        """ Returns an array containing the total count of hbonds per frame """
        return [ len(f) for f in self.timeseries ]

    def _get_bonded_hydrogens(self, atom):
        try:
            return self._hydrogens[atom]
        except:
            hydrogens = []
            for i in range(3):
                try:
                    next_atom = self.u.atoms[atom.number+1+i]
                except:
                    break
                else:
                    if next_atom.name.startswith('H'):
                        hydrogens.append(next_atom)
            self._hydrogens[atom] = hydrogens
            return hydrogens

    def _update_selections(self):
        """ Get all pairs of donors and acceptors between selections 1 and 2 """
        self._s1 = self.u.selectAtoms(self.selection1)
        self._s2 = self.u.selectAtoms(self.selection2)        
        if self.selection1_type in ('donor', 'both'):
            # selection 1 is donor
            # get all donor atoms within selection 1
            self._s1_donor_atoms = self._s1.selectAtoms(' or '.join([ 'name %s' % i for i in self.donors ]))
            # get all acceptor atoms within the selection 2
            self._s2_acceptor_atoms = self._s2.selectAtoms(' or '.join([ 'name %s' % i for i in self.acceptors ]))
        else:
            self._s1_donor_atoms = []
            
        if self.selection1_type in ('acceptor', 'both'):
            # selection 1 is acceptor
            # get all acceptor atoms within selection 1
            self._s1_acceptor_atoms = self._s1.selectAtoms(' or '.join([ 'name %s' % i for i in self.acceptors ]))
            # get all donor atoms within selection 2
            self._s2_donor_atoms = self._s2.selectAtoms(' or '.join([ 'name %s' % i for i in self.donors ]))
        else:
            self._s1_acceptor_atoms = []
    
    def _calc_angle(self, v1, v0, v2):
        """Calculate the angle (in degrees) between two atoms
        """
        v1 = v0.pos-v1.pos
        v2 = v0.pos-v2.pos
        v1 /= numpy.linalg.norm(v1)
        v2 /= numpy.linalg.norm(v2)

        if v1.tolist() == v2.tolist():
            return 0.0
        return numpy.arccos(numpy.dot(v1, v2) / (numpy.linalg.norm(v1)*numpy.linalg.norm(v2))) * 180 / numpy.pi

    def _calc_eucl_distance(self, a1, a2):
        """Calculate the Euclidean distance between two atoms.
        """
        v1 = a1.pos
        v2 = a2.pos
        return numpy.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)

# TODO: add this once hbonds is integrated into MDAnalysis
# if __name__=='__main__':
#   import optparse
#     u = MDAnalysis.Universe(PSF, PDB, permissive=True)
#     h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u, 'protein', 'resname TIP3', distance=3.0, angle=120.0)
#     results = h.run()