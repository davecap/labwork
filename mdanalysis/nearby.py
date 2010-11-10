# Nearby Count Analysis
from __future__ import with_statement

import os
import numpy

import MDAnalysis
from MDAnalysis.core.distances import self_distance_array
from MDAnalysis.core.AtomGroup import AtomGroup
import MDAnalysis.KDTree.NeighborSearch as NS

import logging
logger = logging.getLogger('nearbycount')

class NearbyCountAnalysis(object):
    """ Count nearby residues of a given type from a given residue selection
    """
    
    def __init__(self, selection1='protein', selection2='all', cutoff=3.0):
        """Calculate hydrogen bonds between two selections.

        :Arguments:
          *selection1*
            Selection string for first selection
          *selection2*
            Selection string for second selection
          *cutoff*
            Distance cutoff
            
        The timeseries accessible as the attribute :attr:`NearbyCountAnalysis.timeseries`.
        """

        self.selection1 = selection1
        self.selection2 = selection2
        self.cutoff = cutoff
        
        if not (self.selection1 and self.selection2):
            raise Exception('NearbyCountAnalysis: invalid selections')

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

    def process(self, frame):
        """ Process a single trajectory frame """
        self._s2_ns = NS.AtomNeighborSearch(self._s2)
        res = self._s2_ns.search_list(self._s1, self.cutoff)
        logger.debug("Got %d nearby" % len(res))
        self.timeseries.append(res)
        return res

    def results(self):
        """ Returns an array containing the total count of hbonds per frame """
        return [ len(f) for f in self.timeseries ]

    def _update_selections(self):
        self._s1 = self.u.selectAtoms(self.selection1)
        self._s2 = self.u.selectAtoms(self.selection2)
