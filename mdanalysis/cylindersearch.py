from __future__ import with_statement

import os
import numpy
from numpy.linalg import norm

import MDAnalysis
from MDAnalysis.core.distances import self_distance_array
from MDAnalysis.core.AtomGroup import AtomGroup
import MDAnalysis.KDTree.NeighborSearch as NS

import logging
logger = logging.getLogger('cylindersearch')

class CylinderSearch(object):
    """ Count nearby residues of a given type from a given residue selection
    """
    
    def __init__(self, a, b, search, radius=10.0, extension=0.0, level='R', update_selections=True):
        """Calculate hydrogen bonds between two selections.

        :Arguments:
          *a*
            Selection for point A of cylinder core
          *b*
            Selection for point B of cylinder core
          *search*
            Selection for searchable residues/atoms within the cylinder
          *radius*
            Radius of the cylinder
          *extension*
            Extension of the cylinder to search on either side of the two points
          *level*
            R or A (Residue or Atom) for searching. Residues' center of mass are used to calculate distances.
          *update_selections*
            Update selections for points A and B at each frame
            
        The timeseries accessible as the attribute :attr:`CylinderSearch.timeseries`.
        """

        self.selection_a = a
        self.selection_b = b
        self.selection_search = search
        self.radius = radius
        self.extension = extension
        self.level = level
        self.update_selections = update_selections

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
        self.ns = NS.AtomNeighborSearch(self.u.selectAtoms(self.selection_search))
        self.timeseries = []  # final result

    def process(self, frame):
        """ Process a single trajectory frame """
        if self.update_selections:
            self._update_selections()
        
        res = []
        # find all selection within r of A and B
        near_a = set(self.ns.search(self.a, self.height+self.extension, level=self.level))
        near_b = set(self.ns.search(self.b, self.height+self.extension, level=self.level))
        near_both = near_a & near_b
        for r in near_both:
            if self.level == 'R':
                point = r.centerOfMass() # center of mass of the found residue
                name = '%s:%s' % (r.name, r.id)
            else:
                point = r.pos
                name = '%s:%s:%s' % (r.resname, r.resid, r.name)
            distance_to_vector = self._point_distance(point)
            if distance_to_vector <= self.radius:
                distance_to_a = norm(self.a-point)
                distance_to_b = norm(self.b-point)
                cylinder_offset = numpy.sqrt(distance_to_a**2-distance_to_vector**2)
                if distance_to_b > self.height:
                    cylinder_offset = -1*cylinder_offset
                res.append((name, cylinder_offset))
        self.timeseries.append(numpy.array(res))

    def results(self):
        """ Returns an array containing the total count of hbonds per frame """
        return self.timeseries
        
    def _update_selections(self):
        self.a = self.u.selectAtoms(self.selection_a).centerOfMass()
        self.b = self.u.selectAtoms(self.selection_b).centerOfMass()
        self.height = norm(self.b-self.a)
        self.vector = self.b-self.a

    def _point_distance(self, point):
        return norm(numpy.cross(point-self.a, point-self.b))/self.height
