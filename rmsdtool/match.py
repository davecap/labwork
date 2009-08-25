#!/usr/bin/env python

"""
Set of routines to calculate the RMSD between two molecular structures.
The module can be run from the command line using PDB files as input.

Input:
    - set of PDB files
    - region of the protein
Output:
    - top residues that have high RMSDs

"""

DEFAULT_THRESHOLD = 1

import math
import numpy
import vector3d, util, molecule, polymer

import Bio.PDB
import xpdb

def rmsd(crds1, crds2):
  """Returns RMSD between 2 sets of [nx3] numpy array"""

  assert(crds1.shape[1] == 3)
  assert(crds1.shape == crds2.shape)

  n_vec = numpy.shape(crds1)[0]
  correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
  v, s, w = numpy.linalg.svd(correlation_matrix)
  is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w)) < 0.0
  if is_reflection:
    s[-1] = - s[-1]
  E0 = sum(sum(crds1 * crds1)) + \
       sum(sum(crds2 * crds2))
  rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
  rmsd_sq = max([rmsd_sq, 0.0])
  return numpy.sqrt(rmsd_sq)


def optimal_superposition(crds1, crds2):
  """Returns best-fit rotation matrix as [3x3] numpy matrix"""
  assert(crds1.shape[1] == 3)
  assert(crds1.shape == crds2.shape)
  correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
  v, s, w = numpy.linalg.svd(correlation_matrix)
  is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w)) < 0.0
  if is_reflection:
    v[-1,:] = -v[-1,:]
  return numpy.dot(v, w)


def get_i_residue(residues, tag):

  def get_tag(residue):
    tag = ""
    if residue.chain_id != " " and residue.chain_id != "":
      tag += residue.chain_id + ":"
    tag += str(residue.num)
    if residue.insert:
      tag += residue.insert
    return tag  

  # clean up tag
  tag = tag.strip()
  if tag[0] == ":":
    tag = tag[1:]
  if not tag[0].isdigit() and tag[1].isdigit():
    tag = tag[0] + ":" + tag[1:]

  for i, residue in enumerate(residues):
    if tag.lower() == get_tag(residue).lower():
      return i
  raise "Can't find residue", tag
  
def get_superposable_atoms(polymer, segments, atom_types=['CA', 'N', 'C', 'CB']):
    result = []
    allowed_i = []
    residues = polymer.residues()
    #if no segments provided, take the whole backbone
    if len(segments) == 0:
        print "no segments given, using backbone"
        for i, residue in enumerate(residues):
            result.extend([a for a in residue.atoms() if a.type in atom_types])
    else:
        print segments
        for res_num_i, res_num_j in segments:
            i = get_i_residue(residues, str(res_num_i))
            j = get_i_residue(residues, str(res_num_j))
            allowed_i.extend(range(i,j))
        for i, residue in enumerate(residues):
            if i in allowed_i:
                result.extend([a for a in residue.atoms() if a.type in atom_types])
    return result

def get_crds(atoms):
  crds = numpy.zeros((len(atoms), 3), float)
  for i, a in enumerate(atoms):
    crds[i,0] = a.pos.x
    crds[i,1] = a.pos.y
    crds[i,2] = a.pos.z
  return crds


def calculate_superposition_matrix(atoms1, atoms2):

  def convert_to_matrix3d(numpy_matrix3d):
    result = vector3d.Matrix3d()
    for i in range(3):
      for j in range(3):
        result.setElem(i, j, numpy_rotation[j, i])
    return result

  numpy_rotation = optimal_superposition(get_crds(atoms1), get_crds(atoms2))
  return convert_to_matrix3d(numpy_rotation)
    

def sum_rmsd(atoms1, atoms2):
  sum_squared = 0.0
  for atom1, atom2 in zip(atoms1, atoms2):
    sum_squared += vector3d.pos_distance(atom1.pos, atom2.pos)**2
  return math.sqrt(sum_squared/float(len(atoms1)))

def get_best_alignment(pdb1, pdb2, segments, atom_types):
  """Returns rmsd and filename of transformed pdb2."""
  polymer1 = polymer.Polymer(pdb1)
  atoms1 = get_superposable_atoms(polymer1, segments, atom_types)
  polymer2 = polymer.Polymer(pdb2)
  atoms2 = get_superposable_atoms(polymer2, segments, atom_types)

  center1 = molecule.get_center(atoms1)
  polymer1.transform(vector3d.Translation(-center1))
  polymer2.transform(vector3d.Translation(-molecule.get_center(atoms2)))
  polymer2.transform(calculate_superposition_matrix(atoms1, atoms2))

  rmsd = sum_rmsd(atoms1, atoms2)
  
  temp_pdb2 = util.fname_variant(pdb2)
  polymer2.transform(vector3d.Translation(center1))
  polymer2.write_pdb(temp_pdb2)
  
  return rmsd, temp_pdb2

def get_raw_rmsd(pdb1, pdb2, segments, atom_types):
    polymer1 = polymer.Polymer(pdb1)
    atoms1 = get_superposable_atoms(polymer1, segments, atom_types)
    polymer2 = polymer.Polymer(pdb2)
    atoms2 = get_superposable_atoms(polymer2, segments, atom_types)
    
    return sum_rmsd(atoms1, atoms2)

def get_rmsd(pdb1, pdb2, segments, atom_types):
    polymer1 = polymer.Polymer(pdb1)
    atoms1 = get_superposable_atoms(polymer1, segments, atom_types)
    polymer2 = polymer.Polymer(pdb2)
    atoms2 = get_superposable_atoms(polymer2, segments, atom_types)

    #do alignment
    center1 = molecule.get_center(atoms1)
    polymer1.transform(vector3d.Translation(-center1))
    polymer2.transform(vector3d.Translation(-molecule.get_center(atoms2)))

    crds1 = get_crds(atoms1)
    crds2 = get_crds(atoms2)
    return rmsd(crds1, crds2)


def segments_str(segments):
  residues = []
  for i, j in segments:
    if i == j:
      residues.append(str(i))
    else:
      residues.append("%s-%s" % (i,j))
  return ', '.join(residues)
  

def biopython_atom_to_vector3d(atom):
    return vector3d.Vector3d(atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2])

#passing biopython Atom class
def bp_sum_rmsd(residue1, residue2, atom_types=['CA', 'N', 'C', 'CB', 'CG']):
    atoms1 = [ a for a in residue1.get_list() if a.get_name() in atom_types ]
    atoms2 = [ a for a in residue2.get_list() if a.get_name() in atom_types ]
    #convert Atom to vector3d
    #v3d_atoms1 = [ biopython_atom_to_vector3d(a) for a in atoms1 ]
    #v3d_atoms2 = [ biopython_atom_to_vector3d(a) for a in atoms2 ]
    
    sum_squared = 0.0
    for atom1, atom2 in zip(atoms1, atoms2):
        sum_squared += vector3d.pos_distance(biopython_atom_to_vector3d(atom1), biopython_atom_to_vector3d(atom2))**2
    return math.sqrt(sum_squared/float(len(atoms1)))

def residue_key(residue):
    return "%s:%s:%s" % (residue.get_segid().strip(), residue.get_resname(), residue.get_id()[1])

def get_target_residues_from_structure(structure, segids):
    residues = []
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                #NOTE: there is whitespace in the .get_segid() value
                if residue.get_segid().strip() in segids:
                    residues.append(residue)
                    #print residue.get_id()[1]
                    #print residue.get_resname()
                    #print residue.get_segid()
    return residues

if __name__ == '__main__':
    
    import sys, os, optparse
    
    usage = """
        %prog [options] <PDB 0> <PDB 1> <PDB 2> ... <PDB N>"
    
    Copyright (c) 2007 Bosco Ho
    Modified by David Caplan, 2009
    
    Calculate the RMSD between two or more structures.
    
    PDB 0 is the starting structure. All following structures (PDB 1 ... PDB N) are compared to it.
    
    Optional segments can be specified. This will perform a simple analysis on those residues and rank them by RMSD.

    segments: a string that encodes the residues to be matched.
    
    e.g, "[('A:5', 'A:10'), ('B:3', 'B:19')]". For
    convenience, the ":" character is optional, and quotes are not
    needed if there are no chain identifiers. Put insertions at the end 
    of the residue tag "A:335E"

    """
    
    parser = optparse.OptionParser(usage)
    
    parser.add_option("-r", action="store_true", dest="no_rotation", default=False,
                        help="Calculates direct RMSD without any rotations [default: %default]")
    parser.add_option("-n", action="store_true", dest="no_save", default=True,
                        help="Calculates RMSD without saving rotated structure [default: %default]")
    parser.add_option("-s", "--segments", dest="segments", default="[]",
                        help="Segments to be compared [default: whole backbone]")
    parser.add_option("-t", "--threshold", dest="threshold", default=DEFAULT_THRESHOLD, 
                        help="Default RMSD threshold [default: %default]")
    options, args = parser.parse_args()

    #if options.index_range and len(args) != 2:
    #    parser.error("When using an indx range, only <start> <end> are required as arguments")

    if len(args) < 2:
        parser.error("Provide some PDB files")
        
    pdb0 = args[0]
    pdb1_to_n = args[1:]

    segments = eval(options.segments)
    # segments1 = eval(args[2])
    # s = segments_str(segments1)
    # if len(args) > 3:
    #     segments2 = eval(args[3])
    #     s += 'to ' + segments_str(segments2)
    # else:
    #     segments2 = segments1
    #     print "Aligning CA atoms of residues:", s

    # open PDB0 file

    # read
    pdb0_structure = xpdb.get_structure(pdbid='pdb0', pdbfile=pdb0)
    pdb0_residues = get_target_residues_from_structure(pdb0_structure, segids=["A"])
    rmsd_data = { }
    for r in pdb0_residues:
        #print "%s:%s:%s" % (r.get_segid().strip(), r.get_resname(), r.get_id()[1])
        rmsd_data[residue_key(r)] = { 'pdb0_residue': r, 'rmsds': [] }
    
    # write PDB file
    # sloppyio = xpdb.SloppyPDBIO()
    # sloppyio.set_structure(structure)
    # sloppyio.save('new_big_fat.pdb')
    
    for i, pdbi in enumerate(pdb1_to_n):
        pdbi_structure = xpdb.get_structure(pdbid='pdb%d'%i, pdbfile=pdbi)
        pdbi_residues = get_target_residues_from_structure(pdbi_structure, segids=["A"])
        
        for i, r in enumerate(pdbi_residues):
            if not rmsd_data.has_key(residue_key(r)):
                raise "Can't find residue key %s! Make sure structures have the same residues." % residue_key(r)
            else:
                rmsd = bp_sum_rmsd(rmsd_data[residue_key(r)]['pdb0_residue'], r)
                #print "residue: %s, rmsd: %f" % (residue_key(r), rmsd)
                #bp_sum_rmsd(residue1, residue2, atom_types=['CA', 'N', 'C', 'CB', 'CG']):
                rmsd_data[residue_key(r)]['rmsds'].append(rmsd)
    
    for r, d in rmsd_data.iteritems():
        narray = numpy.array(d['rmsds'])
        print "%s -> min: %f max: %f mean: %f" % (r, narray.min(), narray.max(), narray.mean())
    
    # if options.no_rotation:
    #     print "No rotations"
    #     rmsd = get_raw_rmsd(pdb0, pdb1_to_n[0], segments, ['CA'])
    # elif options.no_save:
    #     rmsd = get_rmsd(pdb0, pdb1_to_n[0], segments, ['CA'])
    # else:
    #     rmsd, temp_pdb = get_best_alignment(pdb0, pdb1_to_n[0], segments, ['CA'])
    #     print "Optimal superposition of %s written to: %s" % (pdb1_to_n[0], temp_pdb)
    #print "RMSD: %.3f" % rmsd
