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

import math
import numpy
import vector3d, util, molecule, polymer


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
  
  
def get_superposable_atoms(polymer, segments, 
           atom_types=['CA', 'N', 'C', 'CB']):
  result = []
  allowed_i = []
  residues = polymer.residues()
  for res_num_i, res_num_j in segments:
    i = get_i_residue(residues, str(res_num_i))
    j = get_i_residue(residues, str(res_num_j))
    allowed_i.extend(range(i,j))
  for i, residue in enumerate(residues):
    if i in allowed_i:
      result.extend([a for a in residue.atoms()
                     if a.type in atom_types])
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
  

def get_raw_rmsd(pdb1, pdb2, segments1, segments2, atom_types):
  polymer1 = polymer.Polymer(pdb1)
  polymer2 = polymer.Polymer(pdb2)
  atoms1 = get_superposable_atoms(polymer1, segments1, atom_types)
  atoms2 = get_superposable_atoms(polymer2, segments2, atom_types)
  return sum_rmsd(atoms1, atoms2)


def get_best_alignment(pdb1, pdb2, segments1, segments2, atom_types):
  """Returns rmsd and filename of transformed pdb2."""
  polymer1 = polymer.Polymer(pdb1)
  atoms1 = get_superposable_atoms(polymer1, segments1, atom_types)
  polymer2 = polymer.Polymer(pdb2)
  atoms2 = get_superposable_atoms(polymer2, segments2, atom_types)

  center1 = molecule.get_center(atoms1)
  polymer1.transform(vector3d.Translation(-center1))
  polymer2.transform(vector3d.Translation(-molecule.get_center(atoms2)))
  polymer2.transform(calculate_superposition_matrix(atoms1, atoms2))

  rmsd = sum_rmsd(atoms1, atoms2)
  
  temp_pdb2 = util.fname_variant(pdb2)
  polymer2.transform(vector3d.Translation(center1))
  polymer2.write_pdb(temp_pdb2)
  
  return rmsd, temp_pdb2


def get_rmsd(pdb1, pdb2, segments1, segments2, atom_types):
  polymer1 = polymer.Polymer(pdb1)
  atoms1 = get_superposable_atoms(polymer1, segments1, atom_types)
  polymer2 = polymer.Polymer(pdb2)
  atoms2 = get_superposable_atoms(polymer2, segments2, atom_types)

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
  


if __name__ == '__main__':

  import sys, os, getopt
  
  opts, args = getopt.getopt(sys.argv[1:], "nrs")
  flags = [o for o,a in opts]

  usage = """

  Copyright (c) 2007 Bosco Ho

  Calculates the CA rmsd between 2 PDB structures, and generates the
  optimal superposition of the 2nd structure to the 1st

  Usage: python match.py [-nrs] pdb1 pdb2 segments1 [segments2]

  -n Calculates direct RMSD without any rotations

  -r Calculates RMSD without saving rotated structure

  -s Show aligned structures with pymol

  segments1: a string that encodes the residues to be matched from
  pdb1. if segments2 is not given, it is assumed that the residues
  listed in segments1 will be used for pdb2

  segments2: string encoding residues from pdb2

  format of the string: e.g, "[('A:5', 'A:10'), ('B:3', 'B:19')]". For
  convenience, the ":" character is optional, and quotes are not
  needed if there are no chain identifiers. Put insertions at the end 
  of the residue tag "A:335E"

  """

  if len(args) < 3:
    print usage 
    sys.exit(0)

  pdb1 = args[0]
  pdb2 = args[1]

  segments1 = eval(args[2])
  s = segments_str(segments1)
  if len(args) > 3:
    segments2 = eval(args[3])
    s += 'to ' + segments_str(segments2)
  else:
    segments2 = segments1
  print "Aligning CA atoms of residues:", s

  if '-n' in flags:
    print "No rotations"
    rmsd = get_raw_rmsd(pdb1, pdb2, segments1, segments2, ['CA'])
  elif '-r' in flags:
    rmsd = get_rmsd(pdb1, pdb2, segments1, segments2, ['CA'])
  else:
    rmsd, temp_pdb = get_best_alignment(pdb1, pdb2, segments1, segments2, ['CA'])
    print "Optimal superposition of %s written to: %s" % (pdb2, temp_pdb)

  print "RMSD: %.3f" % rmsd

  if '-s' in flags:
    # insert your own command here
    # cmd = 'show.py %s %s > /dev/null &' % (pdb1, temp_pdb)
    # os.system(cmd)
    pass







