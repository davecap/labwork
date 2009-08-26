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

# def segments_str(segments):
#   residues = []
#   for i, j in segments:
#     if i == j:
#       residues.append(str(i))
#     else:
#       residues.append("%s-%s" % (i,j))
#   return ', '.join(residues)
  
# Biopython code
  
def biopython_atom_to_vector3d(atom):
    return vector3d.Vector3d(atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2])

#passing biopython Atom class
def bp_sum_rmsd(residue1, residue2, atom_types=['CA', 'N', 'C', 'CB', 'CG']):
    atoms1 = [ a for a in residue1.get_list() if a.get_name() in atom_types ]
    atoms2 = [ a for a in residue2.get_list() if a.get_name() in atom_types ]
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
    parser.add_option("-s", "--segments", dest="segments", default="['A']",
                        help="Segments to be compared [default: whole backbone]")
    parser.add_option("-t", "--threshold", dest="threshold", default=DEFAULT_THRESHOLD, 
                        help="Default RMSD threshold [default: %default]")
    options, args = parser.parse_args()

    if len(args) < 2:
        parser.error("Provide some PDB files")
        
    pdb0 = args[0]
    pdb1_to_n = args[1:]

    segments = eval(options.segments)
    
    # read PDB0
    pdb0_structure = xpdb.get_structure(pdbid='pdb0', pdbfile=pdb0)
    pdb0_residues = get_target_residues_from_structure(pdb0_structure, segids=segments)
    rmsd_data = { }
    for r in pdb0_residues:
        #print "%s:%s:%s" % (r.get_segid().strip(), r.get_resname(), r.get_id()[1])
        rmsd_data[residue_key(r)] = { 'pdb0_residue': r, 'rmsds': [] }
    
    # write PDB file
    # sloppyio = xpdb.SloppyPDBIO()
    # sloppyio.set_structure(structure)
    # sloppyio.save('new_big_fat.pdb')
    
    for i, pdbi in enumerate(pdb1_to_n):
        #print "Comparing pdb%d" % i
        pdbi_structure = xpdb.get_structure(pdbid='pdb%d'%i, pdbfile=pdbi)
        pdbi_residues = get_target_residues_from_structure(pdbi_structure, segids=segments)
        
        for i, r in enumerate(pdbi_residues):
            if not rmsd_data.has_key(residue_key(r)):
                raise "Can't find residue key %s! Make sure structures have the same residues." % residue_key(r)
            else:
                rmsd = bp_sum_rmsd(rmsd_data[residue_key(r)]['pdb0_residue'], r)
                #print "residue: %s, rmsd: %f" % (residue_key(r), rmsd)
                #bp_sum_rmsd(residue1, residue2, atom_types=['CA', 'N', 'C', 'CB', 'CG']):
                rmsd_data[residue_key(r)]['rmsds'].append(rmsd)
    
    print "RESIDUE,MIN,MAX,MEAN"
    
    for r, d in rmsd_data.iteritems():
        narray = numpy.array(d['rmsds'])
        print "%s,%f,%f,%f" % (r, narray.min(), narray.max(), narray.mean())
    
    # if options.no_rotation:
    #     print "No rotations"
    #     rmsd = get_raw_rmsd(pdb0, pdb1_to_n[0], segments, ['CA'])
    # elif options.no_save:
    #     rmsd = get_rmsd(pdb0, pdb1_to_n[0], segments, ['CA'])
    # else:
    #     rmsd, temp_pdb = get_best_alignment(pdb0, pdb1_to_n[0], segments, ['CA'])
    #     print "Optimal superposition of %s written to: %s" % (pdb1_to_n[0], temp_pdb)
    #print "RMSD: %.3f" % rmsd
