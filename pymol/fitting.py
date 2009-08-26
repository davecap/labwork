#! /usr/bin/env python
# Copyright (c) 2005 Robert L. Campbell

from pymol import cmd

def fitting(obj1,select1,obj2,select2):
  """
DESCRIPTION

  "fitting" allows the superpositioning of object1 onto object2 using
  the atoms in selection1 and selection2 (side chains are ignored).  The
  residue names, residue numbers chain identifiers, segment identifiers,
  and alt ids of selection1 are changed to match those in selection2,
  temporarily.  This allows the normal "fit" command to work.  They are
  reset after "fit" is run and two new selections are created showing
  the selected atoms.

  Be careful when creating your selection strings. Within the
  selections, do not include the object name because the chain, residue
  name, residue number etc. of selection1 of object1 are converted to
  match those in selection2.  If the object names are included in the
  selections, no atoms will be selected since an atom cannot exist in
  both object1 and object2 at the same time.

  It is important that the beginning residue numbers specify the
  aligned residues, but the ending numbers are not critical.  The
  shorter of the two selections is used in the fit calculation.

USAGE 

  fitting object1, selection1, object2, selection2

  DO NOT include object names in selections!

EXAMPLES

  fitting 1xuu, c. a & (i. 296-309 or i. 335-340), 1ame, i. 8-21 or i. 47-52

  """

  list_m = []
  list_n = []

  backbone = 'n. n+ca+c+o &! r. hoh+wat'
  select1 = '(%s) & %s' % (select1,backbone)
  select2 = '(%s) & %s' % (select2,backbone)
  m=cmd.get_model("%s & %s" % (obj1,select1))
  n=cmd.get_model("%s & %s" % (obj2,select2))

# for the atoms to be used in fit:
# store id, chain, resn, resi, name, segi, alt
  for at in m.atom:
    list_m.append((at.id,at.chain,at.resn,at.resi,at.name,at.segi, at.alt))
  for at in n.atom:
    list_n.append((at.id,at.chain,at.resn,at.resi,at.name,at.segi, at.alt))

  if len(m.atom) <= len(n.atom):
    total = len(m.atom)
  else:
    total = len(n.atom)

# set a new segi for the atoms to be used in fit command and to allow resetting later
  seg_fit="1fit"

# change the chain,resn,resi,segi and alt of select1 to match select2
  for i in range(total):
    cmd.do("alter %s & id %s, chain='%s'" % (obj1,list_m[i][0],list_n[i][1]))
    cmd.do("alter %s & id %s, resn='%s'" % (obj1,list_m[i][0],list_n[i][2]))
    cmd.do("alter %s & id %s, resi=%s" % (obj1,list_m[i][0],list_n[i][3]))
    cmd.do("alter %s & id %s, segi='%s'" % (obj1,list_m[i][0],seg_fit))
    cmd.do("alter %s & id %s, alt='%s'" % (obj1,list_m[i][0],list_n[i][6]))
# change the segid for obj2 and select2
    cmd.do("alter %s & id %s, segi='%s'" % (obj2,list_n[i][0],seg_fit))

  print "Fitting %s and %s\n     to %s and %s" % (obj1,select1,obj2,select2)

  print "Altered to:"
  print "%s & %s & segi %s\n" % (obj1,select2,seg_fit),
  print "%s & %s & segi %s\n" % (obj2,select2,seg_fit),
  print "--------------------------------------------\n"
  rms = cmd.fit("%s & %s & segi %s" % (obj1,select2,seg_fit),"%s & %s & segi %s" % (obj2,select2,seg_fit) ,quiet=0)

  cmd.delete("%s_fitting" % obj1)
  cmd.delete("%s_fitting" % obj2)
# create new objects to show the fit atoms
  cmd.create("%s_fitting" % obj1, "%s & %s & segi %s" % (obj1,select2,seg_fit))
  cmd.create("%s_fitting" % obj2, "%s & %s & segi %s" % (obj2,select2,seg_fit))

# reset chain,resn,resi,segi & alt of obj1 & select1 from stored list
  for atoms_m in list_m:
    cmd.do("alter %s & id %s, chain='%s'" % (obj1,atoms_m[0],atoms_m[1]))
    cmd.do("alter %s & id %s, resn='%s'" % (obj1,atoms_m[0],atoms_m[2]))
    cmd.do("alter %s & id %s, resi=%s" % (obj1,atoms_m[0],atoms_m[3]))
    cmd.do("alter %s & id %s, segi='%s'" % (obj1,atoms_m[0],atoms_m[5]))
    cmd.do("alter %s & id %s, alt='%s'" % (obj1,atoms_m[0],atoms_m[6]))
# reset segi of obj2 & select2 from stored list
  for atoms_n in list_n:
    cmd.do("alter %s & id %s, segi='%s'" % (obj2,atoms_n[0],atoms_n[5]))

  print "RMSD for fitting selection %s of %s onto \n                 selection %s of %s = %6.3f" % (select1, obj1, select2, obj2, rms)

cmd.extend("fitting",fitting)
