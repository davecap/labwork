# Classes for molecule topology components
#
#
# Author: David Caplan <david@davidacaplan.com>
# Date: December 2009

import sys

# Atoms types:
#   <symbol>: <mass>, [unique index]
class AtomType(object):
    symbol = ''
    mass = -1
    index = -1
    comment = ''
    
    def __init__(self, symbol, mass, index=0, comment=None):
        self.symbol = symbol
        self.mass = float(mass)
        self.index = int(index)
        self.comment = comment

    def __repr__(self):
        return self.symbol

    def __str__(self):
        return self.symbol

    def charmm(self):
        # MASS atom-type-code atom-type-name mass
        if self.index < 0:
            sys.stderr.write('Warning: Atom type %s has no unique index' % self.symbol)
        return 'MASS %d %s %f %s' % (self.index, self.symbol, self.mass, self.comment)
    
    def gromacs(self):
        output = '%s  %f' % (self.symbol, self.mass)
        if self.comment:
            output = '%s    ;%s' % (output, self.comment)
        return output
    

class DummyAtomType(AtomType):
    def __init__(self):
        self.symbol = 'DUM'
        self.mass = -9999
        self.index = 9999
        self.comment = 'Dummy Atom Type'

class Atom(object):
    atom_type = None
    symbol = ''
    charge = 0.0
    comment = None
    
    def __init__(self, atom_type, symbol, charge, comment=None):
        self.atom_type = atom_type
        self.symbol = symbol
        self.charge = float(charge)
        self.comment = comment
    
    def gromacs(self, group=0):
        # #f.write("N N -0.280 0")
        
        output = '%s %s %f %d' % (self.symbol, self.atom_type, self.charge, group)
        if self.comment:
            output = '%s    ;%s' % (output, self.comment)
        return output
    
    def __repr__(self):
        return "%s" % (self.symbol)
            
    def __str__(self):
        return "%s" % (self.symbol)


class DummyAtom(Atom):
    def __init__(self, symbol):
        self.symbol = symbol
        self.charge = -9999.9999
        self.comment = 'Dummy Atom'
        self.atom_type = DummyAtomType()

class Bond(object):
    atoms = []
    comment = None
    
    # bond types
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 4
    
    def __init__(self, atom1, atom2, kind=1, comment=None):
        self.atoms = [atom1, atom2]
        self.kind = kind
        self.comment = comment
    
    def gromacs(self):
        output = '%s %s' % (self.atoms[0], self.atoms[1])
        if self.comment:
            output = '%s    ;%s' % (output, self.comment)
        return output
    
    def __repr__(self):
        return "(%s, %s)" % (self.atoms[0], self.atoms[1])

    def __str__(self):
        return "(%s, %s)" % (self.atoms[0], self.atoms[1])

class Angle(object):
    atoms = []
    comment = None

    def __init__(self, atom1, atom2, atom3, comment=None):
        self.atoms = [atom1, atom2, atom3]
        self.comment = comment

    def gromacs(self):
        output = '%s %s %s' % (self.atoms[0], self.atoms[1], self.atoms[2])
        if self.comment:
            output = '%s    ;%s' % (output, self.comment)
        return output

    def __repr__(self):
        return "(%s, %s, %s)" % (self.atoms[0], self.atoms[1], self.atoms[2])

    def __str__(self):
        return "(%s, %s, %s)" % (self.atoms[0], self.atoms[1], self.atoms[2])

class Dihedral(object):
    atoms = []
    
    def __init__(self, atom1, atom2, atom3, atom4, comment=None):
        self.atoms = [atom1, atom2, atom3, atom4]
        self.comment = comment
    
    def gromacs(self):
        output = '%s %s %s %s' % (self.atoms[0], self.atoms[1], self.atoms[2], self.atoms[3])
        if self.comment:
            output = '%s    ;%s' % (output, self.comment)
        return output
    
    def __repr__(self):
        return "(%s, %s, %s, %s)" % (self.atoms[0], self.atoms[1], self.atoms[2], self.atoms[3])
    
    def __str__(self):
        return "(%s, %s, %s, %s)" % (self.atoms[0], self.atoms[1], self.atoms[2], self.atoms[3])
        
class Improper(Dihedral):
    pass

class CMAP(object):
    # CMAPs contain 2 dihedrals
    dihedrals = []
    
    def __init__(self, dihedral1, dihedral2):
        self.dihedrals = [dihedral1, dihedral2]
        
    def __repr__(self):
        return "(%s, %s)" % (self.dihedrals[0], self.dihedrals[1])
    


class Donor(object):
    hydrogen = None
    heavy_atom = None
    antecedent1 = None
    antecedent2 = None
    
    def __init__(self, hydrogen, heavy_atom, antecedent1=None, antecedent2=None):
        self.hydrogen = hydrogen
        self.heavy_atom = heavy_atom
        self.antecedent1 = antecedent1
        self.antecedent2 = antecedent2
    
class Acceptor(object):
    atom1 = None
    atom2 = None
    # Unused
    atom3 = None

    def __init__(self, atom1, atom2, atom3=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3

# Internal Coordinates (CHARMM)
#{ BILD or BUILd } name name name name bond angle phi angle bond

class IC(object):
    dihedral = None
    properties = []
    
    def __init__(self, dihedral, properties):
        self.dihedral = dihedral
        self.properties = properties

# Residue types:
class Residue(object):
    name = ''
    total_charge = 0.0
    comment = ''
    
    groups = []
    atoms = dict()
    bonds = []
    angles = []
    dihedrals = []
    impropers = []
    cmaps = []
    donors = []
    acceptors = []
    ics = []
    
    def __init__(self, name, total_charge, comment=None):
        self.name = name
        self.total_charge = total_charge
        self.comment = comment
        self.groups = []
        self.bonds = []
        self.atoms = dict()
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        self.cmaps = []
        self.donors = []
        self.acceptors = []
        self.ics = []

        if comment:
            self.comment = '! %s' % comment

    def __str__(self):
        return self.name
    
    def gromacs(self):
        output = '[ %s ]' % (self.name)
        if self.comment:
            output = '%s    ;%s' % (output, self.comment)
        return output
    
    def find_or_create_atom(self, symbol):
        try:
            atom = self.atoms[symbol]
            return atom
        except KeyError:
            sys.stderr.write("Warning: Could not find atom by symbol: %s, using dummy atom\n" % symbol)
            return DummyAtom(symbol)

class PatchResidue(Residue):
    pass