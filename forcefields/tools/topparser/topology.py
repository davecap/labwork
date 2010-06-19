# Various topology formats, useful for converting from one to another
#
# Author: David Caplan <david@davidacaplan.com>
# Date: December 2009
#

from pyparsing import *
import string
import sys

from molecule import *

class Topology(object):

    atom_types = dict()
    residues = []
    
    def __init__(self):
        atom_types = dict()
        residues = []
        
    def summary(self):
        print "Topology Summary:"
        print " %d Atom Types" % len(self.atom_types)
        
        for res in self.residues:
            if res.__class__.__name__ == 'Residue':
                print "Residue: %s" % res
            else:
                print "Patch Residue: %s" % res
            
            print " %d Atoms" % len(res.atoms)
            print " %d Groups" % len(res.groups)
            print " %d Bonds" % len(res.bonds)
            print " %d Angles" % len(res.angles)
            print " %d Dihedrals" % len(res.dihedrals)
            print " %d Impropers" % len(res.impropers)
            
            # for g in res.groups:
            #     print g
            # print " Bonds:"
            # for b in res.bonds:
            #     print b
            # print " Angles:"
            # for b in res.angles:
            #     print b
            #print " Dihedrals:"
            #for b in res.dihedrals:
            #    print b.atom_types()
            print " Impropers:"
            for b in res.impropers:
                print b.atom_types()
    
    def read(self, infile):
        pass
    
    def write(self, outfile):
        pass

class CharmmTopology(Topology):
    rtf_bnf = None
    current_residue = None
    
    def __init__(self):
        pass
        
    def read(self, infile):
        self.current_residue = None
        self.atom_types = dict()
        self.atoms = dict()
        self.residues = []
        
        fields = self.getBNF().parseFile(infile)
        return fields
    
    def write(self, outfile):
        pass
        
    def process_parsed(self, s, loc, toks):
        # toks contains the parsed token
        # print toks
        pass

    def repeat(self, element):
        return element + ZeroOrMore(~LineEnd() + element)

    # ['MASS', ['223', 'CU', '63.546']]
    def add_atom_type(self, s, loc, toks):
        new_type = AtomType(toks[1][1], toks[1][2], toks[1][0])

        if new_type.symbol in self.atom_types:
            sys.stderr.write("Warning: Atom type duplicate found: %s\n" % new_type.symbol)
        self.atom_types[new_type.symbol] = new_type

    def find_or_create_atom_type(self, name):
        if name not in self.atom_types:
            sys.stderr.write("Warning: Referenced atom type (%s) not found, creating dummy type\n" % name)
            self.atom_types[name] = AtomType(name, '0.0', '-1')
        return self.atom_types[name]

    # ['RESI', ['CU', '2.0000']]
    def add_residue(self, s, loc, toks):
        new_residue = Residue(toks[1][0], toks[1][1])

        self.residues.append(new_residue)
        #print "Changing current residue from %s to %s" % (self.current_residue, new_residue)
        self.current_residue = new_residue

    # ['PRES', ['CHO', '0.0000']]
    def add_presidue(self, s, loc, toks):
        new_presidue = PatchResidue(toks[1][0], toks[1][1])

        self.residues.append(new_presidue)
        self.current_residue = new_presidue

    # ['GROUp', ['ATOM', ['CAD', 'C', '2.0000']]]
    # or [['ATOM', ['CAD', 'C', '2.0000']]]
    def add_group(self, s, loc, toks):
        if self.current_residue is None:
            raise "Group found without any residue!"
        atoms = []
        # first element may be the "GROUp" label, delete it if it is there
        if len(toks) > 0 and type(toks[0]) is str:
            del toks[0]
        for atom in toks:
            # get all the atom types and create new atom objects
            atom_type = self.find_or_create_atom_type(atom[1][1])
            new_atom = Atom(atom_type, atom[1][0], atom[1][2], comment=None)
            self.current_residue.atoms[new_atom.symbol] = new_atom
            atoms.append(new_atom)
        #print "Adding group to current residue %s: %s" % (self.current_residue, atoms)
        self.current_residue.groups.append(atoms)

    # ['BOND', ['FE', 'NA'], ['FE', 'NB'], ['FE', 'NC'], ['FE', 'ND'], ['NA', 'C1A']]
    def add_bond(self, s, loc, toks):
        if self.current_residue is None:
            raise "Bond found without any residue!"

        # first element is the "BOND" label
        for bond in toks[1:]:
            # get all the atom types and create new atom objects
            atom1 = self.current_residue.find_or_create_atom(bond[0])
            atom2 = self.current_residue.find_or_create_atom(bond[1])

            new_bond = Bond(atom1, atom2)
            self.current_residue.bonds.append(new_bond)

    def add_double_bond(self, s, loc, toks):
        if self.current_residue is None:
            raise "Bond found without any residue!"

        # first element is the "DOUBLE" label
        for bond in toks[1:]:
            # get all the atom types and create new atom objects
            atom1 = self.current_residue.find_or_create_atom(bond[0])
            atom2 = self.current_residue.find_or_create_atom(bond[1])

            new_bond = Bond(atom1, atom2, kind=Bond.DOUBLE)
            self.current_residue.bonds.append(new_bond)

    def add_triple_bond(self, s, loc, toks):
        if self.current_residue is None:
            raise "Bond found without any residue!"

        # first element is the "TRIPLE" label
        for bond in toks[1:]:
            # get all the atom types and create new atom objects
            atom1 = self.current_residue.find_or_create_atom(bond[0])
            atom2 = self.current_residue.find_or_create_atom(bond[1])

            new_bond = Bond(atom1, atom2, kind=Bond.TRIPLE)
            self.current_residue.bonds.append(new_bond)

    def add_aromatic_bond(self, s, loc, toks):
        if self.current_residue is None:
            raise "Bond found without any residue!"

        # first element is the "TRIPLE" label
        for bond in toks[1:]:
            # get all the atom types and create new atom objects
            atom1 = self.current_residue.find_or_create_atom(bond[0])
            atom2 = self.current_residue.find_or_create_atom(bond[1])

            new_bond = Bond(atom1, atom2, kind=Bond.AROMATIC)
            self.current_residue.bonds.append(new_bond)

    # ['ANGL', ['HCO', 'O', 'CA']]
    def add_angle(self, s, loc, toks):
        if self.current_residue is None:
            raise "Angle found without any residue!"

        # first element is the "ANGL" label
        for angle in toks[1:]:
            # get all the atom types and create new atom objects
            atom1 = self.current_residue.find_or_create_atom(angle[0])
            atom2 = self.current_residue.find_or_create_atom(angle[1])
            atom3 = self.current_residue.find_or_create_atom(angle[2])

            new_angle = Angle(atom1, atom2, atom3)
            self.current_residue.angles.append(new_angle)

    # ['DIHE', ['HCO', 'O', 'CA', 'C']]
    def add_dihedral(self, s, loc, toks):
        if self.current_residue is None:
            raise "Dihedral found without any residue!"

        # first element is the "DIHE" label
        for dihe in toks[1:]:
            # get all the atom types and create new atom objects
            atom1 = self.current_residue.find_or_create_atom(dihe[0])
            atom2 = self.current_residue.find_or_create_atom(dihe[1])
            atom3 = self.current_residue.find_or_create_atom(dihe[2])
            atom4 = self.current_residue.find_or_create_atom(dihe[3])

            new_dihe = Dihedral(atom1, atom2, atom3, atom4)
            self.current_residue.dihedrals.append(new_dihe)

    # ['IMPR', ['HCO', 'O', 'CA', 'C']]
    def add_improper(self, s, loc, toks):
        if self.current_residue is None:
            raise "Improper found without any residue!"

        # first element is the "IMPR" label
        for dihe in toks[1:]:
            # get all the atom types and create new atom objects
            atom1 = self.current_residue.find_or_create_atom(dihe[0])
            atom2 = self.current_residue.find_or_create_atom(dihe[1])
            atom3 = self.current_residue.find_or_create_atom(dihe[2])
            atom4 = self.current_residue.find_or_create_atom(dihe[3])

            new_dihe = Improper(atom1, atom2, atom3, atom4)
            self.current_residue.impropers.append(new_dihe)

    # ['IC', ['O', 'CA', '*C', 'HCO', '0.0000', '0.0000', '180.0000', '0.0000', '0.0000']]
    def add_ic(self, s, loc, toks):
        pass

    def getBNF(self):
        if self.rtf_bnf is None:
            #
            # useful types
            #

            integer = Word( nums )
            dot = Literal(".")
            float_val = Combine(Optional(Literal("-") ^ Literal("+")) + integer + dot + Optional(integer))

            name = Word(alphanums)
            text = OneOrMore(name)
            
            iupac = Combine(Optional(Literal("-") ^ Literal("+") ^ Literal("*")) + name)

            comment = Combine(Literal("!") + SkipTo(Word("\r\n")))
            EOL = (Suppress(comment) | LineEnd().suppress())
            bond = Group(iupac + iupac)
            angle = Group(iupac + iupac + iupac)
            dihedral = Group(iupac + iupac + iupac + iupac)

            #
            # keyword definitions
            #

            residue_k = CaselessKeyword("RESI") ^ CaselessKeyword("RESIdue")
            presidue_k = CaselessKeyword("PRES") ^ CaselessKeyword("PRESIdue")

            angle_k = CaselessKeyword("ANGL") ^ CaselessKeyword("ANGLe") ^ CaselessKeyword("THET") ^ CaselessKeyword("THETa")
            dihedral_k = CaselessKeyword("DIHEdral") ^ CaselessKeyword("DIHE") ^ CaselessKeyword("PHI")
            improper_k = CaselessKeyword("IMPRoper") ^ CaselessKeyword("IMPR") ^ CaselessKeyword("IMPHi") ^ CaselessKeyword("IMPH")
            donor_k = CaselessKeyword("DONOr") ^ CaselessKeyword("DONO")
            acceptor_k = CaselessKeyword("ACCEptor") ^ CaselessKeyword("ACCE")
            ic_k = CaselessKeyword("IC") ^ CaselessKeyword("BILD") ^ CaselessKeyword("BUILd") ^ CaselessKeyword("BUIL")

            decl_k = CaselessKeyword("DECLare") ^ CaselessKeyword("DECL")
            defaults_k = CaselessKeyword("DEFAults") ^ CaselessKeyword("DEFA")
            autogenerate_k = CaselessKeyword("AUTOgenerate") ^ CaselessKeyword("AUTO")
            group_k = CaselessKeyword("GROUp") ^ CaselessKeyword("GROU")

            delete_k = CaselessKeyword("DELEte") ^ CaselessKeyword("DELE")
            patch_k = CaselessKeyword("PATChing") ^ CaselessKeyword("PATC") ^ CaselessKeyword("PATCH")
            patch_first_k = CaselessKeyword("FIRSt") ^ CaselessKeyword("FIRS")
            patch_last_k = CaselessKeyword("LAST")
            
            mass_k = CaselessKeyword("MASS")
            atom_k = CaselessKeyword("ATOM")
            bond_k = CaselessKeyword("BOND")
            double_k = CaselessKeyword("DOUBLE")
            triple_k = CaselessKeyword("TRIPLE")
            aromatic_k = CaselessKeyword("AROMATIC") 
            cmap_k = CaselessKeyword("CMAP")

            #
            # line definitions
            #

            #
            # document structure definition
            #

            # NOTE: comments are suppressed!
            comment_line = Literal("!") + SkipTo(Word("\r\n"))

            # One or more header lines starting with *
            header_line = Literal("*").setParseAction(replaceWith("HEAD")) + Group(restOfLine + LineEnd().suppress())
            
            # version line
            version_line = integer + integer + Group(restOfLine + LineEnd().suppress())
            
            # one or more MASS: atom-type-code atom-type-name mass
            mass_line = mass_k + Group(integer + name + float_val + Optional(name)) + EOL            

            # DECLare out-of-residue-name
            decl_line = decl_k + Group(iupac) + EOL
            # DEFAults [ FIRSt { name } ] [ LAST { name } ]
            defaults_line = defaults_k + Group(name + name + name + name) + EOL
            # AUTOgenerate [ ANGLes ] [ DIHEdrals ]
            autogenerate_line = autogenerate_k + Group(name + name) + EOL

            # one or more RESIdue or PRESidue name [total-charge]
            residue_line = residue_k + Group(name + Optional(float_val)) + EOL
            presidue_line = presidue_k + Group(name + Optional(float_val)) + EOL

            #   one or more GROUP
            # Note: this is annoying because in the official charmm RTFs, some GROUP lines have unmarked comments:
            #   ie:         GROUP           O1  O2 (-) 
            #   instead of: GROUP       !   O1  O2 (-)
            group_line = group_k + EOL
                        
            #       one or more ATOM iupac atom-type-name charge repeat(exclusion-names)
            atom_line = atom_k + Group(iupac + iupac + float_val) + EOL

            #   one or more BOND repeat(iupac iupac)
            bond_line = bond_k + self.repeat(bond) + EOL
            
            # one or more DOUBLE repeat(iupac iupac)
            double_line = double_k + self.repeat(bond) + EOL
            triple_line = triple_k + self.repeat(bond) + EOL
            aromatic_line = aromatic_k + self.repeat(bond) + EOL

            #   one or more ANGLe or THETa repeat(iupac iupac iupac)
            angle_line = angle_k + self.repeat(angle) + EOL

            #   one or more DIHEdral or PHI repeat(iupac iupac iupac iupac)
            dihedral_line = dihedral_k + self.repeat(dihedral) + EOL

            #   one or more IMPRoper or IMPHi repeat(iupac iupac iupac iupac)
            improper_line = improper_k + self.repeat(dihedral) + EOL

            #   one or more CMAP repeat(iupac iupac iupac iupac iupac iupac iupac iupac)
            cmap_line = cmap_k + self.repeat(dihedral + dihedral) + EOL

            #   one or more DONOr [ hydrogen ] [ heavy-atom ] [ antecedent-1 antecedent-2 ]
            donor_line = donor_k + Group(iupac + iupac + Optional(~LineEnd() + iupac + iupac)) + EOL

            #   one or more ACCEptor iupac [ iupac [iupac] ]
            acceptor_line = acceptor_k + Group(iupac + Optional(~LineEnd() + iupac + Optional(~LineEnd() + iupac))) + EOL

            #   one or more IC, BILD, BUILd name name name name bond angle phi angle bond
            ic_line = ic_k + Group(dihedral + float_val + float_val + float_val + float_val + float_val) + EOL

            # DELEte   { ATOM/BOND/ANGL/DIHE/IMPR }  (iupac x1,2,3,4)  [COMBine iupac]
            delete_line = delete_k + Group(iupac + Optional(~LineEnd() + iupac + Optional(~LineEnd() + iupac + Optional(~LineEnd() + iupac)))) + EOL

            #  PATChing [ FIRSt { name } ] [ LAST { name } ]
            patch_line = patch_k + Group(patch_first_k + name + Optional(~LineEnd() + patch_last_k + name)) + EOL

            # END
            end_line = CaselessKeyword("END") + LineEnd().suppress()

            #
            # section definitions
            #

            group_section = (group_line + OneOrMore(Group(atom_line))).setParseAction(self.add_group)
            property_sections = (angle_line.setParseAction(self.add_angle) \
                                    | bond_line.setParseAction(self.add_bond) \
                                    | double_line.setParseAction(self.add_double_bond) \
                                    | triple_line.setParseAction(self.add_triple_bond) \
                                    | aromatic_line.setParseAction(self.add_aromatic_bond) \
                                    | dihedral_line.setParseAction(self.add_dihedral) \
                                    | improper_line.setParseAction(self.add_improper) \
                                    | cmap_line.setParseAction(self.process_parsed) \
                                    | donor_line.setParseAction(self.process_parsed) \
                                    | acceptor_line.setParseAction(self.process_parsed) \
                                    | ic_line.setParseAction(self.process_parsed) \
                                    | delete_line \
                                    | patch_line)

            residue_section = (residue_line.setParseAction(self.add_residue) \
                                + OneOrMore(Group(group_section)) \
                                + ZeroOrMore(Group(property_sections)))
            
            patch_residue_section = (presidue_line.setParseAction(self.add_presidue) \
                                + ZeroOrMore(Group(group_section) | Group(property_sections)) )
                                            
            optional_section = ZeroOrMore(Group(defaults_line ^ decl_line ^ autogenerate_line))

            self.rtf_bnf = OneOrMore(Group(header_line)) + version_line \
                        + OneOrMore(Group(mass_line.setParseAction(self.add_atom_type))) \
                        + optional_section \
                        + OneOrMore(Group(residue_section) | Group(patch_residue_section)) \
                        + end_line
            self.rtf_bnf.ignore(Combine(Literal("!") + SkipTo(LineEnd())))
        return self.rtf_bnf
    

class GromacsTopology(Topology):
    
    def __init__(self):
        pass
        
    def read(self, infile):
        pass
        
    def write(self, topology, outfile_prefix="out"):
        
        if len(topology.atom_types) > 0:
            # Write the atom types to the ATP
            outfile = "%s.atp" % outfile_prefix
            print "Writing ATP file %s" % outfile
            f = open(outfile, 'w')
            f.write("; Atom types from CHARMM RTF\n\n")
            for a in topology.atom_types.values():
                f.write("%s\n" % a.gromacs())
            f.close()
            print "Done writing ATP file."
        
        # Write the topology to the RTP
        outfile = "%s.rtp" % outfile_prefix
        print "Writing RTP file %s" % outfile
        
        f = open(outfile, 'w')
        
        f.write("[ bondedtypes ]\n")
        f.write("; bonds angles dihedrals impropers\n")
        f.write("1 1 1 2\n")
        f.write("\n")
        
        for res in topology.residues:
            # Skip patch residues
            if res.__class__.__name__ != 'Residue':
                continue
            
            f.write("%s\n" % res.gromacs())
            
            f.write(" [ atoms ]\n")
            f.write(" ; name type charge chargegroup\n")
            # #f.write("N N -0.280 0")
            for i, g in enumerate(res.groups):
                # iterate each atom in each group, and pass it a group index
                for a in g:
                    f.write("   %s\n" % a.gromacs(i))
            
            if len(res.bonds) > 0:
                f.write(" [ bonds ]\n")
                f.write(" ;atom1 atom2 b0 kb\n")
                #f.write("NH N\n")
                for b in res.bonds:
                    f.write("   %s\n" % b.gromacs())
            
            if len(res.angles) > 0:
                f.write(" [ angles ]\n")
                f.write(" ;atom1 atom2 atom3 th0 cth\n")
                #f.write("NH N C\n")
                for a in res.angles:
                    f.write("   %s\n" % a.gromacs())
            
            if len(res.dihedrals) > 0:
                f.write(" [ dihedrals ]\n")
                f.write(" ;atom1 atom2 atom3 atom4 phi0 cp mult\n")
                #f.write("NH N C H\n")
                for d in res.dihedrals:
                     f.write("  %s\n" % d.gromacs())
            
            if len(res.impropers) > 0:
                f.write(" [ impropers ]\n")
                f.write(" ;atom1 atom2 atom3 atom4 q0 cq\n")
                #f.write("-C -CA N -O\n")
                for i in res.impropers:
                     f.write("  %s\n" % i.gromacs())
            f.write("\n")
        f.write("\n")
        f.close()
    
