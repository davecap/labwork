#!/usr/bin/perl -w

# Copyright 2006 Mark James Abraham
#
# Version 1.2
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

# This script takes a CHARMM parameter filename as command line input and 
# produces ffcharmmbon.itp and ffcharmmnb.itp files for GROMACS as output.
# These files are then used by pdb2gmx -ff charmm to produce a .top
# file. They need to be placed by the user in the GMXLIB directory. Then
# the user should iteratively apply grompp and the fix_top_for_charmm.pl
# until grompp returns no WARNINGS - see below.
#
# CHARMM topology files are not converted by this script, or any other by
# this author. You will need to develop your own here, or use those
# available on the web from previous CHARMM-to-GROMACS conversion attempts.
# In the first instance, consult Yuguang Mu's topologies here:
# http://www.gromacs.org/old/topologies/uploaded_force_fields/charmm_gromacs.tar.gz
#
# CHARMM force field files can be obtained from various sources, including
# http://www.pharmacy.umaryland.edu/faculty/amackere/force_fields.htm
# This script was developed in conjunction with version C27. It expects
# par_all27_prot_lipid.prm or par_all27_lipid.prm as input and does
# convert the lipid terms (unlike the 1.0 version) but these have not
# been tested by the author. It is likely to work with other versions
# of the CHARMM force field, however the user is strongly encouraged to verify 
# the form and function of the force field by verifying the equivalence of 
# energy evaluations on the same structure in CHARMM and GROMACS. Be aware 
# that the switch/shift on LJ functions implemented in these two codes are
# different, and the absence of documentation of the implementation in
# CHARMM other than in the source code, and the lack of CHARMM source code
# availability to this author has meant that full reproducibility of
# the energy evaluation was impossible to test.
#
# Life is not as simple as generating new ffcharm*.itp files and then
# using pdb2gmx because some atom types have 
# different LJ parameters when they are in 1-4 interactions. These are
# indicated in ffcharmmnb.itp as commented-out seventh and eighth fields.
# grompp can't deal with this variation without serious hacking in the
# C source code. Also the default angle potential is not correct where
# a Urey-Bradley potential should be used, and it is necessary to 
# express some LJ potentials as Ryckaert-Bellemans functions rather
# than periodic sinusoidal ones, thus nearly all of them are converted
# to R-B. Accordingly, I have written fix_top_for_charmm.pl to deal
# with all of this and correct the .top file as needed. See the description
# there for details of usage and function.
#
# Revision History
#
# Version 1.0
#   Initial release, converts only protein part of par_all27_prot_lipid.prm
#   CHARMM parameter files, but possibly other versions.
#
# Changes in Version 1.1
#   Converts par_all27_prot_lipid.prm and par_all27_lipid.prm correctly,
#   but lipid parts not tested. 

# Changes in Version 1.2 
#   Converts a July 2005 par_silicates.inp successfully.

use strict;
use IO::File;
#use Data::Dumper;   # this line is needed if you want the debugging output

my $energy_conversion_factor = 4.184; # convert kcal to kJ
use constant PI => 4 * atan2 1, 1;
sub deg2rad { PI * $_[0] / 180.0 }
sub rad2deg { 180.0 * $_[0] / PI }

# This subroutine uses the two hashes below to effectively create a 
# lookup table for masses based on regular expressions that match
# element type strings in CHARMM .prm files
sub find_mass {
    my $mass_hash = shift;
    my $element_name = shift;
    my $element_mass = undef;
    foreach my $element_regexp ( keys %$mass_hash ) {
	if ( $element_name =~ $element_regexp ) {
	    $element_mass = $$mass_hash{$element_regexp};
	    last;
	}
    }
    return $element_mass;
}
    
# These hashes contain regular expressions that return the mass for the atom 
# that are intended to be passed by reference to the above subroutine. 
# %masses should be tested first so that "CLA" matches chlorine and not 
# carbon, etc.
my %masses = ( qr'^HE' => 4.003,
	       qr'^NE' => 20.180,
	       qr'^FE' => 55.847,
	       qr'^ZN' => 65.370,
	       qr'^F[NA1-3]' => 18.998,
	       qr'^SOD' => 22.990,
	       qr'^POT' => 39.102,
	       qr'^CLA' => 35.450,
	       qr'^CAL' => 40.080,
	       qr'^MG' => 24.305,
	       qr'^CES' => 132.900,
	       qr'^DUM' => 0.0,
	       qr'^AL' => 26.982
	       );
my %singlemasses = ( qr'^C' => 12.011,
		     qr'^H' => 1.008,
		     qr'^N' => 14.007,
		     qr'^O' => 15.999,
		     qr'^S' => 32.060,
		     qr'^P' => 30.974 
		     );

sub convert_dihedrals;

foreach my $filename (@ARGV) {
    my $filehandle = new IO::File($filename, 'r');
    # We need to store the information in the CHARMM .prm file before
    # writing them to the right GROMACS files, and we do this on-the-fly
    # by keeping lists of the strings already formatted for output into
    # the GROMACS file.
    my @bonds;
    my @ureybradleys;
    my @angles;
    my @dihedrals;
    my @impropers;
    my @atomtypes;

    my $line;
    my $i = 0;
    while (! eof $filehandle) {
	$line = <$filehandle>;
	# skip comments and whitespace lines
	next if ( $line =~ /^[\*\!]/ || $line =~ /^\s*$/ );
	# The next code sections depend on the ordering of the CHARMM
	# .prm file, with BONDS before ANGLES, etc.
	my $lipid_terminator = qr'^!lipid section';
	my $ignore_line = qr'^\s*([\*\!].*)?$';
	if ( $line =~ /^BONDS/ ) {
	    my $terminator = qr'^ANGLES';
	    while ($line = <$filehandle>) {
#		print STDERR $line;
		last if ( $line =~ $lipid_terminator || $line =~ $terminator );
		next if ( $line =~ $ignore_line );
		my @terms = split ' ', $line;
#		print STDERR "$line";
#		print Dumper(@terms);
		# There's an Angstrom to nm conversion that is accounted
		# for twice, and a conversion involving a kilo prefix, and
		# a stray factor of two.
		push @bonds, sprintf("% 5s % 5s % 5d %.4g %.4g\n", $terms[0], $terms[1], 1, $terms[3] / 10, $terms[2] * $energy_conversion_factor * 1000 * 2 / 10 );
#		print Dumper(@bonds);

	    }
	    while ( $line !~ $terminator ) {
		die "Ran out of input file!" if ( eof $filehandle );
		$line = <$filehandle>;
	    }
	}
	if ( $line =~ /^ANGLES/ ) {
	    # Angles are messy because in CHARMM the Urey-Bradley 
	    # angle potential is denoted only by the number of
	    # terms in the .prm file... We can tell this
	    # is a U-B angle or not by the presence of a number in the
	    # sixth element on the line.
	    my $terminator = qr'^DIHEDRALS';
	    while ($line = <$filehandle>) {
		last if ( $line =~ $lipid_terminator || $line =~ $terminator );
		next if ( $line =~ $ignore_line );
		my @terms = split ' ', $line;
#		print "$line@terms[5]\n";
#		print Dumper(@terms);
		if ( $#terms >= 5 && ord $terms[5] >= ord '0' && ord $terms[5] <= ord '9' ) {
#		    print "Found Urey-Bradley\n$line";
		    # There's an Angstrom to nm conversion that is accounted
		    # for twice, a conversion involving a kilo 
		    # prefix, and two stray factors of two.
		    push @ureybradleys, sprintf("% 5s % 5s % 5s % 5d %.4g %.4g %.4g %.4g\n", $terms[0], $terms[1], $terms[2], 5, $terms[4], $terms[3] * $energy_conversion_factor * 2, $terms[6] / 10, $terms[5] * $energy_conversion_factor * 1000 * 2 / 10);
		} else {
		    # There's just a stray factor of two to account for here.
		    push @angles, sprintf("% 5s % 5s % 5s % 5d %.4g %.4g\n", $terms[0], $terms[1], $terms[2], 1, $terms[4], $terms[3] * $energy_conversion_factor * 2);
		}
#		print Dumper(@angles);
	    }
	    while ( $line !~ $terminator ) {
		die "Ran out of input file!" if ( eof $filehandle );
		$line = <$filehandle>;
	    }
	}
#	print Dumper(@bonds);
	if ( $line =~ /^DIHEDRALS/ ) {
	    # see extensive comment below on dihedral conversions
	    my %current;
	    my @temp_dihedrals;
	    my $terminator = qr'^IMPROPER';
	    while ($line = <$filehandle>) {
		last if ( $line =~ $lipid_terminator || $line =~ $terminator );
		next if ( $line =~ $ignore_line );
		my @terms = split ' ', $line;
#		print "$line";
#		print Dumper(@terms);
		$current{'delta'} = deg2rad $terms[6];
		$current{'k'} = $terms[4] * $energy_conversion_factor;
		$current{'n'} = $terms[5];
		$current{'string'} = sprintf "% 5s % 5s % 5s % 5s ", $terms[0], $terms[1], $terms[2], $terms[3];
#		print Dumper(@dihedrals);
		push @temp_dihedrals, { %current };
	    }
	    # Magic conversion here because the CHARMM function type
	    # does not have an immediate conversion to a GROMACS type in
	    # all cases.
	    @dihedrals = convert_dihedrals \@temp_dihedrals;
	    while ( $line !~ $terminator ) {
		die "Ran out of input file!" if ( eof $filehandle );
		$line = <$filehandle>;
	    }
	}
	if ( $line =~ /^IMPROPER/ ) {
	    my $terminator = qr'^(CMAP)|(NONBONDED)';
	    while ($line = <$filehandle>) {
		last if ( $line =~ $lipid_terminator || $line =~ $terminator );
		next if ( $line =~ $ignore_line );
		my @terms = split ' ', $line;
#		print "$line";
#		print Dumper(@terms);
		# There's just a stray factor of two to account for here.
		push @impropers, sprintf("% 5s % 5s % 5s % 5s % 5d %.4g %.4g\n", $terms[0], $terms[1], $terms[2], $terms[3], 2, $terms[6], $terms[4] * $energy_conversion_factor * 2);
	    }
	    while ( $line !~ $terminator ) {
		die "Ran out of input file!" if ( eof $filehandle );
		$line = <$filehandle>;
	    }
	}
	if ( $line =~ /^CMAP/ ) {
	    my $terminator = qr'^NONBONDED';
	    while ( $line !~ $terminator ) {
		die "Ran out of input file!" if ( eof $filehandle );
		$line = <$filehandle>;
	    }
	}
#	print Dumper(@impropers);
	if ( $line =~ /^NONBONDED/ ) {
	    my $terminator = qr'^(HBOND)|(NBFIX)|(END)';
	    while ($line = <$filehandle>) {
		last if ( $line =~ $lipid_terminator || $line =~ $terminator );
		next if ( $line =~ $ignore_line || $line =~ /^cutnb/ );
		my @terms = split ' ', $line;
#		print "$line";
#		print Dumper(@terms);
		# It was quite difficult to work out how to convert the
		# LJ terms from one to the other, because the documentation
		# for both CHARMM and GROMACS was unclear, and they both
		# use different cut-off schemes that made comparing the
		# results of the calculations on the same structure
		# a near-futile exercise. The conversions implied below
		# seem to allow the reproduction of CHARMM in GROMACS,
		# however.
		my $four_epsilon = $terms[2] * $energy_conversion_factor;
		my $sigma = $terms[3] * 2 / 10;
		my $flag_14 = 0;
		my ($four_epsilon_14, $sigma_14);
		if ( $#terms >= 4 && $terms[4] ne '!') {
		    # we have special 1-4 vdW terms
		    $flag_14 = 1;
		    $four_epsilon_14 = $terms[5] * $energy_conversion_factor;
		    $sigma_14 = $terms[6] * 2 / 10;
		}
		# search for the element mass from the name and then
		# store the string suitable for outputing one line for
		# this atom type
		foreach my $mass_hash ( \%masses, \%singlemasses ) {
		    my $element_mass = find_mass $mass_hash, $terms[0];
		    if ( defined $element_mass ) {
			if ( $flag_14 ) {
			    push @atomtypes, sprintf("% 5s %f   0.0     A %.4g %.4g ; %.4g %.4g\n", $terms[0], $element_mass, 2*$four_epsilon * $sigma**6, $four_epsilon * $sigma**12, 2*$four_epsilon_14 * $sigma_14**6, $four_epsilon_14 * $sigma_14**12 );
			} else {
			    push @atomtypes, sprintf("% 5s %f   0.0     A %.4g %.4g\n", $terms[0], $element_mass, 2*$four_epsilon * $sigma**6, $four_epsilon * $sigma**12);
			}
			last;
		    } elsif ( $mass_hash == \%singlemasses ) {
			die "Didn't find mass for atom type $terms[0]\n";
		    }
		}

#		print Dumper(@atomtypes);
	    }
	    while ( $line !~ $terminator ) {
		die "Ran out of input file!" if ( eof $filehandle );
		$line = <$filehandle>;
	    }
	}
#	print Dumper(@atomtypes);
    }
    my $ffbon = new IO::File("ffcharmmbon.itp", 'w');
    print $ffbon "[ bondtypes ]\n ; i    j    func   b0      kb\n";
    foreach my $line ( @bonds ) { print $ffbon $line };
    print $ffbon "\n[ angletypes ] ; normal potential harmonic in theta\n  ;  i    j    k   func    th0         kth\n";
    foreach my $line ( @angles ) { print $ffbon $line };
    print $ffbon "[ angletypes ] ; Urey-Bradley potentials\n ; i     j     k   func    th0         kth      b0       kb\n";
    foreach my $line ( @ureybradleys ) { print $ffbon $line };
    print $ffbon "\n[ dihedraltypes ]\n  ;  i    j    k     l   func    phi0     kphi      n\n";
    foreach my $line ( @dihedrals ) { print $ffbon $line };
    print $ffbon "; Now the improper dihedrals\n";
    foreach my $line ( @impropers ) { print $ffbon $line };
    undef $ffbon;
    my $ffnb = new IO::File("ffcharmmnb.itp", 'w');
    print $ffnb "[ atomtypes ]\n;name mass     charge ptype     c6       c12\n";
    foreach my $line ( @atomtypes ) { print $ffnb $line };
    undef $ffbon;
}

sub construct_normal_dihedral {
    my $atoms = shift;
    my $phi = shift;
    my $k = shift;
    my $n = shift;
    return sprintf "%s% 5d %.4g %.4g %d\n", $atoms, 1, $phi, $k, $n;
}

sub construct_RB_dihedral {
    my $atoms = shift;
    my $RB_list = shift;
    # need to change the dihedral convention going from periodic to RB in
    # GROMACS (see 4.2 of GROMACS manual) by applying RB[i] *= (-1)^i
    foreach my $value ( $$RB_list[1], $$RB_list[3], $$RB_list[5] ) {
	$value = 0-$value;
    }
    return sprintf "%s% 5d %.4g %.4g %.4g %.4g %.4g %.4g\n", $atoms, 3, $$RB_list[0], $$RB_list[1], $$RB_list[2], $$RB_list[3], $$RB_list[4], $$RB_list[5];
}

# We need some elaborate functionality to convert the CHARMM dihedral type
# of k * (1 + cos(n * xi - delta ) ) functions summed over n into something
# GROMACS can implement. While the above functional form exists in
# GROMACS, you can't have more than one function of this type, and
# CHARMM has a number of dihedral interactions that require more than
# one such function. However for delta = 0 or pi and n <= 5, then the above
# cosine function can be expanded in powers of cos xi, and the coefficients
# of the expansion can be summed in this conversion and presented to
# GROMACS as a ready-made Ryckaert-Bellemans dihedral. In practice, this
# works because CHARMM only uses such delta and n values for atom type 
# combinations that need multiple functions of the above form. Warnings
# are issued when delta is some other value, and the algorithm dies if
# n is > 6. In order to simplify GROMACS logfile output so that it only
# has to report one sort of dihedral term for most simulations, all 
# dihedral terms with n <= 5 are expressed as R-B, even when not necessary.
# Dihedrals with n=6 are left in periodic form, since it is not possible
# to convert these to R-B form when the summation is limited to the
# fifth power of cos xi.

sub convert_dihedrals {
    my $dihedrals = shift;
    my %dihedral_hash;
    my @new_dihedrals;
    push @new_dihedrals, "; First come the dihedrals with sin delta != 0 and those with\n;  multiplicity >= 6, then the remainder, converted to Ryckaert-Bellemans\n;  form, followed by comments of the contributing underlying periodic form(s).\n";
    foreach my $current ( @$dihedrals ) {
#	print Dumper($current);
	# never test a floating point number against a floating point
	# number, always allow a small error range
	if ( abs sin $$current{'delta'} > 1e-12 ) {
	    print "Sin of delta not zero, hope this is the only dihedral of\n type $$current{'string'}\n (sin delta = " . sin($$current{'delta'}) . ")\n";
	    push @new_dihedrals, construct_normal_dihedral($$current{'string'}, rad2deg($$current{'delta'}), $$current{'k'}, $$current{'n'});
	    next;
	}
	if ( ! defined $dihedral_hash{$$current{'string'}} && $$current{'n'} != 6 ) {
	    $dihedral_hash{$$current{'string'}} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
	}
	my $RB_list = $dihedral_hash{$$current{'string'}};
#	print Dumper($RB_list);
#	print Dumper($$RB_list[0]);
	my $cosphi = cos $$current{'delta'};
	if ( $$current{'n'} == 1 ) {
	    $$RB_list[0] += $$current{'k'};
	    $$RB_list[1] += $$current{'k'}*$cosphi;
	} elsif ( $$current{'n'} == 2 ) {
	    $$RB_list[0] += $$current{'k'}*(1-$cosphi);
	    $$RB_list[2] += 2*$$current{'k'}*$cosphi;
	} elsif ( $$current{'n'} == 3 ) {
	    $$RB_list[0] += $$current{'k'};
	    $$RB_list[1] += -3*$$current{'k'}*$cosphi;
	    $$RB_list[3] += 4*$$current{'k'}*$cosphi;
	} elsif ( $$current{'n'} == 4 ) {
	    $$RB_list[0] += $$current{'k'}*(1+$cosphi);
	    $$RB_list[2] += -8*$$current{'k'}*$cosphi;
	    $$RB_list[4] += 8*$$current{'k'}*$cosphi;
	} elsif ( $$current{'n'} == 5 ) {
	    $$RB_list[0] += $$current{'k'};
	    $$RB_list[1] += 5*$$current{'k'}*$cosphi;
	    $$RB_list[3] += -20*$$current{'k'}*$cosphi;
	    $$RB_list[5] += 16*$$current{'k'}*$cosphi;
	} elsif ( $$current{'n'} == 6 ) {
	    push @new_dihedrals, construct_normal_dihedral($$current{'string'}, rad2deg($$current{'delta'}), $$current{'k'}, $$current{'n'});
	} else {
	    print "Not implemented dihedral of multiplicity $$current{'n'}\n";
	}
    }
    foreach my $key ( keys %dihedral_hash ) {
#	my $RB_list = $dihedral_hash{$key};
	push @new_dihedrals, construct_RB_dihedral($key, $dihedral_hash{$key});
	foreach my $current ( @$dihedrals ) {
	    if ( $$current{'string'} eq $key ) {
		# also output the equivalent "normal dihedrals" that were
		# summed to produce this R-B dihedral to help the user
		push @new_dihedrals, '; ' . substr(construct_normal_dihedral($$current{'string'}, rad2deg($$current{'delta'}), $$current{'k'}, $$current{'n'}), 2);
	    }
	}
    }
    return @new_dihedrals;
}
