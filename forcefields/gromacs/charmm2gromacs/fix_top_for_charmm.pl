#!/usr/bin/perl -w

# Copyright 2006 Mark James Abraham
#
# Version 1.1
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

use strict;
use IO::File;
#use Data::Dumper;   # this line is needed if you want the debugging output

# This script prepares a .top file produced by "pdb2gmx -ff charmm" for
# actual use as input to grompp. It requires that grompp be used 
# iteratively - one grompp is done by this script in order to produce 
# the incorrect .top file, and a bunch of WARNINGS from which 
# this script generates the corrected .top file. The user then should 
# run grompp again to check the output for the absence of WARNINGS.

# In the indicated .top file, the default pdb2gmx angle potentials 
# are with replaced with Urey-Bradley potentials where they are 
# appropriate, the normal periodic dihedrals are replaced with 
# Ryckaert-Bellemans dihedrals and 1-4 LJ parameters are included
# for those interactions that involve atoms that have such 
# parameters specialized for 1-4 interactions.

# This script requires the command line arguments you would give to
# grompp, which are passed through unmodified. ffcharmmnb.itp needs
# to be in either the directory indicated by the GMXLIB 
# environment variable, or in its absence, the current working
# directory. The output will be in the file indicated by the -p
# command line argument, or in its absence "topol.top". Any previous
# copy of this file will be backed up in the normal GROMACS manner.

# Revision History
#
# Version 1.0
#   Initial release
#

print "grompp @ARGV -maxwarn 999 2>&1 1>/dev/null\n";
my @output = `grompp @ARGV -maxwarn 999 2>&1 1>/dev/null`;

my @warnings = grep (/(^WARNING)|(^  No default)/, @output);
#print Dumper(@warnings);

if ( $#warnings == -1 ) {
    die "Success! .top file did not need fixing for CHARMM forcefield.\n";
}

# first parse ffcharmnb.itp to get the atom types that need special 1-4
# interactions
my $nb_filename = (defined $ENV{GMXLIB} ? "$ENV{GMXLIB}/" : "./") . 'ffcharmmnb.itp';
unless ( -r $nb_filename ) {
    die "Couldn't find ffcharmnb.itp as $nb_filename!";
}
my $nb_file = new IO::File($nb_filename, 'r');
my %special_LJ_params;
my %LJ_params;
while ( my $line = <$nb_file> ) {
    #skip commented lines
    next if ( $line =~ /^[\[;]/ );
    my @words = split ' ', $line;
    if ( $#words >= 6 ) {
        print @words;
	# we've got a line with the special 1-4 parameters
#	print "Got $line\n";
	$special_LJ_params{$words[0]} = [ $words[7], $words[8] ];
    }
    # and keep these because we need them too
    $LJ_params{$words[0]} = [ $words[4], $words[5] ];
}
undef $nb_file;
#print Dumper(%special_LJ_params);

"@ARGV" =~ /-p (\S+) /;
my $filename = $1;
if ( ! defined $filename ) {
    $filename = 'topol.top';
}

# The next two lines generate the filename of the backup .top that
# GROMACS would generate if this were a GROMACS utility overwriting
# a .top file, i.e. #topol.top.32# if there were already 31 of them.
my @filenames = glob "#$filename\.*#";
my $saved_filename = "#$filename." . ($#filenames+2) . "#";
#print Dumper($saved_filename);
#exit;
`cp $filename $saved_filename`;
my $input = new IO::File($saved_filename, 'r');
my $output = new IO::File($filename, 'w');

my $i = 0;
my $line_to_seek;
# We rely on the grompp warnings being in numerical order
if ( $warnings[$i] =~ /line (\d+)\]/ ) {
    $line_to_seek = $1;
} else {
    $line_to_seek = -1;
}
my $line_number = 1;
my ($inatoms, $inpairs) = (0, 0);
my @atom_types  = undef;
while ( my $line = <$input> ) {
#    print "\$line_to_seek $line_to_seek warning $warnings[$i]\n";
    if ( $line_number == $line_to_seek && $i != $#warnings ) {
	$i++;
#	printf "\$i %d \$#warnings %d\n", $i, $#warnings;
	if ( $warnings[$i] =~ /^  No default (.*) types/ ) {
	    my $interaction_type = $1;
	    if ( $interaction_type eq 'Angle' ) {
		$line =~ s/^(\s+\d+\s+\d+\s+\d+\s+)(\d+)/${1}5/;
		print STDERR "Replaced Angle with U-B\n";
	    } elsif ( $interaction_type eq 'Ryckaert-Bell.' ) {
		$line =~ s/^(\s+\d+\s+\d+\s+\d+\s+\d+\s+)(\d+)/${1}1/;
		print STDERR "Replaced periodic with R-B\n";
	    } else {
		print STDERR "Didn't recognise interaction type '$interaction_type' on warning:\n$warnings[$i]";
	    }
	    $i++;
	    if ( $i <= $#warnings ) {
		if ( $warnings[$i] =~ /line (\d+)\]/ ) {
		    $line_to_seek = $1;
		} else {
		    $line_to_seek = -1;
		}
	    } else {
		$line_to_seek = -1;
	    }
	} else {
	    $i++;
	    $line_to_seek = -1;
	}
    }
    if ( $inatoms ) {
	if ( $line =~ /^\[/ ) {
	    $inatoms--;
	} else {
	    if ( $line !~ /^;/ && $line =~ /\S/ ) {
		my @words = split ' ', $line;
		push @atom_types, $words[1];
	    }
	}
    }
    if ( $inpairs ) {
	if ( $line =~ /^\[/ ) {
	    $inpairs--;
	} else {
	    if ( $line !~ /^;/ && $line =~ /\S/ ) {
		my @words = split ' ', $line;
#		print Dumper($atom_types[$words[0]], $atom_types[$words[1]]);
		my ($first_LJ, $second_LJ);
		if ( defined $special_LJ_params{$atom_types[$words[0]]} ) {
#		    print "Found special LJ on first atom\n";
		    $first_LJ = $special_LJ_params{$atom_types[$words[0]]};
		} else {
		    $first_LJ = $LJ_params{$atom_types[$words[0]]};
		}
		if ( defined $special_LJ_params{$atom_types[$words[1]]} ) {
#		    print "Found special LJ on second atom\n";
		    $second_LJ = $special_LJ_params{$atom_types[$words[1]]};
		} else {
		    $second_LJ = $LJ_params{$atom_types[$words[1]]};
		}
		chop $line;
#		print STDERR Dumper($$first_LJ[0], $$second_LJ[0]);
		$line .= sprintf " %e %e ; %s %s\n", sqrt($$first_LJ[0] * $$second_LJ[0]), sqrt($$first_LJ[1] * $$second_LJ[1]), $atom_types[$words[0]], $atom_types[$words[1]];
	    }
	}
    }
    # these are deliberately placed after the tests for inatoms and
    # in pairs
    if ( $line =~ /^\[ atoms/ ) {
	$inatoms++;
    }
    if ( $line =~ /^\[ pairs/ ) {
	$inpairs++;
    }
    print $output $line;
    $line_number++;
}
#print STDERR Dumper(%LJ_params);
