#!/usr/bin/perl
use strict;
use warnings;

my $sum = 0;
my $N = 0;
my @X;
my $SS = 0;
my $avg;
while( my $line = <stdin> ){
     chomp $line;
     next if ($line =~ m/^[ \#]/);
     last if($line eq "");
    #print $line;
     push(@X, $line);
     $N ++;
     $sum += $line;
}
$avg = $sum/$N;
foreach my $X (@X){
     $SS += ($X - $avg)**2;
}
my $stdev = sqrt($SS/($N-1));
print $stdev/sqrt($N);
print "\n";
