#!/usr/bin/perl

use strict;

my %imiss;
my %removed;

my $threshold = $ARGV[1];
my $input_file = $ARGV[0];

open IMISS, '<', $ARGV[0].".imiss"
        or die "Cannot open genotypes file (".$ARGV[0].".imiss): $!\n";
while(<IMISS>){
	s/^\s+//;
    my @fields = split /\s+/, $_;
    $imiss{$fields[0]}{$fields[1]} = $fields[5];
}

open GENOME, '<', $ARGV[0].".genome"
        or die "Cannot open genotypes file (".$ARGV[0].".genome): $!\n";
while(<GENOME>){
    s/^\s+//;
    my @fields = split /\s+/, $_;
 	if($fields[11] > $threshold){
 		if($imiss{$fields[0]}{$fields[1]}>$imiss{$fields[2]}{$fields[3]}){
 			unless($removed{$fields[0]}{$fields[1]}){
 				print "$fields[0] $fields[1]\n";
 				$removed{$fields[0]}{$fields[1]} = 1;
 			}
 		}
 		elsif($imiss{$fields[0]}{$fields[1]}<$imiss{$fields[2]}{$fields[3]}){
 			unless($removed{$fields[2]}{$fields[3]}){
 				print "$fields[2] $fields[3]\n";
 				$removed{$fields[2]}{$fields[3]} = 1;
 			}
 		}
 		else{
 			unless($removed{$fields[0]}{$fields[1]}){
 				print "$fields[0] $fields[1]\n";
 				$removed{$fields[0]}{$fields[1]} = 1;
 			}
 		}
 	}
}
    
	

