#! /usr/bin/perl

## To use this script, type : perl convertfasta2nex.pl fastafilename > nexusfilename
# S. Hedtke Sept. 2011
## rev. Jul. 2013

my $filename=$ARGV[0];
chomp $filename;
my $numtaxa=`grep -c ">" $filename`;
chomp $numtaxa;
my $numchar=0;
my @seqs;
my @species;
my $seq;
my $maxlength;


open (DATA,$filename);

while (my $t=<DATA>) {
	chomp $t;
	if ($t=~/^\>/) { 
		chomp $t; 
		my @poo=split(/\>/,$t); 
		unless ($species[0] eq '') {push @seqs,$seq;} 
		push @species,$poo[1]; $seq=''; $numchar=0;
		next; }
	else { $numchar=$numchar+length $t; $seq=$seq.$t; if ($maxlength<$numchar) {$maxlength=$numchar;} }

}
push @seqs,$seq;
print "#NEXUS\nBegin data;\ndimensions ntax=$numtaxa nchar=$maxlength;\nformat datatype=dna missing=-;\nmatrix\n";

for (my $i=0;$i<scalar@species;$i++) {
	my $len=length$seqs[$i];
	if ($len<$maxlength) { my $add=$maxlength-$len; for (my $j=0; $j<$add; $j++) {$seqs[$i]=$seqs[$i].'-';}}
	print "\n$species[$i]\t$seqs[$i]";
}

print "\n;\nend;";
exit;

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#   See <http://www.gnu.org/licenses/>.

