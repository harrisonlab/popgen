#!/usr/bin/perl -w

use Bio::SeqIO;
use Cwd;

#Get current directory
my $dir = getcwd;
#Find all subdirs with Busco runs
my @dirs = grep { -d } glob "$dir/run*";

#Read in all the target single copy BUSCO loci to be extracted
open (LIST, "align_input_list.txt");

while (<LIST>)
{
	$id = $_;
	chomp $id;
	$ls{$id} = 1;
}
	
close LIST;

#Enter all subdirs with BUSCO runs to load the genome assembly FASTA sequences into a hash
for my $d (@dirs)
{
	chdir($d) or die "Failed to chdir to $dir: $!";
	%cdna=();
	undef @fasta;
	undef @full_table;
	
	my @fasta = <*.fa*>;
	my @full_table = <full_table*>;
	
	$seqin = Bio::SeqIO->new(-file => "<$fasta[0]", '-format' => 'Fasta');
	while (my $seq = $seqin->next_seq)
	{
		$locus = $seq->id;
		$dna = $seq->seq;
		$cdna{$locus} = $dna;
	}

# Loop over the results table to find the assembly gene ids for BUSCO single copy genes, look them up in the hash and save their FASTA
# sequences to individual files

	open (FT, $full_table[0]);
	while (<FT>)
	{
		my $line = $_;
		chomp $line;
		my @lines = split (/\t+/, $line); 	
		# BUSCO ID
		$a = $lines[0];
		# Assembly gene ID
		$b = $lines[2];
		
		if (exists $ls{$a})
		{
			$fhs = "$dir/$a.fasta";
			open (FS, ">>$fhs");
			print FS (">$b\n$cdna{$b}\n");
			close FS;
		}
	}
	
	close FT;
}
