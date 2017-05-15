#! /usr/bin/perl -w
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::SearchIO::Writer::HSPTableWriter;
$contigs = $ARGV[0];
$db = $ARGV[1];
$seed = int(rand(10000000));
$format = "\"6 std qlen slen nident qframe sframe sstrand\"";

my $in = Bio::SeqIO->new(-file => $contigs, -format => 'Fasta');

while ($seq = $in->next_seq())
 {
		$id = $seq->id;
		$sequence = $seq->seq;
		open (TEMP, '>', "temp_$seed.fasta");
		print TEMP "\>$id\n$sequence\n";
		close TEMP;
		system "blastx -outfmt $format -num_threads 1 -max_target_seqs 1 -evalue 0.0000000001 -query temp_$seed.fasta -db $db >> $contigs\_vs\_$db";
}
