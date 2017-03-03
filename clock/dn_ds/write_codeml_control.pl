#!/usr/local/bin/perl  -w

@phylipa= <*.phylip>;
foreach my $phyl (@phylipa) 
{

$phyl =~ m/^([a-zA-Z0-9_]+)/;
$locus = $1;
open (OUT, '>', "$locus.ctl");  
print OUT "seqfile = $phyl\n";
print OUT "outfile = $locus.out\n";
print OUT "noisy = 0\n";
print OUT "verbose = 0 \n";
print OUT "runmode = -2\n";
print OUT "seqtype = 1\n";
print OUT "CodonFreq = 3\n";
print OUT "model = 0\n";
print OUT "NSsites = 0\n";
print OUT "icode = 0\n";
print OUT "fix_kappa = 0\n";
print OUT "kappa = 1\n";
print OUT "fix_omega = 0\n";
print OUT "omega = 0.5\n";
print OUT "cleandata = 1\n";
close OUT;


$phyl =~ m/^([a-zA-Z0-9_]+)/;
$locus = $1;
open (OUT2, '>', "$locus\_fixed.ctl");  
print OUT2 "seqfile = $phyl\n";
print OUT2 "outfile = $locus\_fixed.out\n";
print OUT2 "noisy = 0\n";
print OUT2 "verbose = 0 \n";
print OUT2 "runmode = -2\n";
print OUT2 "seqtype = 1\n";
print OUT2 "CodonFreq = 3\n";
print OUT2 "model = 0\n";
print OUT2 "NSsites = 0\n";
print OUT2 "icode = 0\n";
print OUT2 "fix_kappa = 0\n";
print OUT2 "kappa = 1\n";
print OUT2 "fix_omega = 1\n";
print OUT2 "omega = 1\n";
print OUT2 "cleandata = 1\n";
close OUT2;

}