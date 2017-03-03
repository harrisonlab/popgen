#!/usr/local/bin/perl  -w

open (RESULTS, '>', "codeml_table_nonfixed.txt");
open (RESULTS_FIX, '>', "codeml_table_fixed.txt");
print RESULTS "locus\tt\tN\tS\tdN\tdS\tdNdS\tlnL\n";
print RESULTS_FIX "locus\tt\tN\tS\tdN\tdS\tdNdS\tlnL\n";

@codemla= <*.out>;
foreach my $codeml (@codemla) 
{
	open (FILE, "$codeml");
	while (<FILE>)
	{
		$codeml =~ m/^([a-zA-Z0-9_]+)/;
		$locus = $1;
		next if 1 .. /pairwise comparison, codon frequencies: Fcodon./;
		$line = $_;
		chomp $line;
		if ($line =~ m/\S+/)
		{
			if ($line =~ m/lnL/)
			{
				@temp = split(/=/, $line);
				$likelihood = $temp[1];
			}
			
			
			if ($line =~ m/dN/)
			{
			@temp = split(/\s+/,$line);
			$time = $temp[1];
			$s = $temp[3];
			$n = $temp[5];
			$dnds = $temp[7];
			$dn = $temp[10];
			$ds = $temp[13];
				
			}
		
		}
	
	}
	
	close FILE;
	unless ($codeml =~ m/fixed/)
	{ 
		print RESULTS "$locus\t$time\t$n\t$s\t$dn\t$ds\t$dnds\t$likelihood\n";
	}
	else 	
	{
		print RESULTS_FIX "$locus\t$time\t$n\t$s\t$dn\t$ds\t$dnds\t$likelihood\n";
	}
}
	
close RESULTS;
close RESULTS_FIX;