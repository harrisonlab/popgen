#!/usr/local/bin/perl  -w

#Parses codeml output for pairwise Ka/Ks calculation in a 
#given folder, and calculates LRT based
#on comparing to the model fixed at omega=0. LRT >3.841
#implies significance at the 0.05 threshold in one-tailed test.

#Parsed output file: codeml_table.txt

open (RESULTS, '>', "codeml_table.txt");
print RESULTS "locus\tt\tN\tS\tdN\tdS\tdNdS\tlnL\tLRT\n";

@codemla= <*.out>;
@codemla = grep !/fixed/, @codemla;
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
			@temp = split(/=/,$line);
			foreach my $n (@temp)
			{ $n =~ s/[^0-9\.]//g;  }
			$time = $temp[1];
			$s = $temp[2];
			$n = $temp[3];
			$dnds = $temp[4];
			$dn = $temp[5];
			$ds = $temp[6];
				
			}
		
		}
	
	}
	
	close FILE;

	$codeml =~ s/\.out/_fixed\.out/;
	#print $codeml;
	
	open (FILE2, $codeml);
	while (<FILE2>)
	{
		next if 1 .. /pairwise comparison, codon frequencies: Fcodon./;
		$line = $_;
		chomp $line;
		#print $line;
		if ($line =~ m/lnL/)
			{
				#print "aaaa";
				@temp = split(/=/, $line);
				$likelihood2 = $temp[1];
				$lrt = -2 * ($likelihood2 - $likelihood)
			}
	}
	#print $likelihood2;
	close FILE2;

	print RESULTS "$locus\t$time\t$n\t$s\t$dn\t$ds\t$dnds\t$likelihood\t$lrt\n";
	

}
	
close RESULTS;
