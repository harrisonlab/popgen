#!/usr/bin/perl

# Converts an aligned fasta (aa or dna) seq file to phylip format

my $usage = "Usage: $0 [-h] [-v] [-c numChar] [infile]\n" .
    "  -h: help\n" .
    "  -c: long seq names are shortened to 10 char, default: 7 chars from the\n".
    "      beggining is combined with the last 3 chars.  You can change the\n".
    "      behavior by this option.  E.g., -c 10 will shorten the long name\n" .
    "      by taking the first 10 characters of the name.\n".
    "  -v:  verbose (print name conversion in STDERR)\n" .
    " infile should be an aligned fasta, " .
    "STDIN is used if no infile is given\n";

use IO::File;
use Getopt::Std;
getopts('hvc:') || die "$usage\n";

die "$usage\n" if (defined ($opt_h));

my $totNumChar = 100;  # number of characters allowed for name in phylip
my $numFrontChar = 100; # When the name is too long, this amount of characters
                      # are used from the beginning of the name, and the rest
                      # are from the end of the name.
if (defined ($opt_c)) {
    if ($opt_c <= $totNumChar && $opt_c >= 0) {
	$numFrontChar = $opt_c;
    } else {
	die "Error: give an integer between 0 and $totNumChar (ends inclusive)".
	    " for -c.\n";
    }
}

my $tmpFH = IO::File->new_tmpfile || die "Unable to make new temp file: $!";

my $firstLine = 1;
my $maxLen = 0;
my $numTaxa = 0;
my $name;

while(<>){

    chomp;
    s/^\s+//; s/\s$//;
    next if (/^$/);

    if (s/^>\s*//) {

	if ($firstLine == 0) {
	    if ($seqLen != $maxLen && $maxLen != 0) {
		warn "WARN: The $numTaxa-th species ($name) have ",
		     "different seq length\n";
		warn "Previous Seq Len: $maxLen, This Seq Len: $seqLen\n";
	    }
	    print $tmpFH "\n";    # end of the previous sequence
	} else {
	    $firstLine = 0;
	}

	$maxLen = $seqLen if ($seqLen > $maxLen); $seqLen = 0;
	$numTaxa ++;

	$name = $_;
	if (CharLen($_) <=10) {
	    printf $tmpFH "%-10s ", $_;
	} else  {
	    $shortName = TrimName($_);
	    print STDERR "$_ => $shortName\n" if (defined ($opt_v));
	    printf $tmpFH "%10s ", $shortName;
	}
    } else {
	$seqLen += CharLen ($_);
	print $tmpFH $_;
    }
}

print $tmpFH "\n";

### print out to the STDOUT
print "$numTaxa $maxLen\n";

seek ($tmpFH, 0, 0) || die "seek: $!";
my $line;
while (defined ($line = $tmpFH->getline())) {
    chomp ($line);
    print "$line";
    #$missingBases = $maxLen - (CharLen($line) - $totNumChar);
    print "\n";
}

sub CharLen {  # returns number of characters in a string
    my $string = shift;
    return scalar(split (//, $string));
}

sub TrimName { # trim a long name
    my $name = shift;
    my @charArray = split (//, $name);

    if ($totNumChar == $numFrontChar) {
	return join ('', splice (@charArray, 0, $numFrontChar));
    }  else {
	return join ('', splice (@charArray, 0, $numFrontChar),
		     splice (@charArray, - ($totNumChar - $numFrontChar)));
    }
}
    

