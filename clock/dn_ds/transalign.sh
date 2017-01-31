#!/bin/bash
#Input: FASTA file containing a set of unaligned nucleotide sequences to be aligned in a codon-aware manner.
#Output: Aligned sequences in the fasta and phylip formats (suffix "_tAlign")
#Gene tree in the Newick format based on the protein alignment.

fasta=$1
transalign=/home/sobczm/bin/transalign
clustalw=/home/sobczm/bin/clustalw1.83/clustalw
emboss=/home/armita/prog/emboss/EMBOSS-4.0.0/bin

perl $transalign/transAlign.pl -d"$fasta" -b1 -gf -fd -if -ope -mb -n20 -sn0 -ri -p$clustalw
output_align=$( echo $fasta | sed -e 's/.fasta/_mafft.fasta/' )
output_tree=$( echo $fasta | sed -e 's/.fasta/_mafft.nwk/' )
$emboss/transeq $fasta -outseq temp
mafft temp >$output_align
FastTree $output_align >$output_tree
