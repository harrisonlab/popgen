#!/bin/bash
input=/home/sobczm/popgen/clock/DNA_genomes
scripts=/home/sobczm/bin/popgen/clock/motif_discovery

#Get V. dahliae JR2 genome, extract promoters at 1 kbp, 2 kbp and 3 kbp using extract_promoter.sh.
cd $input/JR2
wget  ftp://ftp.ensemblgenomes.org/pub/fungi/release-34/fasta/verticillium_dahliaejr2/dna/Verticillium_dahliaejr2.GCA_000400815.2.dna_rm.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-34/gff3/verticillium_dahliaejr2/Verticillium_dahliaejr2.GCA_000400815.2.34.gff3.gz

#Reverse complement the promoter sequences to get the orientation - START codon.
for a in *promoters*.fasta
do
/home/armita/prog/emboss/EMBOSS-4.0.0/bin/revseq $a ${a%.fasta}_revcom.fasta
done

#Look for the He (2005) motif with max. 30 bp gap between the two parts of the motif.
qsub $scripts/sub_glam2scan.sh Verticillium_dahliaejr2.GCA_000400815.2.dna_rm.toplevel_promoters_1000.fasta $input/promoters/extended/he2005_extended.txt
qsub $scripts/sub_glam2scan.sh Verticillium_dahliaejr2.GCA_000400815.2.dna_rm.toplevel_promoters_2000.fasta $input/promoters/extended/he2005_extended.txt
qsub $scripts/sub_glam2scan.sh Verticillium_dahliaejr2.GCA_000400815.2.dna_rm.toplevel_promoters_3000.fasta $input/promoters/extended/he2005_extended.txt