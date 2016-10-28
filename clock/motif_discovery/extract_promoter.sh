#!/bin/bash
input=/home/sobczm/popgen/clock/DNA_genomes/bedtools

gff_file=Trichoderma_atroviride_imi_206040.TRIAT_v2_0.32.gff3

grep "##sequence-region" $gff_file | cut -d ' ' -f4,5,6 | sed 's/ /\t/g' >Tatroviride.scaffolds

genome=Tatroviride.scaffolds
grep -E "#|CDS" Trichoderma_atroviride_imi_206040.TRIAT_v2_0.32.gff3 >Trichoderma_atroviride_imi_206040.TRIAT_v2_0.32_CDS.gff3
cds_gff_file=Trichoderma_atroviride_imi_206040.TRIAT_v2_0.32_CDS.gff3
bedtools sort -i $cds_gff_file > Trichoderma_atroviride_imi_206040.TRIAT_v2_0.32_CDS_sorted.gff3
bedtools flank -i Trichoderma_atroviride_imi_206040.TRIAT_v2_0.32_CDS_sorted.gff3 -g $genome -l 2000 -r 0 >select_intervals
dna=Trichoderma_atroviride_imi_206040.TRIAT_v2_0.dna_rm.toplevel.fa
bedtools getfasta -fi $dna -bed select_intervals -fo genes.2kb.promoters.bed.fa


#The GFF3 (?) format used by Ensembl is prepared as follows.
#There are multiple CDS with the same ID for a given gene and they
#Correspond to exons minus the 5' and 3' UTR regions.

#So use CDS for getting out the upstream region of the gene. For a given protein_id,
#extract the CDS with the lowest coordinate number (sort by coordinates first) and take 3000 bp u
#upstream from that, for positive strand.
#For negative strand, do so for the highest coordinate number.


#For Fus2, use Andy's gff2fasta and extract 1, 2, 3 kbp upstream:
input=/home/sobczm/popgen/clock/DNA_genomes
scripts=/home/sobczm/bin/popgen/clock/motif_discovery
cd $input
$scripts/gff2fasta.pl Fus2_canu_contigs_hardmasked.fa Fus2_final_genes_appended.gff3 outFus2
#And change values in the script to extract 1 and 2 kbp upstream.
