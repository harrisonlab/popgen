#!/bin/bash
input=/home/sobczm/popgen/clock/DNA_genomes/bedtools

gff_file=Trichoderma_atroviride_imi_206040.TRIAT_v2_0.32.gff3
cds_gff_file=Trichoderma_atroviride_imi_206040.TRIAT_v2_0.32_CDS.gff3
cds_gff_file_sorted=Trichoderma_atroviride_imi_206040.TRIAT_v2_0.32_CDS_sorted.gff3
genome=Tatroviride.scaffolds
dna=Trichoderma_atroviride_imi_206040.TRIAT_v2_0.dna_rm.toplevel.fa

grep "##sequence-region" $gff_file | cut -d ' ' -f4,5,6 | sed 's/ /\t/g' >$genome
grep -E "#|CDS" $gff_file >$cds_gff_file

bedtools sort -i $cds_gff_file > $cds_gff_file_sorted
bedtools flank -i $cds_gff_file_sorted -g $genome -l 2000 -r 0 >select_intervals
bedtools getfasta -fi $dna -bed select_intervals -fo genes.2kb.promoters.bed.fa


#The GFF3 (?) format used by Ensembl is prepared as follows.
#There are multiple CDS with the same ID for a given gene and they
#correspond to exons minus the 5' and 3' UTR regions.

#So use CDS for getting out the upstream region of the gene. For a given protein_id,
#extract the CDS with the lowest coordinate number (sort by coordinates first) and take 3000 bp
#upstream from that, for positive strand.
#For negative strand, do so for the highest coordinate number.

#Important, bedtools extracts the sequence for negative strand, as it is, without
#reverse complement necessary for scan of the promoter region!


#For Fus2, use Andy's gff2fasta and extract 1, 2, 3 kbp upstream:
input=/home/sobczm/popgen/clock/DNA_genomes
scripts=/home/sobczm/bin/popgen/clock/motif_discovery
cd $input
$scripts/gff2fasta.pl Fus2_canu_contigs_hardmasked.fa Fus2_final_genes_appended.gff3 outFus2
#And change values in the script to extract 1 and 2 kbp upstream.
