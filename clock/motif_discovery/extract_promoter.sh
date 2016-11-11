#!/bin/bash
input=/home/sobczm/popgen/clock/DNA_genomes
scripts=/home/sobczm/bin/popgen/clock/motif_discovery
#cd $input
for gff_file in *.32.gff3
do
dna="${gff_file%.32.gff3}.dna_rm.toplevel.fa"
genome="${gff_file%.32.gff3}"
gff_file_sorted="${gff_file%.gff3}_sorted.gff3"
#Cut from gff3 file description to get a bedtools input file with scaffolds' length
grep "##sequence-region" $gff_file | cut -d ' ' -f4,5,6 | sed 's/ /\t/g' >$genome
bedtools sort -i $gff_file > $gff_file_sorted
#loop over desired promoter region sizes to be extracted
interval=( 1000 2000 )
    for n in "${interval[@]}"
    do
#The GFF3 (?) format used by Ensembl is prepared as follows.
#There are multiple CDS with the same ID for a given gene and they
#correspond to exons minus the 5' and 3' UTR regions.
#So use CDS for getting out the upstream region of the gene. For a given protein_id,
#extract the CDS with the lowest coordinate number (sort by coordinates first) and take 3000 bp
#upstream from that, for positive strand.
#For negative strand, do so for the highest coordinate number.
    pos_strand="${gff_file_sorted%.gff3}_plus_strand_$n.gff3"
    neg_strand="${gff_file_sorted%.gff3}_minus_strand_$n.gff3"
    coord_table="${gff_file_sorted%.gff3}.$n"
    bedtools_pos="${pos_strand%.gff3}.fasta"
    bedtools_neg="${neg_strand%.gff3}.fasta"
    python $scripts/get_promoter_gff3.py $gff_file_sorted $n
#For some reason, the start position of the extracted sequence in the fasta
#header output by getfasta is incorrect, and -1 than real position.
#HOWEVER, the extracted sequence is as per the gff file.
    bedtools getfasta -fi $dna -bed $pos_strand -fo $bedtools_pos
    bedtools getfasta -fi $dna -bed $neg_strand -fo $bedtools_neg
    python $scripts/substitute_names.py $bedtools_pos $coord_table
    python $scripts/substitute_names.py $bedtools_neg $coord_table
#Important, bedtools extracts the sequence for negative strand, as it is, without
#reverse complement necessary for scan of the promoter region!
#After FASTA sequence extraction need to reverse complement sequences
#in the file containing minus strand sequences:
    sub_output_neg="${bedtools_neg%.fasta}_n.fasta"
    sub_output_pos="${bedtools_pos%.fasta}_n.fasta"
    revcom_output_neg="${bedtools_neg%.fasta}_n_revcom.fasta"
    /home/armita/prog/emboss/EMBOSS-4.0.0/bin/revseq $sub_output_neg \
    $revcom_output_neg
#Concatenate all the genes into one file and turn into single-line fasta
    file="${dna%.fa*}_promoters_$n.fasta"
    cat $revcom_output_neg >> $file
    cat $sub_output_pos >> $file
    awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file>temp && mv temp $file
    done
done

#For Fus2, use Andy's gff2fasta and extract 1, 2 kbp upstream:
input=/home/sobczm/popgen/clock/DNA_genomes
scripts=/home/sobczm/bin/popgen/clock/motif_discovery
cd $input
$scripts/gff2fasta.pl Fus2_canu_contigs_hardmasked.fa Fus2_final_genes_appended.gff3 outFus2
#And change values in the script to extract 1 and 2 kbp upstream.
