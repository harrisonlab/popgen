#!/bin/bash
input=/home/sobczm/popgen/renseq/input/transcriptomes
#Extract contigs containing putative R genes from both strands of analysis.
#Exclude sequences with domain matches in plus and minus strands as likely result of merging of different transcripts during the assembly process.
#Extract the target sequences and reverse complement sequences on the negative strand.
#When faced with many isoforms of the same gene, use all of them for bait design


#Convert input FASTA into single line per sequences
file=$input/KIM/GBRQ01_1_fsa_nt_combined_kim.fasta
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file>temp && mv temp $file
#Repeat for:
file=$input/RAJ/GBJZ01_1_fsa_nt_raj.fasta
file=$input/NZ/GBGJ01_1_fsa_nt_nz.fasta
file=$input/MARIA/Trinity.fasta
file=$input/BRIAN/Trinity_11Oct2016.fasta

#Separate the gene list into those on the positive and negative strands.
for a in $input/*/*gene_list.txt
do
awk -F $"\t" '$2=="+" {print $1}' $a >"${a%.*}_pos.txt"
awk -F $"\t" '$2=="-" {print $1}' $a >"${a%.*}_neg.txt"
done

cd $input/KIM
#genes on positive strand
bioawk -cfastx 'BEGIN{while((getline k <"kim_gene_list_pos.txt")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"$seq}' GBRQ01_1_fsa_nt_combined_kim.fasta >> kim_rgenes.fasta
#first, reverse complement genes on a negative strand, and then save with the rest.
bioawk -cfastx 'BEGIN{while((getline k <"kim_gene_list_neg.txt")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"revcomp($seq)}' GBRQ01_1_fsa_nt_combined_kim.fasta >> kim_rgenes.fasta

cd $input/RAJ
#genes on positive strand
bioawk -cfastx 'BEGIN{while((getline k <"raj_gene_list_pos.txt")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"$seq}' GBJZ01_1_fsa_nt_raj.fasta >> raj_rgenes.fasta
#first, reverse complement genes on a negative strand, and then save with the rest.
bioawk -cfastx 'BEGIN{while((getline k <"raj_gene_list_neg.txt")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"revcomp($seq)}' GBJZ01_1_fsa_nt_raj.fasta >> raj_rgenes.fasta

cd $input/NZ
#genes on positive strand
bioawk -cfastx 'BEGIN{while((getline k <"nz_gene_list_pos.txt")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"$seq}' GBGJ01_1_fsa_nt_nz.fasta >> nz_rgenes.fasta
#first, reverse complement genes on a negative strand, and then save with the rest.
bioawk -cfastx 'BEGIN{while((getline k <"nz_gene_list_neg.txt")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"revcomp($seq)}' GBGJ01_1_fsa_nt_nz.fasta >> nz_rgenes.fasta

###Create a file with all R gene sequences
cd $input
for a in $input/*/*_rgenes.fasta
do
cat $a >> all_rgenes.fasta
done


##Check if the genes contain retrotransposon domains, if so get rid of them
##!!! That, and repeat-masking may not be enough for good enough specificity.
# Enquire with MYbaits if supplied transcriptomes can be used to assess specificity
# Alternatively, blast the baits to the transcriptome one by one
## Six-frame translation of sequences

## InterProScan

##Check how many gene clusters are recovered to estimate the true number of recovered R genes across all the transcriptomes
