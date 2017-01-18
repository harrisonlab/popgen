#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
#Generate trees out of aligned protein sequences of all detected RLP, RLK and SRLK genes. 


##Create input files containing all RLP, RLK and SRLK sequences
#RLP
#tree for contigs from a selected assembly: cornell

#RLK
#tree for contigs from a selected assembly: cornell

#SRLK
cd /home/sobczm/popgen/renseq/input/transcriptomes/ultimate_baits/srlk
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for n in "${names[@]}"
do
bioawk -cfastx 'BEGIN{while((getline k <"$n_final_list_withslg_no_transposon_sorted")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"$seq}' /home/sobczm/popgen/renseq/input/transcriptomes/rgaugury/$n/onion_${n}_protein.fa >${n}_srlk_prot.fa
done 
#Concatenate all sequences into one file 
for a in *srlk_prot.fa; do cat $a >>all_srlk_prot.fa; done

##Mafft alignment
qsub $scripts/sub_mafft.sh all_srlk_prot.fa

##FastTree generation (check if accepts stop codons)
qsub $scripts/sub_fasttree.sh all_srlk_prot_aligned.fa