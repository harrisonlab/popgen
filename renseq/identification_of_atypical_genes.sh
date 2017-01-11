#!/bin/bash
input=/home/sobczm/popgen/renseq/input/transcriptomes/plant_rgenes/common_onion/domains
scripts=/home/sobczm/bin/popgen/renseq
#Looking for a subset of RLKs: SRLK (S-receptor like kinase genes)
#among PFAM results. Modeled after additional_domains.sh
############Xing2013
#B_lectin(PF01453, IPR001480)-SLG(S_locus_glycop: PF00954,
# IPR000858)-PAN_APPLE (PF00024, PF14295, PF08276, SM00473, IPR003609)
#-TM-Kinase domain
############Chen2006 and Catanzariti2015
#SLG missing
#Just any lectin? -> keyword search
#Pan or Apple
#TM-kinase domain

####First, identify lectin domain-containing genes.
cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cd $input/$name
grep "lectin" onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | awk -v \
OFS="\t" '$1=$1' >${name}_lectin.txt
done

###Secondly, identify SLG domain containing genes.
cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cd $input/$name
grep "PF00954" onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | awk -v \
OFS="\t" '$1=$1' >${name}_slg.txt
done

#### Thirdly, identify genes containing Pan/Apple domain
cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cd $input/$name
grep "PF00024\|PF14295\|PF08276\|SM00473" onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | awk -v \
OFS="\t" '$1=$1' >${name}_panapple.txt
done

####Fourthly, identify genes containing a kinase domain
cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cd $input/$name
grep -i "kinase" onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | awk -v \
OFS="\t" '$1=$1' >${name}_kinase.txt
done

####Identify the genes containing potential transposon domains: we don't want them.
####Transposon InterProScan domain IDS taken from ultimate_bait_design.sh
####and changed to Pfam IDs.
cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cd $input/$name
grep -i "PF00078\|PF13482\|PF10551\|PF03221\|PF00078\|PF14214\|PF05699\|transpos*\|integrase" onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | awk -v \
OFS="\t" '$1=$1' | sed 's/_frame.*//g' >${name}_transposon.txt
done


##Find overlap of sequence names between the four categories above,
#or without the SLG domain, as per definition from Chen2006 and Catanzariti2015.
##withSLG
cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cd $input/$name
cut -f1 ${name}_lectin.txt | sed 's/_frame.*//g' >temp_ref
for files in ${name}_slg.txt ${name}_panapple.txt ${name}_kinase.txt
do
    cut -f1 $files | sed 's/_frame.*//g' >temp
    grep -Fx -f temp temp_ref >final_list_withslg
    cat final_list_withslg >temp_ref
done
#Eliminate the gene_ids matching those in the list with potential transposon genes

done


##noSLG
cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cd $input/$name
cut -f1 ${name}_lectin.txt >temp_ref
for files in ${name}_panapple.txt ${name}_kinase.txt
do
    cut -f1 $files >temp
    grep -Fx -f temp temp_ref >final_list_noslg
    cat final_list_noslg >temp_ref
done
done



#####Lastly, run Phobius on selected sequences to idenitfy transembrane domains
#####and check for any transposon domains




####Bait design
