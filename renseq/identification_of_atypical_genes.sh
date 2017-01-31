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

cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cd $input/$name
####First, identify lectin domain-containing genes.
grep "lectin" onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | awk -v \
OFS="\t" '$1=$1' >${name}_lectin.txt
###Secondly, identify SLG domain containing genes.
grep "PF00954" onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | awk -v \
OFS="\t" '$1=$1' >${name}_slg.txt
#### Thirdly, identify genes containing Pan/Apple domain
grep "PF00024\|PF14295\|PF08276\|SM00473" onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | awk -v \
OFS="\t" '$1=$1' >${name}_panapple.txt
####Fourthly, identify genes containing a kinase domain
grep -i "kinase" onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | awk -v \
OFS="\t" '$1=$1' >${name}_kinase.txt
####Identify the genes containing potential transposon domains: we don't want them.
####Transposon InterProScan domain IDS taken from ultimate_bait_design.sh
####and changed to Pfam IDs.
grep -i "PF00078\|PF13482\|PF10551\|PF03221\|PF00078\|PF14214\|PF05699\|transpos*\|integrase" onion_167_TAIR10.protein.fa_pfamscan-04-08-2014.out | awk -v \
OFS="\t" '$1=$1' | sort >${name}_transposon.txt
done

##Find overlap of sequence names between the four categories above,
#or without the SLG domain, as per definition from Chen2006 and Catanzariti2015.
##withSLG
cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cd $input/$name
cut -f1 ${name}_lectin.txt | sort >temp_ref
for files in ${name}_slg.txt ${name}_panapple.txt ${name}_kinase.txt
do
    cut -f1 $files | sort >temp
    grep -Fx -f temp temp_ref >final_list_withslg
    cat final_list_withslg >temp_ref
done
#Eliminate the gene_ids matching those in the list with potential transposon genes
comm -2 final_list_withslg ${name}_transposon.txt | sort | uniq >${name}_final_list_withslg_no_transposon
done

##noSLG
cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cd $input/$name
cut -f1 ${name}_lectin.txt | sort  >temp_ref
for files in ${name}_panapple.txt ${name}_kinase.txt
do
    cut -f1 $files | sort >temp
    grep -Fx -f temp temp_ref >final_list_noslg
    cat final_list_noslg >temp_ref
done
#Eliminate the gene_ids matching those in the list with potential transposon genes
comm -2 final_list_noslg ${name}_transposon.txt | sort | uniq >${name}_final_list_noslg_no_transposon
done

##Copy the selected lists to a directory for bait design.
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cd $input/$name
cp -r ${name}_final_list_noslg_no_transposon ${name}_final_list_withslg_no_transposon \
/home/sobczm/popgen/renseq/input/transcriptomes/ultimate_baits/srlk
done

input2=/home/sobczm/popgen/renseq/input/transcriptomes/ultimate_baits/srlk
trans=/home/sobczm/popgen/renseq/input/transcriptomes
#Copy assemblies to one directory and rename
cp $trans/BRIAN/Trinity_11Oct2016.fasta $input2/brian_assembly.fasta
cp $trans/CORNELL/cornell_Trinity.fasta $input2/cornell_assembly.fasta
cp $trans/H6/h6_Trinity.fasta $input2/h6_assembly.fasta
cp $trans/HAN/han_Trinity.fasta $input2/han_assembly.fasta
cp $trans/KIM/GBRQ01_1_fsa_nt_combined_kim.fasta $input2/kim_assembly.fasta
cp $trans/LIU/liu_Trinity.fasta $input2/liu_assembly.fasta
cp $trans/MARIA/Trinity.fasta $input2/maria_assembly.fasta
cp $trans/NZ/GBGJ01_1_fsa_nt_nz.fasta $input2/nz_assembly.fasta
cp $trans/RAJ/GBJZ01_1_fsa_nt_raj.fasta $input2/raj_assembly.fasta
cp $trans/SP3B/sp3b_Trinity.fasta $input2/sp3b_assembly.fasta
cp $trans/SUN/sun_Trinity.fasta $input2/sun_assembly.fasta

#Prepare lists with genes on positive and negative strands
cd $input2
for a in *no_transposon
do
python $scripts/add_strand_column2.py $a
done

#Get unique IDs only in target gene lists
for a in *.strand
do
cp $a temp
sort temp | uniq >$a
done

#Separate the gene list into those on the positive and negative strands.
for a in *.strand
do
awk -F $"\t" '$2=="+" {print $1}' $a >"${a%.*}_pos.txt"
awk -F $"\t" '$2=="-" {print $1}' $a >"${a%.*}_neg.txt"
done

#Extract select contigs in each assembly
#genes on positive strand
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for n in "${names[@]}"
do
bioawk -cfastx 'BEGIN{while((getline k <"$n_final_list_noslg_no_transpos_pos.txt")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"$seq}' ${n}_assembly.fasta >> ${n}_noslg.fasta
#first, reverse complement genes on a negative strand, and then save with the rest.
bioawk -cfastx 'BEGIN{while((getline k <"$n_final_list_noslg_no_transpos_neg.txt")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"revcomp($seq)}' ${n}_assembly.fasta >> ${n}_noslg.fasta
done

####Bait design
#Merge the sequences into one file:
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for name in "${names[@]}"
do
cat ${name}_noslg.fasta >> all_noslg.fasta
cat ${name}_withslg.fasta >> all_withslg.fasta
done

##Cluster how many gene clusters are recovered to estimate the true number of recovered SLRK genes across all the transcriptomes

usearch=/home/sobczm/bin/usearch/usearch9.0.2132_i86linux32
#386 unique sequences
fasta=all_noslg.fasta
#34 clusters
$usearch -cluster_fast $fasta -id 0.5 -sort length -uc all_noslg_0.5.clusters
#25 clusters
$usearch -cluster_fast $fasta -id 0.4 -sort length -uc all_noslg_0.4.clusters
#23 clusters
$usearch -cluster_fast $fasta -id 0.3 -sort length -uc all_noslg_0.3.clusters

usearch=/home/sobczm/bin/usearch/usearch9.0.2132_i86linux32
#343 unique sequences
fasta=all_withslg.fasta
#33 clusters
$usearch -cluster_fast $fasta -id 0.5 -sort length -uc all_withslg_0.5.clusters
#25 clusters
$usearch -cluster_fast $fasta -id 0.4 -sort length -uc all_withslg_0.4.clusters
#18 clusters
$usearch -cluster_fast $fasta -id 0.3 -sort length -uc all_withslg_0.3.clusters

#bait design
n=5
python $scripts/create_baits.py --inp all_withslg.fasta --coverage $n \
--out all_withslg_genes_baits_"$n"x.fasta

#Check for presence of Ns and repetitive sequences
fasta=all_withslg_genes_baits_"$n"x.fasta
#Hard-mask repetitive regions with Ns
$usearch -fastx_mask $fasta -qmask dust -fastaout "${fasta%.*}"_masked.fasta -hardmask
#Remove the baits containing Ns.
python $scripts/remove_N_fasta.py all_withslg_genes_baits_"$n"x_masked.fasta

fasta=all_withslg_genes_baits_"$n"x_masked_noN.fasta
id_threshold=0.90
$usearch -cluster_fast $fasta -id $id_threshold -sort length -leftjust -rightjust \
-consout ${fasta%.*}_"$id_threshold"_nogap_clust.fasta -uc ${fasta%.*}_"$id_threshold"_nogap.clusters
#38695 unique baits clustered into 22226 clusters at 90% identity.
