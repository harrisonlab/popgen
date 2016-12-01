#!/bin/bash
trans=/home/sobczm/popgen/renseq/input/transcriptomes
input=/home/sobczm/popgen/renseq/input/transcriptomes/ultimate_baits
scripts=/home/sobczm/bin/popgen/renseq

#Copy assemblies to one directory and rename
cd $input
mkdir assemblies
cp $trans/BRIAN/Trinity_11Oct2016.fasta $input/assemblies/brian_assembly.fasta
cp $trans/CORNELL/cornell_Trinity.fasta $input/assemblies/cornell_assembly.fasta
cp $trans/H6/h6_Trinity.fasta $input/assemblies/h6_assembly.fasta
cp $trans/HAN/han_Trinity.fasta $input/assemblies/han_assembly.fasta
cp $trans/KIM/GBRQ01_1_fsa_nt_combined_kim.fasta $input/assemblies/kim_assembly.fasta
cp $trans/LIU/liu_Trinity.fasta $input/assemblies/liu_assembly.fasta
cp $trans/MARIA/Trinity.fasta $input/assemblies/maria_assembly.fasta
cp $trans/NZ/GBGJ01_1_fsa_nt_nz.fasta $input/assemblies/nz_assembly.fasta
cp $trans/RAJ/GBJZ01_1_fsa_nt_raj.fasta $input/assemblies/raj_assembly.fasta
cp $trans/SP3B/sp3b_Trinity.fasta $input/assemblies/sp3b_assembly.fasta
cp $trans/SUN/sun_Trinity.fasta $input/assemblies/sun_assembly.fasta

names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
#Convert input FASTA into single line per sequences
for file in *.fasta
do
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file>temp && mv temp $file
done

#Get unique IDs only in target gene lists
for a in $input/lists/*/*.txt
do
cp $a temp
sort temp | uniq >$a
done

#Separate the gene list into those on the positive and negative strands.
for a in $input/lists/*/*.txt
do
awk -F $"\t" '$2=="+" {print $1}' $a >"${a%.*}_pos.txt"
awk -F $"\t" '$2=="-" {print $1}' $a >"${a%.*}_neg.txt"
done

#Extract NBS and TIR domains containing contigs in each assembly
#genes on positive strand
for n in "${names[@]}"
do
bioawk -cfastx 'BEGIN{while((getline k <"$input/lists/nbs/$n_nbs_pos.txt")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"$seq}' $input/assemblies/${n}_assembly.fasta >> $input/nbs/${n}_nbs.fasta
#first, reverse complement genes on a negative strand, and then save with the rest.
bioawk -cfastx 'BEGIN{while((getline k <"$input/lists/nbs/${n}_nbs_neg.txt")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"revcomp($seq)}' $input/assemblies/${n}_assembly.fasta >> $input/nbs/${n}_nbs.fasta
done
#Extract alternative domain (Avr cleaveage) containing contigs in each assemnly
for n in "${names[@]}"
do
bioawk -cfastx 'BEGIN{while((getline k <"$input/lists/avr/$n_avr_pos.txt")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"$seq}' $input/assemblies/${n}_assembly.fasta >> $input/avr/${n}_avr.fasta
#first, reverse complement genes on a negative strand, and then save with the rest.
bioawk -cfastx 'BEGIN{while((getline k <"$input/lists/avr/${n}_avr_neg.txt")>0)i[k]=1}{if(i[$name])print \
">"$name"\n"revcomp($seq)}' $input/assemblies/${n}_assembly.fasta >> $input/avr/${n}_avr.fasta
done

## Six-frame translation of sequences
for a in $input/avr/*.fasta
do
b=$( echo $a | sed -e 's/.fasta/\_prot.fasta/' )
java -jar $scripts/Translate6Frame.jar -i $a -o $b
done

for a in $input/nbs/*.fasta
do
b=$( echo $a | sed -e 's/.fasta/\_prot.fasta/' )
java -jar $scripts/Translate6Frame.jar -i $a -o $b
done

#Remove all the asterisks (stop codons) in the protein sequences, as InterProScan
#does not accept sequences containing stop codons.
for a in $input/avr/*prot.fasta
do
b=$( echo $a | sed -e 's/.fasta/\_nostop.fasta/' )
sed 's/\*//g' $a > $b
qsub $scripts/sub_interproscan.sh $b
done

for a in $input/nbs/*prot.fasta
do
b=$( echo $a | sed -e 's/.fasta/\_nostop.fasta/' )
sed 's/\*//g' $a > $b
qsub $scripts/sub_interproscan.sh $b
done

#Check the six datasets for the presence of transposon domains. Keywords:
#IPR000477 # Reverse transcriptase (RT) catalytic domain profile
#IPR012337 # Ribonuclease H domain
#IPR018289 # MULE transposase domain
#IPR006600 # Tc5 transposase DNA-binding domain
#IPR000477 # Reverse transcriptase (RNA-dependent)
#IPR025476 # Helitron helicase-like domain
#IPR008906 # hAT family C-terminal dimerisation region
#transpos*
#integrase
#One group of contigs containing integrase domains was discovered. Filter the input FASTA
#file to remove them.
cd $input/avr
for name in "${names[@]}"
do
cat ${name}_avr_prot_nostop.fasta.tsv | sort >temp_file
grep 'IPR000477\|IPR004875\|IPR025476\|IPR012337\|IPR018289\|IPR006600\|IPR000477\|IPR008906\|transpos*\|integrase' temp_file | sed -e 's/_frame.*//' | sed -e 's/_[[:digit:]].*//' | cut -d" " -f1 >${name}_no_transposon.txt
python $scripts/keep_out_genes.py ${name}_no_transposon.txt ${name}_avr.fasta
done

cd $input/nbs
for name in "${names[@]}"
do
cat ${name}_nbs_prot_nostop.fasta.tsv | sort >temp_file
grep 'IPR000477\|IPR004875\|IPR025476\|IPR012337\|IPR018289\|IPR006600\|IPR000477\|IPR008906\|transpos*\|integrase' temp_file | sed -e 's/_frame.*//' | sed -e 's/_[[:digit:]].*//' | cut -d" " -f1 >${name}_no_transposon.txt
python $scripts/keep_out_genes.py ${name}_no_transposon.txt ${name}_nbs.fasta
done

#Found 3 contig groups containing potential transposons and removed those from further analysis

###Create a file with all R gene sequences
names_co=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" )
names_wo=( "liu" "sun" )

cd $input/avr

for name in "${names_co[@]}"
do
cat ${name}_avr_filtered.fasta >> ../all_rgenes.fasta
cat ${name}_avr_filtered.fasta >> ../common_onion_rgenes.fasta
cat ${name}_avr_filtered.fasta >> ../all_avr_genes.fasta
cat ${name}_avr_filtered.fasta >> ../common_onion_avr_genes.fasta
done

for name in "${names_wo[@]}"
do
cat ${name}_avr_filtered.fasta >> ../all_rgenes.fasta
cat ${name}_avr_filtered.fasta >> ../welsh_onion_rgenes.fasta
cat ${name}_avr_filtered.fasta >> ../all_avr_genes.fasta
cat ${name}_avr_filtered.fasta >> ../welsh_onion_avr_genes.fasta
done

cd $input/nbs

for name in "${names_co[@]}"
do
cat ${name}_nbs_filtered.fasta >> ../all_rgenes.fasta
cat ${name}_nbs_filtered.fasta >> ../common_onion_rgenes.fasta
cat ${name}_nbs_filtered.fasta >> ../all_nbs_genes.fasta
cat ${name}_nbs_filtered.fasta >> ../common_onion_nbs_genes.fasta
done

for name in "${names_wo[@]}"
do
cat ${name}_nbs_filtered.fasta >> ../all_rgenes.fasta
cat ${name}_nbs_filtered.fasta >> ../welsh_onion_rgenes.fasta
cat ${name}_nbs_filtered.fasta >> ../all_nbs_genes.fasta
cat ${name}_nbs_filtered.fasta >> ../welsh_onion_nbs_genes.fasta
done

##Check how many gene clusters are recovered to estimate the true number of recovered R genes across all the transcriptomes
usearch=/home/sobczm/bin/usearch/usearch9.0.2132_i86linux32
fasta=all_rgenes.fasta
#79 clusters
$usearch -cluster_fast $fasta -id 0.5 -sort length -uc all_rgenes_0.5.clusters
#67 clusters
$usearch -cluster_fast $fasta -id 0.4 -sort length -uc all_rgenes_0.5.clusters
#57 clusters
$usearch -cluster_fast $fasta -id 0.3 -sort length -uc all_rgenes_0.5.clusters

fasta=all_nbs_genes.fasta
#48 clusters
$usearch -cluster_fast $fasta -id 0.5 -sort length -uc all_nbs_0.5.clusters
#41 clusters
$usearch -cluster_fast $fasta -id 0.4 -sort length -uc all_nbs_0.5.clusters
#37 clusters
$usearch -cluster_fast $fasta -id 0.3 -sort length -uc all_nbs_0.5.clusters

#bait Design
n=4
python $scripts/create_baits.py --inp $input/all_rgenes.fasta --coverage $n \
--out all_rgenes_baits_"$n"x.fasta
#Results in around 35,000 baits

#Check for presence of Ns and repetitive sequences
fasta=all_rgenes_baits_"$n"x.fasta
#Hard-mask repetitive regions with Ns
$usearch -fastx_mask $fasta -qmask dust -fastaout "${fasta%.*}"_masked.fasta -hardmask
#Remove the baits containing Ns.
python $scripts/remove_N_fasta.py all_rgenes_baits_"$n"x_masked.fasta
#Results in retention of around 33835 baits, ie. around 4% baits discarded
#Cluster the baits at 95% identity
fasta=all_rgenes_baits_"$n"x_masked_noN.fasta
id_threshold=0.95
$usearch -cluster_fast $fasta -id $id_threshold -sort length \
-consout ${fasta%.*}_"$id_threshold"_clust.fasta -uc ${fasta%.*}_"$id_threshold".clusters
#Contains 20259 unique sequences subsequently clustered into 1458 clusters
id_threshold=0.99
#1829 clusters
#to fix the ends for no end gaps, use -leftjust and -rightjust
#99%
$usearch -cluster_fast $fasta -id $id_threshold -sort length -leftjust -rightjust \
-consout ${fasta%.*}_"$id_threshold"_nogap_clust.fasta -uc ${fasta%.*}_"$id_threshold"_nogap.clusters
#Results in 19981 bait sequences
#95%
$usearch -cluster_fast $fasta -id $id_threshold -sort length -leftjust -rightjust \
-consout ${fasta%.*}_"$id_threshold"_nogap_clust.fasta -uc ${fasta%.*}_"$id_threshold"_nogap.clusters
#Results in 19589 bait sequences, and at 90% in 19459 sequences

#Suggestion for the future: Use 5x coverage for initial bait design
n=5
python $scripts/create_baits.py --inp $input/all_rgenes.fasta --coverage $n \
--out all_rgenes_baits_"$n"x.fasta
#Check for presence of Ns and repetitive sequences
fasta=all_rgenes_baits_"$n"x.fasta
#Hard-mask repetitive regions with Ns
$usearch -fastx_mask $fasta -qmask dust -fastaout "${fasta%.*}"_masked.fasta -hardmask
#Remove the baits containing Ns.
python $scripts/remove_N_fasta.py all_rgenes_baits_"$n"x_masked.fasta

fasta=all_rgenes_baits_"$n"x_masked_noN.fasta
id_threshold=0.90
$usearch -cluster_fast $fasta -id $id_threshold -sort length -leftjust -rightjust \
-consout ${fasta%.*}_"$id_threshold"_nogap_clust.fasta -uc ${fasta%.*}_"$id_threshold"_nogap.clusters
#24701 unique baits clustered into clusters at 90% identity
