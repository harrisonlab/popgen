#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/apple_2017

#Re-prediction of R genes and new bait design using the latest Golden Delicious genome from Daccord et al. (2017)

#Download genome, genome annotation and proteins.
wget https://iris.angers.inra.fr/gddh13/downloads/GDDH13_1-1_formatted.fasta.bz2
wget https://iris.angers.inra.fr/gddh13/downloads/gene_models_20170612.gff3.bz2
wget https://iris.angers.inra.fr/gddh13/downloads/GDDH13_1-1_prot.fasta.bz2
wget https://iris.angers.inra.fr/gddh13/downloads/GDDH13_1-1_TE.gff3.bz2

#Extract CDS from mRNA.
awk '$3 == "CDS" {print $0}' gene_models_20170612.gff3 >gene_models_20170612_cds.gff3
##Extract the sequences of the genes to be used as blast database (CDS + 5'/3' UTR)
genome_seq=GDDH13_1-1_formatted.fasta
gffread=/home/sobczm/bin/gffread/gffread/gffread
$gffread -w GDDH13_1-1_formatted_cds.fasta -g $genome_seq gene_models_20170612_cds.gff3
sed -i 's/mRNA://g' GDDH13_1-1_formatted_cds.fasta 
#Split the protein sequences into small files with 100 protein sequences and run 50 at a time.
file=GDDH13_1-1_prot.fasta
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file>temp && mv temp $file
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%100==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $file
for file in myseq*.fa
do
    Jobs=$(qstat | grep 'sub_rgaugu' | wc -l)
    while [ $Jobs -gt 100 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_rgaugu' | wc -l)
    done
qsub $scripts/sub_rgaugury.sh $file
done

#Run Phobius and final prediction
rm *RLKorRLP.domain.prediction.txt
rgaugury=/home/sobczm/bin/rgaugury/RGAugury.pl
for file in myseq*.fa
do
perl $rgaugury -p $file 
done

#Parse the output of RGAugury runs.
name=daccord
for f in ./rgaugury/*.TMCC.candidates.lst
do
    cat $f >> ${name}.TMCC.candidates.lst
done
for f in ./rgaugury/*.RLP.candidates.lst
do
    cat $f >> ${name}.RLP.candidates.lst
done
for f in ./rgaugury/*.RLK.candidates.lst
do
    cat $f >> ${name}.RLK.candidates.lst
done
for f in ./rgaugury/*.RGA.candidates.lst
do
    cat $f >> ${name}.RGA.candidates.lst
done
for f in ./rgaugury/*.NBS.candidates.lst
do
    cat $f >> ${name}.NBS.candidates.lst
done


#################SParse the NBS-domain containing genes into
#################Slists with different subtypes and move into a separate folder.

file=${name}.NBS.candidates.lst
cut -f1 $file > ${name}.NBS.all.lst
awk '$2 == "TN" {print $1}' $file > ${name}.NBS.TN.lst
awk '$2 == "CN" {print $1}' $file > ${name}.NBS.CN.lst
awk '$2 == "NL" {print $1}' $file > ${name}.NBS.NL.lst
awk '$2 == "TNL" {print $1}' $file > ${name}.NBS.TNL.lst
awk '$2 == "CNL" {print $1}' $file > ${name}.NBS.CNL.lst
awk '$2 == "NBS" {print $1}' $file > ${name}.NBS.NBS.lst
awk '$2 == "TX" {print $1}' $file > ${name}.NBS.TX.lst
awk '$2 == "OTHER" {print $1}' $file > ${name}.NBS.OTHER.lst

#################Parse the RLK genes into lists with different subtypes 
file=${name}.RLK.candidates.lst
cut -f1 $file > ${name}.RLK.all.lst
awk '$3 == "lrr" {print $1}' $file > ${name}.RLK.LRR.lst
awk '$3 == "lysm" {print $1}' $file > ${name}.RLK.LysM.lst
awk '$3 == "other_receptor" {print $1}' $file > ${name}.RLK.other.lst

#################Parse the RLP genes into lists with different subtypes and move into a separate folder.
file=${name}.RLP.candidates.lst
cut -f1 $file > ${name}.RLP.all.lst
awk '$3 == "lrr" {print $1}' $file > ${name}.RLP.LRR.lst
awk '$3 == "lysm" {print $1}' $file > ${name}.RLP.LysM.lst

#################TM-CC genes
file=${name}.TMCC.candidates.lst
cut -f1 $file > ${name}.TMCC.all.lst

#################SRLK genes
for file in ./rgaugury/myseq*.pfam.local.search.out
do
####First, identify lectin domain-containing genes.
grep "lectin" $file | cut -d ' ' -f1 | awk -F"_frame" '{print $1}'  >>${name}_lectin.txt
###Secondly, identify SLG domain containing genes.
grep "PF00954" $file | cut -d ' ' -f1 | awk -F"_frame" '{print $1}'  >>${name}_slg.txt
#### Thirdly, identify genes containing Pan/Apple domain
grep "PF00024\|PF14295\|PF08276\|SM00473"  $file | cut -d ' ' -f1 | awk -F"_frame" '{print $1}'  >>${name}_panapple.txt
####Fourthly, identify genes containing a kinase domain
grep -i "kinase" $file | cut -d ' ' -f1  >>${name}_kinase.txt
done
#And lastly identify the genes containing a TM domain
for file2 in ./rgaugury/myseq*.phobius.txt
do
threshold=1
awk -v threshold=1 '$2 > threshold' $file2 | cut -d ' ' -f1 | awk -F"_frame" '{print $1}' >>${name}_tmm.txt
done

##Find overlap of sequence names between the five categories above
cut -f1 ${name}_lectin.txt | sort >temp_ref
for files in ${name}_slg.txt ${name}_panapple.txt ${name}_tmm.txt
do
    cut -f1 $files | sort >temp_ls
    grep -Fx -f temp_ls temp_ref >final_list_withslg
    cat final_list_withslg >temp_ref
done
#find overlap with the kinase domain
awk -F"_frame" '{print $1}' ${name}_kinase.txt | sort >temp.kinase
comm -12 temp.kinase final_list_withslg >temp.diff
grep -f temp.diff ${name}_kinase.txt >${name}_slrk.lst

#In all lists, get rid of genes with potential transposon domains.
for file in ./rgaugury/myseq*.pfam.local.search.out
do
grep -i "PF00078\|PF13482\|PF10551\|PF03221\|PF00078\|PF14214\|PF05699\|transpos*\|integr*s" $file | cut -d ' ' -f1 | awk -F"_frame" '{print $1}' >> ${name}.transposon.all
sort ${name}.transposon.all | uniq >temp.list
mv temp.list ${name}.transposon.all
done

#Eliminate the gene_ids matching those in the list with potential transposon genes from all lists
for a in *.lst
do
sort $a | uniq >temp.list
comm -2 temp.list ${name}.transposon.all >temp.diff
mv temp.diff ${a%.lst}.no_transposon
done

#Extract select contigs in each assembly
for list in *.no_transposon
do
python $scripts/keep_list_genes2.py $list GDDH13_1-1_formatted_cds.fasta No
done
cp *.fasta $input/baits

#The plan is to have 1 bait library of NBS-LRR genes (20,000 baits) and one of RLK and RLP genes (20,000 baits).