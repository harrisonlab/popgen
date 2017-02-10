#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/apple

#Parse the output of RGAugury runs.
names=( "Li" "Velasco" )

###Combine the results from small split files
for name in "${names[@]}"
do
cd $input/rgaugury/$name
for f in ./final/*.TMCC.candidates.lst
do
    cat $f >> ${name}.TMCC.candidates.lst
done
for f in ./final/*.RLP.candidates.lst
do
    cat $f >> ${name}.RLP.candidates.lst
done
for f in ./final/*.RLK.candidates.lst
do
    cat $f >> ${name}.RLK.candidates.lst
done
for f in ./final/*.RGA.candidates.lst
do
    cat $f >> ${name}.RGA.candidates.lst
done
for f in ./final/*.NBS.candidates.lst
do
    cat $f >> ${name}.NBS.candidates.lst
done
done

#################SParse the NBS-domain containing genes into
#################Slists with different subtypes and move into a separate folder.

for name in "${names[@]}"
do
cd $input/rgaugury/$name
mkdir lists
file=${name}.NBS.candidates.lst
cut -f1 $file > ./lists/${name}.NBS.all.lst
awk '$2 == "TN" {print $1}' $file > ./lists/${name}.NBS.TN.lst
awk '$2 == "CN" {print $1}' $file > ./lists/${name}.NBS.CN.lst
awk '$2 == "NL" {print $1}' $file > ./lists/${name}.NBS.NL.lst
awk '$2 == "TNL" {print $1}' $file > ./lists/${name}.NBS.TNL.lst
awk '$2 == "CNL" {print $1}' $file > ./lists/${name}.NBS.CNL.lst
awk '$2 == "NBS" {print $1}' $file > ./lists/${name}.NBS.NBS.lst
awk '$2 == "TX" {print $1}' $file > ./lists/${name}.NBS.TX.lst
awk '$2 == "OTHER" {print $1}' $file > ./lists/${name}.NBS.OTHER.lst
done


#################SParse the RLK genes into lists with different subtypes and move into a separate folder.
for name in "${names[@]}"
do
cd $input/rgaugury/$name
file=${name}.RLK.candidates.lst
cut -f1 $file > ./lists/${name}.RLK.all.lst
awk '$3 == "lrr" {print $1}' $file > ./lists/${name}.RLK.LRR.lst
awk '$3 == "lysm" {print $1}' $file > ./lists/${name}.RLK.LysM.lst
awk '$3 == "other_receptor" {print $1}' $file > ./lists/${name}.RLK.other.lst
done

#################SParse the RLP genes into lists with different subtypes and move into a separate folder.
for name in "${names[@]}"
do
cd $input/rgaugury/$name
file=${name}.RLP.candidates.lst
cut -f1 $file > ./lists/${name}.RLP.all.lst
awk '$3 == "lrr" {print $1}' $file > ./lists/${name}.RLP.LRR.lst
awk '$3 == "lysm" {print $1}' $file > ./lists/${name}.RLP.LysM.lst
done

#################STM-CC genes
for name in "${names[@]}"
do
cd $input/rgaugury/$name
file=${name}.TMCC.candidates.lst
cut -f1 $file > ./lists/${name}.TMCC.all.lst
done

#################SRLK genes
for name in "${names[@]}"
do
cd $input/rgaugury/$name/final
for file in myseq*.pfam.local.search.out
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
for file2 in myseq*.phobius.txt
do
threshold=1
awk -v threshold=1 '$2 > threshold' $file2 | cut -d ' ' -f1 | awk -F"_frame" '{print $1}' >>${name}_tmm.txt
done
done

##Find overlap of sequence names between the five categories above
for name in "${names[@]}"
do
cd $input/rgaugury/$name/final
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
grep -f temp.diff ${name}_kinase.txt >$input/rgaugury/$name/lists/${name}_slrk.lst
done


#In all lists, get rid of genes with potential transposon domains.
####Transposon InterProScan domain IDS taken from ultimate_bait_design.sh and changed to Pfam IDs.
######
#To do that, get a list of all proteins containing
#potential transposon domains
for name in "${names[@]}"
do
cd $input/rgaugury/$name/final
for file in myseq*.pfam.local.search.out
do
grep -i "PF00078\|PF13482\|PF10551\|PF03221\|PF00078\|PF14214\|PF05699\|transpos*\|integr*s" $file | cut -d ' ' -f1 | awk -F"_frame" '{print $1}' >> ../lists/${name}.transposon.all
sort ../lists/${name}.transposon.all | uniq >temp.list
mv temp.list ../lists/${name}.transposon.all
done
done

#Eliminate the gene_ids matching those in the list with potential transposon genes from all lists
for name in "${names[@]}"
do
cd $input/rgaugury/$name/lists
for a in *.lst
do
sort $a | uniq >temp.list
comm -2 temp.list ${name}.transposon.all >temp.diff
mv temp.diff ${a%.lst}.no_transposon
done
done