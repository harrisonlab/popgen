#!/bin/bash
input=/home/sobczm/popgen/renseq/input/transcriptomes
scripts=/home/sobczm/bin/popgen/renseq
#Additional tool to annotate non-canonical R genes not containing the NBS domain.
#Using it to annotate the following groups of candidate genes:
#RLK (transmembrane-LRR/LysM-STTK domain)
#RLP (transmembrane-LRR/LysM domain)
#TMCC (transmembrane-coiled coil)

cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
for f in "${names[@]}"
do
fca=$(echo $f | tinput=/home/sobczm/popgen/renseq/input/transcriptomesr '[:lower:]' '[:upper:]')
#Create a seperate directory structure
mkdir -p $input/rgaugury/$f
cp $input/$fca/onion_${f}_protein.fa $input/rgaugury/$f
#Remove stop-codons from six-frame translated contigs
sed -i 's/\*//g' $input/rgaugury/$f/onion_${f}_protein.fa
#Break down the files into smaller chunks and run rgaugury for each transcriptome.
#As this taking too long, going to parallelise in a crude way by splitting
#into small files with 100 protein sequences and running 20 at a time when logged
#into a workernode.
cd $input/rgaugury/$f
file=onion_${f}_protein.fa
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
done
#Need to re-run steps involving Phobius - works only on the head node
cd $input
names=( "kim" "nz" "maria" "raj" "brian" "cornell" "h6" "han" "sp3b" "liu" "sun" )
rgaugury=/home/sobczm/bin/rgaugury/RGAugury.pl
for f in "${names[@]}"
do
cd $input/rgaugury/$f
rm *RLKorRLP.domain.prediction.txt
for file in myseq*.fa
do
perl $rgaugury -p $file
done
done

#Concatenate the output for RLK, RLP, NBS and TMCC genes.
names=( "cornell" "han" "liu" "sun" "h6" "sp3b" "brian" "maria" "raj" "kim" "nz" )
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
for f in ./final/*.RLKorRLP.merged.domains.txt
do
    lines=$(wc -l $f | cut -d" " -f1)
    if [ $lines > 1 ]; then
    sed -n 2,${lines}p $f >> ${name}.RLKorRLP.merged.domains.txt
    fi
done
for f in ./final/*.NBS.candidates.lst
do
    cat $f >> ${name}.NBS.candidates.lst
done
for f in ./final/*pfam*local*
do
    cat $f >> onion_${name}_protein.out
done
done

#Filter the results to obtain the list of putative RLK genes (TM-LRR-STTK or TM-LysM-STTK)
for name in "${names[@]}"
do
cd $input/rgaugury/$name
cp ${name}.RLKorRLP.merged.domains.txt ../RLK
done

cd $input/rgaugury/RLK
for f in *.txt
do
python $scripts/extract_rlk.py $f
done

for f in *.rlk
do
python $scripts/add_strand_column.py $f
done
