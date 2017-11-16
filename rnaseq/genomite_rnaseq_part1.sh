#!/bin/bash
input=/home/sobczm/popgen/rnaseq/genomite
scripts=/home/sobczm/bin/popgen/rnaseq
original_data=/data/seq_data/external/20171110_genomite_rnaseq_cnag

#Mite genome and gff annotation obtained from V. Zhurov.
#Vesca ver 2.0 genome and annotation ver. v2.0.a2 used for the plant side.

#Raw Data QC
for RawData in $original_data/GENOMITE_03/20171107/FASTQ/*.fastq.gz $original_data/GENOMITE_04/20171107/FASTQ/*.fastq.gz
do
Jobs=$(qstat | grep 'sub_run_fa' | wc -l)
while [ $Jobs -gt 10 ]
do
    sleep 10
    printf "."
    Jobs=$(qstat | grep 'sub_run_fa' | wc -l)
done
qsub $scripts/sub_run_fastqc.sh $RawData
done

#Read trimming on quality and adapter
cd $input/Genomite3
for reads in $original_data/GENOMITE_03/20171107/FASTQ/*1.fastq.gz
do
    Jobs=$(qstat | grep 'sub_read_q' | wc -l)
    while [ $Jobs -gt 10 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_read_q' | wc -l)
    done
reads2=$(echo "$reads" | sed 's/1.fastq.gz/2.fastq.gz/')
qsub $scripts/sub_read_qc.sh $reads $reads2
done

cd $input/Genomite4
for reads in $original_data/GENOMITE_04/20171107/FASTQ/*1.fastq.gz
do
    Jobs=$(qstat | grep 'sub_read_q' | wc -l)
    while [ $Jobs -gt 10 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_read_q' | wc -l)
    done
reads2=$(echo "$reads" | sed 's/1.fastq.gz/2.fastq.gz/')
qsub $scripts/sub_read_qc.sh $reads $reads2
done


#Data QC after trimming
cd $input/Genomite3
for RawData in $input/Genomite3/*trim.fq.gz
do
Jobs=$(qstat | grep 'sub_run_fa' | wc -l)
while [ $Jobs -gt 10 ]
do
    sleep 10
    printf "."
    Jobs=$(qstat | grep 'sub_run_fa' | wc -l)
done
qsub $scripts/sub_run_fastqc.sh $RawData
done

cd $input/Genomite4
for RawData in $input/Genomite4/*trim.fq.gz
do
Jobs=$(qstat | grep 'sub_run_fa' | wc -l)
while [ $Jobs -gt 10 ]
do
    sleep 10
    printf "."
    Jobs=$(qstat | grep 'sub_run_fa' | wc -l)
done
qsub $scripts/sub_run_fastqc.sh $RawData
done


#Create tables matching up sample library names, sample conditions to filenames with reads.
cd $input/genomite_samples
Rscript --vanilla genomite_samples.R

#Create Star genome indices.
cd $input/assemblies/mite
STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles tetur_200909.fa --sjdbGTFfile tetur_current.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 --genomeSAindexNbases 12 --runThreadN 10
cd $input/assemblies/strawberry
STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles Fragaria_vesca_v2.0.a1_pseudomolecules.fasta --sjdbGTFfile f.vesca2.0.a2.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 --genomeSAindexNbases 12 --runThreadN 10

#Map reads to the diploid strawberry genome
cd $input/Genomite3 
indexed_assembly=$input/assemblies/strawberry
samples=$input/genomite_samples/genomite3_sample
while read line
do
    if [ $count == 0 ]
    then
        ((count += 1)) 
    else
        set -- $line
        out_dir=''$input'/Genomite3/'$3'_'$4'_'$5'_'$6'h_'$7'/'
        qsub $scripts/sub_star_genomite.sh $indexed_assembly $PWD/$1 $PWD/$2 $out_dir
    fi
done < $samples


#Map reads to the mite genome
cd $input/Genomite4 
indexed_assembly=$input/assemblies/mite
samples=$input/genomite_samples/genomite4_sample
while read line
do
    if [ $count == 0 ]
    then
        ((count += 1)) 
    else
        set -- $line
        out_dir=''$input'/Genomite4/'$3'_'$4'_'$5'_'$6'h_'$7'/'
        echo $out_dir
        qsub $scripts/sub_star_genomite.sh $indexed_assembly $PWD/$1 $PWD/$2 $out_dir
    fi
done < $samples