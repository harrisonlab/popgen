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
cd $input/assemblies/strawberry2
STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles Fragaria_vesca_v2.0.a1_pseudomolecules.fasta --sjdbGTFfile f.vesca2.0.a2.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 --genomeSAindexNbases 12 --runThreadN 10
cd $input/assemblies/strawberry1.1
STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles fvesca_v1.1_all.fa --sjdbGTFfile Fragaria_vesca_v1.1.a2.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 --genomeSAindexNbases 12 --runThreadN 10

#Map reads to the diploid strawberry genome ver 2.0
cd $input/Genomite3 
indexed_assembly=$input/assemblies/strawberry2
samples=$input/genomite_samples/genomite3_sample
while read line
do
    if [ $count == 0 ]
    then
        ((count += 1)) 
    else
        set -- $line
        out_dir=''$input'/Genomite3/strawberry2/'$3'_'$4'_'$5'_'$6'h_'$7'/'
        qsub $scripts/sub_star_genomite.sh $indexed_assembly $PWD/$1 $PWD/$2 $out_dir
    fi
done < $samples


#Map reads to the diploid strawberry genome ver 1.1
cd $input/Genomite3 
indexed_assembly=$input/assemblies/strawberry1.1
samples=$input/genomite_samples/genomite3_sample
while read line
do
    if [ $count == 0 ]
    then
        ((count += 1)) 
    else
        set -- $line
        out_dir=''$input'/Genomite3/strawberry1.1/'$3'_'$4'_'$5'_'$6'h_'$7'/'
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

#Collect mapping stats
#strawberry ver. 2
for File in $(ls ${input}/Genomite3/strawberry2/strawberry*/Log.final.out); do
Sample=$(echo $File | rev | cut -f2 -d '/' | rev | sed 's/vesca_//');
InputReads=$(cat $File | grep 'Number of input reads' | cut -f2);
ReadNumU=$(cat $File | grep 'Uniquely' | grep 'number' | cut -f2);
ReadPercU=$(cat $File | grep 'Uniquely' | grep '%' | cut -f2);
ReadNumM=$(cat $File | grep 'multiple' | grep 'Number' | cut -f2);
ReadPercM=$(cat $File | grep 'multiple' | grep '%' | cut -f2);
Mismatch=$(cat $File | grep 'Mismatch rate per base' | grep '%' | cut -f2);
echo -e "$Sample""\t""$InputReads""\t" "$ReadNumU""\t""$ReadPercU""\t""$ReadNumM""\t""$ReadPercM""\t""$Mismatch";  
done >strawberry_mapping_stats_ver2.txt

#strawberry ver. 1.1
for File in $(ls ${input}/Genomite3/strawberry1.1/strawberry*/Log.final.out); do
Sample=$(echo $File | rev | cut -f2 -d '/' | rev | sed 's/vesca_//');
InputReads=$(cat $File | grep 'Number of input reads' | cut -f2);
ReadNumU=$(cat $File | grep 'Uniquely' | grep 'number' | cut -f2);
ReadPercU=$(cat $File | grep 'Uniquely' | grep '%' | cut -f2);
ReadNumM=$(cat $File | grep 'multiple' | grep 'Number' | cut -f2);
ReadPercM=$(cat $File | grep 'multiple' | grep '%' | cut -f2);
Mismatch=$(cat $File | grep 'Mismatch rate per base' | grep '%' | cut -f2);
echo -e "$Sample""\t""$InputReads""\t" "$ReadNumU""\t""$ReadPercU""\t""$ReadNumM""\t""$ReadPercM""\t""$Mismatch";  
done >strawberry_mapping_stats_ver1.1.txt


#mite
for File in $(ls ${input}/Genomite4/mite*/Log.final.out); do
Sample=$(echo $File | rev | cut -f2 -d '/' | rev | sed 's/vesca_//');
InputReads=$(cat $File | grep 'Number of input reads' | cut -f2);
ReadNumU=$(cat $File | grep 'Uniquely' | grep 'number' | cut -f2);
ReadPercU=$(cat $File | grep 'Uniquely' | grep '%' | cut -f2);
ReadNumM=$(cat $File | grep 'multiple' | grep 'Number' | cut -f2);
ReadPercM=$(cat $File | grep 'multiple' | grep '%' | cut -f2);
Mismatch=$(cat $File | grep 'Mismatch rate per base' | grep '%' | cut -f2);
echo -e "$Sample""\t""$InputReads""\t" "$ReadNumU""\t""$ReadPercU""\t""$ReadNumM""\t""$ReadPercM""\t""$Mismatch";  
done >mite_mapping_stats.txt


#Count reads with HTSeq
cd $input
#strawberry ver. 2.0
strandedness="reverse"
gff=/home/sobczm/popgen/rnaseq/genomite/assemblies/strawberry2/f.vesca2.0.a2.gff3
for input_dir in ${input}/Genomite3/strawberry2/strawberry*
do
output=$(basename $input_dir)
qsub $scripts/sub_htseq.sh $strandedness ${input_dir}/Aligned.sortedByCoord.out.bam $gff ./htseq_out/${output}.out
done

#Only 10-20% mapped reads counted to annotated feature. Something wrong with ver 2 genome annotations. That is why going to use strawberry 1.1 genome below - on the plus side, it shows higher mapping rate as well. 

#First, need to fix formatting in the GFF file to be compatible with HTSeq.
sed 's/comment="overlaps bad segment.*//' $input/assemblies/strawberry1.1/Fragaria_vesca_v1.1.a2.gff3 > $input/assemblies/strawberry1.1/Fragaria_vesca_v1.1.a2.fixed.gff3

strandedness="reverse"
gff=/home/sobczm/popgen/rnaseq/genomite/assemblies/strawberry1.1/Fragaria_vesca_v1.1.a2.fixed.gff3
for input_dir in ${input}/Genomite3/strawberry1.1/strawberry*
do
output=$(basename $input_dir)
qsub $scripts/sub_htseq.sh $strandedness ${input_dir}/Aligned.sortedByCoord.out.bam $gff ./htseq_out/${output}_ver1.1.out
done

#mite
strandedness="reverse"
gff=/home/sobczm/popgen/rnaseq/genomite/assemblies/mite/tetur_current.gff3
for input_dir in ${input}/Genomite4/mite*
do
output=$(basename $input_dir)
qsub $scripts/sub_htseq.sh $strandedness ${input_dir}/Aligned.sortedByCoord.out.bam $gff ./htseq_out/${output}.out
done

#Sample subFeatureCount script run to obtain the gene lengths etc. used to format input for DESeq.
cd $input/htseq_out
gff=/home/sobczm/popgen/rnaseq/genomite/assemblies/strawberry1.1/Fragaria_vesca_v1.1.a2.fixed.gff3
qsub $scripts/sub_featureCounts.sh /home/sobczm/popgen/rnaseq/genomite/Genomite3/strawberry1.1/strawberry_nostress_nomite_3h_2/Aligned.sortedByCoord.out.bam $gff fc_strawberry_nostress_nomite_3h_2_ver1

gff=/home/sobczm/popgen/rnaseq/genomite/assemblies/mite/tetur_current.gff3
qsub $scripts/sub_featureCounts.sh /home/sobczm/popgen/rnaseq/genomite/Genomite4/mite_nostress_nonadapted_3h_5/Aligned.sortedByCoord.out.bam $gff fc_mite_nostress_nonadapted_3h_5

python $scripts/htseq2featurecount.py fc_strawberry_nostress_nomite_3h_2_ver1_featurecounts.txt strawberry_nostress_nomite_3h_2_ver1.1.out
#Parse HTSeq output to obtain FeatureCount-like output tables.
#strawberry ver1.1
for strawberry in *ver1.1.out
do
python $scripts/htseq2featurecount.py  fc_strawberry_nostress_nomite_3h_2_ver1_featurecounts.txt $strawberry
done
#mite
for mite in mite*.out
do
python $scripts/htseq2featurecount.py  fc_mite_nostress_nonadapted_3h_5_featurecounts.txt $mite
sed -i '1s/.out//' ${mite%.out}_fc.out
done

#QC of RNA-Seq results with DESeq2
#strawberry
cd $input/htseq_out/strawberry
Rscript --vanilla $scripts/DeSeq_Genomite3.R

#mite
cd $input/htseq_out/mite
Rscript --vanilla $scripts/DeSeq_Genomite4.R

#Collapse technical reps 
cd $input/genomite_samples
Rscript --vanilla join_sample_reps_ids.R

#FeatureCount-like output
cd $input/htseq_out/strawberry
python $scripts/merge_tech_reps.py $input/genomite_samples/strawberry_tech_replicates.txt _ver1.1_fc.out

cd $input/htseq_out/mite
python $scripts/merge_tech_reps.py $input/genomite_samples/mite_tech_replicates.txt _fc.out

#Use those merged samples to carry out basic DeSeq QC.
cd $input/htseq_out/strawberry_merged
python $scripts/print_merged_samples_table.py $input/genomite_samples/strawberry_tech_replicates.txt >strawberry_deseq_samples.txt
Rscript --vanilla $scripts/DeSeq_Genomite3.R

cd $input/htseq_out/mite_merged
python $scripts/print_merged_samples_table.py $input/genomite_samples/mite_tech_replicates.txt >mite_deseq_samples.txt
Rscript --vanilla $scripts/DeSeq_Genomite4.R

##Get merged HTSeq output
cd $input/htseq_out/strawberry_htseq_merged
python $scripts/merge_tech_reps_htseq.py $input/genomite_samples/strawberry_tech_replicates.txt _ver1.1.out

#Remove lines starting with "__" and move the merged files to a separate folder
for a in *AD*.out
do
sed -i '/^__/ d' $a
done

mkdir merged
mv *AD*.out ./merged

cd $input/htseq_out/mite_htseq_merged
python $scripts/merge_tech_reps_htseq.py  $input/genomite_samples/mite_tech_replicates.txt .out

#Remove lines starting with "__" and "mrna:" preceding gene name
for a in *AD*.out 
do
sed -i '/^__/ d' $a
sed -i 's/mRNA://' $a
done

mkdir merged
mv *AD*.out ./merged

##Generate Targets table as per V. Zhurov scripts - for DEG in R.
cd $input/htseq_out/mite_htseq_merged/merged
python $scripts/targets.py >Targets.txt

cd $input/htseq_out/strawberry_htseq_merged/merged
<<<<<<< HEAD
python $scripts/targets.py >vesca11_targets.txt

#Run Vlad's script for DEG and visualisation - all original samples
cd $input/htseq_out/mite_htseq_merged/merged
Rscript --vanilla $scripts/1_mite_data_normalisation.R 
Rscript --vanilla $scripts/2_mite_dge.R
Rscript --vanilla $scripts/3_mite_pca.R 

cd $input/htseq_out/strawberry_htseq_merged/merged
Rscript --vanilla $scripts/1_strawberry_data_normalisation.R 
Rscript --vanilla $scripts/2_strawberry_dge.R
Rscript --vanilla $scripts/3_strawberry_pca.R 
#Run Vlad's script for DEG and visualisation - selection of samples
=======
python $scripts/targets.py >Targets.txt
>>>>>>> c90f1c53e4146d309bd97c9769f8155fd8313044
