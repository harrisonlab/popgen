#!/bin/bash
input=/home/sobczm/popgen/rnaseq/genomite
scripts=/home/sobczm/bin/popgen/rnaseq
original_data=/data/seq_data/external/20171110_genomite_rnaseq_cnag

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
python $scripts/targets.py >Targets.txt