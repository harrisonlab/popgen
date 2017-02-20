#!/bin/bash
input=/home/sobczm/popgen/summary_stats
scripts=/home/sobczm/bin/popgen/summary_stats

input_hap=/home/sobczm/popgen/other/passey/Ash_farm_172_pacbio_contigs_unmasked_no118.vcf
input_dip=/home/groups/harrisonlab/project_files/phytophthora_fragariae/summary_stats/95m_contigs_unmasked_UK123.vcf

#Create a directory for each individual BayeScan runS
mkdir -p $input/outliers/prfrag
mkdir -p $input/outliers/vinequalis

#Copy input files for the analysis
cp -r $input_hap $input/outliers/vinequalis
cp -r $input_dip $input/outliers/prfrag
cd $input/outliers

####################Filter SNPs to retain only biallelic SNPs, otherwise not compatible with PGDSpider
#and Bayescan. Furthermore, keep only SNPs with max 5% missing genotypes.
#!! It may be necessary to also pre-filter for samples which were poorly sequenced/
#aligned beforehand to avoid removal of too many potentially informative SNPs.
for filename in $input/outliers/*/*.vcf
do
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf $filename --max-missing 0.95 --min-alleles 2 --max-alleles 2 --recode --out ${filename%.vcf}_bi_filtered
done


#############PGDSpider conversion from VCF to Bayescan format


#First need to prepare a simple file with custom assignment of each individual sample
#name to one of the populations identified using PCA, STRUCTURE analysis etc. 
#Format: Sample_name Population_name
#But first: double check what sample names are present in the following VCF file:
grep "#CHROM" $filename | head -1
filename=$input/outliers/prfrag/95m_contigs_unmasked_UK123_bi_filtered.recode.vcf

#Resulting example file population assignment file: 
cat $input/outliers/prfrag/95m_contigs_unmasked_UK123_bi_filtered.recode.lst

#Now, need to prepare the configuration file. Copy the conversion script for diploid
#organisms
pgdspid=/home/sobczm/bin/PGDSpider_2.1.0.3
config=vcf_to_bayescan_diploid.spid
cp $pgdspid/$config ./
#Now, need to edit the conversion script to include the population assignment file
pop_assignment_file=$input/outliers/pfrag/95m_contigs_unmasked_UK123_bi_filtered.recode.lst
sed -i 's,^\(VCF_PARSER_POP_FILE_QUESTION=\).*,\1'"$pop_assignment_file"',' $config

java -jar $pgdspid/PGDSpider2-cli.jar -inputfile $filename -inputformat VCF -outputfile ${filename%.vcf}.geno -outputformat GESTE_BAYE_SCAN -spid $config
dos2unix ${filename%.vcf}.geno

#Bayescan run.
qsub $scripts/sub_bayescan.sh ${filename%.vcf}.geno

#Plot Bayescan results in R 
Rscript --vanilla $script/plot_bayescan.R 