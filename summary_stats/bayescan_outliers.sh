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

#################### 1) Filter SNPs to retain only biallelic SNPs, otherwise not compatible with PGDSpider
#and Bayescan. Furthermore, keep only SNPs with max 5% missing genotypes.
### It may be necessary to also pre-filter for samples which were poorly sequenced/
#aligned beforehand to avoid removal of too many potentially informative SNPs.
vcftools=/home/sobczm/bin/vcftools/bin
for filename in $input/outliers/*/*.vcf
do
$vcftools/vcftools --vcf $filename --max-missing 0.95 --mac 1 --min-alleles 2 --max-alleles 2 --recode --out ${filename%.vcf}_bi_filtered
done


############# 2) PGDSpider conversion from VCF to Bayescan format
#First need to prepare a simple file with custom assignment of each individual sample
#name to one of the populations identified using PCA, STRUCTURE analysis etc. 
#Format: Sample_name Population_name
#But first: double check what sample names are present in the following VCF file:
cd $input/outliers/prfrag
filename=$input/outliers/pfrag/95m_contigs_unmasked_UK123_bi_filtered.recode.vcf
grep "#CHROM" $filename | head -1

#Resulting example file population assignment file: 
cat $input/outliers/pfrag/95m_contigs_unmasked_UK123_bi_filtered.recode.lst

#Now, need to prepare the configuration file. Copy the conversion script
pgdspid=/home/sobczm/bin/PGDSpider_2.0.5.2
#For diploid:
config=vcf_to_bayescan_diploid.spid
#For haploid:
#config=vcf_to_bayescan_haploid.spid
cp $pgdspid/$config ./
#Now, need to edit the conversion script to include the population assignment file
pop_assignment_file=$input/outliers/pfrag/95m_contigs_unmasked_UK123_bi_filtered.recode.lst
sed -i 's,^\(VCF_PARSER_POP_FILE_QUESTION=\).*,\1'"$pop_assignment_file"',' $config

#Convert to Bayescan input format with PGDSpider.
#For bigger input VCF file (>50 MB) use the script to convert to BayeScan input (will take a couple of hours):
#qsub $scripts/sub_pgdspider.sh $filename $config
#dos2unix ${filename%.vcf}.geno
java -jar -Xmx1024m -Xms512m $pgdspid/PGDSpider2-cli.jar -inputfile $filename -inputformat VCF -outputfile ${filename%.vcf}.geno -outputformat GESTE_BAYE_SCAN -spid $config
dos2unix ${filename%.vcf}.geno

############# 3) Bayescan analysis
#Bayescan run and plot the results in R.
qsub $scripts/sub_bayescan.sh ${filename%.vcf}.geno


##### Repeat of steps 2-3, but for haploid organism (V. inequalis)
#####
cd $input/outliers/vinequalis
filename=$input/outliers/vinequalis/Ash_farm_172_pacbio_contigs_unmasked_no118_bi_filtered.recode.vcf
grep "#CHROM" $filename | head -1

pgdspid=/home/sobczm/bin/PGDSpider_2.0.5.2
config=vcf_to_bayescan_haploid.spid
cp $pgdspid/$config ./

pop_assignment_file=$input/outliers/vinequalis/Ash_farm_172_pacbio_contigs_unmasked_no118_bi_filtered.recode.lst
sed -i 's,^\(VCF_PARSER_POP_FILE_QUESTION=\).*,\1'"$pop_assignment_file"',' $config

#For bigger input VCF file (>50 MB) use the script to convert to BayeScan input (will take a couple of hours)
qsub $scripts/sub_pgdspider.sh $filename $config
dos2unix ${filename%.vcf}.geno

qsub $scripts/sub_bayescan.sh ${filename%.vcf}.geno
