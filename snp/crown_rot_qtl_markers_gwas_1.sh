#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/crown_rot_gwas
cd $input

#First, running plink and TASSEL ithout imputation of genotypes - as we don't really know if we can impute them correctly, anyway - same problem as phasing.

#Copy over the VCF file with target individuals (retaining a couple of EMxFE individuals, but no sample duplicates).
cp /home/sobczm/popgen/snp/snp_chip/crown_rot_local/QTL_positions/sample_ids_crown_rot.out_nodup.vcf $input
#Convert VCF to plink and retain only informative SNPs with MAF of at least 0.05. 
#First, plink does not allow chromosome names to start with a letter, so fix that, and sort by coordinates.
cat sample_ids_crown_rot.out_nodup.vcf | sed 's/LG//' | sed 's/Unknown/0/' | awk 'NR<3{print $0;next}{print $0| "sort -k1,2"}'  >sample_ids_crown_rot.out_nodup_fix.vcf
#Total genotyping rate is 0.481915.
#49992 variants removed due to minor allele threshold(s)
#45073 variants and 106 people pass filters and QC.
plink --vcf sample_ids_crown_rot.out_nodup_fix.vcf --maf 0.05 --recode --out sample_ids_crown_rot.out_nodup_fix_min05
#Change sex to male from Unknown
awk '{$5 = "1"; print}' sample_ids_crown_rot.out_nodup_fix_min05.ped > temp
mv temp sample_ids_crown_rot.out_nodup_fix_min05.ped 
#Substitute the missing phenotype values for mean crown rot scores.
python $scripts/add_phenotype_ped.py sample_ids_crown_rot.out_nodup_fix_min05.ped crown_rot_scores.txt >sample_ids_crown_rot.out_nodup_fix_min05_pheno.ped
cp sample_ids_crown_rot.out_nodup_fix_min05.map sample_ids_crown_rot.out_nodup_fix_min05_pheno.map

##Following the tutorial from Anderson et al. (2010) Nature Protocols 9 (5) 1564-1573 doi:10.1038/nprot.2010.116
input_file=sample_ids_crown_rot.out_nodup_fix_min05_pheno
plink --file $input_file --missing --out raw-GWA-data
#Identification of individuals with elevated missing data rates or outlying heterozygosity rate
plink --file $input_file --het --out raw-GWA-data
#Creates the file raw-GWA-data.het, in which the third column denotes the observed number of homozygous genotypes [O(Hom)] and the fifth column denotes the number of nonmissing genotypes [N(NM)] per individual.
#Calculate the observed heterozygosity rate per individual using the formula (N(NM) − O(Hom))/N(NM). Create a graph in which the observed heterozygosity rate per individual is plotted on the x axis and the proportion of missing SNPs per individuals is plotted on the y axis.
R CMD BATCH $scripts/imiss-vs-het.Rscript
#The Graph suggests one outlier sample with very low heterozygosity, which is suspcious and needs to be removed. Adding individual 880 to the list of individuals to be removed.
##NB, it is the F. chiloensis sample so not surprising.
#Saving sample and family id of that sample in a list
echo '880 880' >>to_remove.txt

#The fourth column in the .imiss file (N_MISS) denotes the number of missing SNPs and the sixth column (F_MISS) denotes the proportion of missing SNPs per individual.

#Generate pairwise IBS for all pairs of individuals in the study
plink --file $input_file --genome --out raw-GWA-data
#identify all pairs of individuals with an IBD > 0.185 (third-degree relatives or closer)
perl $scripts/run-IBD-QC.pl raw-GWA-data
#The code looks at the individual call rates stored in raw-GWA-data.imiss and outputs the IDs of the individual with the lowest call rate to ‘fail-IBD-QC.txt’ for subsequent removal.
#NB: The IBD values are given in PI_HAT column, and in this case exceeding that threshold for an individual in at least one combination so useless here.

#To estimate population stratification, PLINK offers tools to clusterindividuals into homogeneous subsets (which is achieved
#through complete linkage agglomerative clustering based onpair-wise IBS distance) and to perform classical MDS to visualize
#substructure and provide quantitative indices of population genetic variation that can be used as covariates in subsequent
#association analysis to control for stratification, instead of using discrete clusters.
plink --file $input_file --read-genome raw-GWA-data.genome --cluster --mds-plot 4 --silent --out mds

#Convert wild spacing to tabs
cat mds.mds | awk '{$1=$1;print}' OFS='\t' >temp
mv temp mds.mds
#MDS plot based on Dim 1 and 2
R CMD BATCH $scripts/plot_plink_mds.R
#EmxFe samples cluster together but keeping them in the analysis.

#Remove select individuals from the analysis.
plink --file $input_file --remove to_remove.txt --make-bed --out clean-inds-GWA-data
#After removal total genotyping rate in remaining samples is 0.54107.

#calculate the missing genotype rate for each marker
plink --bfile clean-inds-GWA-data --missing --out clean-inds-GWA-data
#The results of this analysis can be found in clean-inds-GWA-data.lmiss.

#Plot a histogram of the missing genotype rate to identify a threshold for extreme genotype failure rate. 
#This can be carried out using the data in column 5 of the clean-inds-GWA-data.lmiss file 
#Plot the missing rate.
#Convert wild spacing to tabs
cat clean-inds-GWA-data.lmiss | awk '{$1=$1;print}' OFS='\t' >temp
mv temp clean-inds-GWA-data.lmiss 
R CMD BATCH $scripts/plot_missing_genotypes_plink.R
#Lots of missing data at around 50% missing data - probably those genotyped on istraw35 versus those on istraw90. Also many markers with negligable calling rates (<80%).

##Remove SNPs with more than 80% missing data ('relaxed' criterion)
plink --bfile clean-inds-GWA-data --geno 0.8 --make-bed --out clean-GWA-data_relaxed
#9026 variants removed due to missing genotype data (--geno).
#36047 variants and 104 people pass filters and QC.

##Remove SNPs with more than 50% missing data ('stringent' criterion)
plink --bfile clean-inds-GWA-data --geno 0.5 --make-bed --out clean-GWA-data_stringent
#26872 variants removed due to missing genotype data (--geno).
#18201 variants and 104 people pass filters and QC.

##For GWAS analysis, following the tutorial from Renteria et al. (2013)
#DOI 10.1007/978-1-62703-447-0_8
#Using allelic GWAS, with addditive model
plink --bfile clean-GWA-data_relaxed --assoc --qt-means --allow-no-sex --adjust --ci 0.95 --out clean-GWA-data_relaxed_1
plink --bfile clean-GWA-data_stringent --assoc --qt-means --allow-no-sex --adjust --ci 0.95 --out clean-GWA-data_stringent_1
#QQ plots
Rscript --vanilla $scripts/qq.plink.R clean-GWA-data_relaxed_1.qassoc "QQ plot"
Rscript --vanilla $scripts/qq.plink.R clean-GWA-data_stringent_1.qassoc "QQ plot" 
#QQ plots - see if using FDR (BH) p-values results in a better QQ plot.
cat clean-GWA-data_relaxed_1.qassoc.adjusted | sed 's/FDR_BH/P/' > clean-GWA-data_relaxed_1.qassoc.adjusted.qq
Rscript --vanilla $scripts/qq.plink.R clean-GWA-data_relaxed_1.qassoc.adjusted.qq "QQ plot" 
#Much worse plot - don't try again.
#QQ plots - see if only retaining SNPs with high % genotypes (95%) results in a better QQ plot.
plink --bfile clean-inds-GWA-data --geno 0.05 --make-bed --out clean-GWA-data_very_stringent
plink --bfile clean-GWA-data_very_stringent --assoc --qt-means  --allow-no-sex --adjust --ci 0.95 --out clean-GWA-data_very_stringent
Rscript --vanilla $scripts/qq.plink.R clean-GWA-data_very_stringent.qassoc "QQ plot"

#The same problem here.
#Dominant - linear egression
plink --bfile clean-GWA-data_relaxed --assoc --qt-means --linear dominant --allow-no-sex --adjust --ci 0.95 --out clean-GWA-data_relaxed_2
plink --bfile clean-GWA-data_stringent --assoc --qt-means --linear dominant --allow-no-sex --adjust --ci 0.95 --out clean-GWA-data_stringent_2
#Recessive - linear egression
plink --bfile clean-GWA-data_relaxed --assoc --qt-means --linear recessive --allow-no-sex --adjust --ci 0.95 --out clean-GWA-data_relaxed_3
plink --bfile clean-GWA-data_stringent --assoc --qt-means --linear recessive --allow-no-sex --adjust --ci 0.95 --out clean-GWA-data_stringent_3

#Correct for population stratification - include the MDS results as covariates.
#NB Covariates can only be used with the linear and logistic commands, so will use linear regression as a test for association.
awk '{print $1,$2,$4,$5}' mds.mds > mds_covar.txt
plink --bfile clean-GWA-data_relaxed --linear --allow-no-sex --covar mds_covar.txt --adjust --ci 0.95 --out clean-GWA-data_relaxed_strat
#QQ plots
Rscript --vanilla $scripts/qq.plink.R clean-GWA-data_relaxed_strat.assoc.linear "QQ plot"
plink --bfile clean-GWA-data_stringent --linear --allow-no-sex --covar mds_covar.txt --adjust --ci 0.95 --out clean-GWA-data_stringent_strat
#QQ plots
Rscript --vanilla $scripts/qq.plink.R clean-GWA-data_stringent_strat.assoc.linear "QQ plot"

##Manhattan plots
for results in  *qassoc
do
cat $results | awk '{$1=$1;print}' OFS='\t' >temp
cut -f2,1,3,9 temp >${results}_man
Rscript --vanilla $scripts/manhattan.R ${results}_man
done

for results in *assoc.linear
do
cat $results | awk '{$1=$1;print}' OFS='\t' >temp
cut -f2,1,3,12 temp >${results}_man
Rscript --vanilla $scripts/manhattan.R ${results}_man
done

#Convert all PDFs to PNG.
for my_pdf in *.pdf
do
convert -verbose -density 500 "${my_pdf}" "${my_pdf%.*}.png"
done

####Create filtered VCF files to be used as input in TASSEL.
#Remove individual 880
vcftools=/home/sobczm/bin/vcftools/bin
vcflib=/home/sobczm/bin/vcflib/bin
$vcflib/vcfremovesamples sample_ids_crown_rot.out_nodup_fix.vcf 880 >sample_ids_crown_rot.out_nodup_fix_no880.vcf
#Remove SNPs with less than 0.05 MAF.
##Remove SNPs with more than 80% missing data ('relaxed' criterion)
$vcftools/vcftools --vcf sample_ids_crown_rot.out_nodup_fix_no880.vcf --maf 0.05 --max-missing 0.20 --recode --out sample_ids_crown_rot.out_nodup_fix_no880_relaxed
##Remove SNPs with more than 50% missing data ('stringent' criterion)
$vcftools/vcftools --vcf sample_ids_crown_rot.out_nodup_fix_no880.vcf --maf 0.05 --max-missing 0.50 --recode --out sample_ids_crown_rot.out_nodup_fix_no880_stringent
#Sort VCF files.
$vcftools/vcf-sort sample_ids_crown_rot.out_nodup_fix_no880_relaxed.recode.vcf >sample_ids_crown_rot.out_nodup_fix_no880_relaxed_sorted.vcf
$vcftools/vcf-sort sample_ids_crown_rot.out_nodup_fix_no880_stringent.recode.vcf >sample_ids_crown_rot.out_nodup_fix_no880_stringent_sorted.vcf

#Analyze the location of candidate SNPs associated with crown rot resistance.
cd $input/candidates
#Are they within 1mbp of one of the QTLs?
for gff in Fa*1mbp.gff
do
intersectBed -wa -a gwas_snp_loci.bed -b $gff >gwas_snp_loci_vs_${gff%.gff}
done

#What genomic feature do they overlap with in vesca?
#No genes
intersectBed -wb -a gwas_snp_loci.bed -b Fragaria_vesca_v1.1.a2.mrna.gff3 >gwas_snp_loci_vs_Fragaria_vesca_v1.1.a2.mrna

#Are they within 1mbp of one of the predicted R genes in vesca?
for a in vesca*.gff3
do
intersectBed -wb -a gwas_snp_loci_1mbp.bed -b $a >charlotte_1mbp2_${a%.gff3}_qtl_overlap
done

for a in vesca*.gff3
do
intersectBed -wb -a  gwas_snp_loci_1ookbp.bed -b $a >gwas_snp_loci_100kbp_${a%.gff3}_overlap
done

for a in vesca*.gff3
do
intersectBed -wb -a gwas_snp_loci.bed -b $a >gwas_snp_loci_${a%.gff3}_overlap
done

###Call the presence and absence of those SNP in the correct set of strawberry breeding material that is need for decisions at the meeting - just 3 genotyped cultivars EM2157 (DB sample id: 2194, crown rot score: 4.4)
#EM2192 (DB sample id: 2013, crown rot score 3.95) EM2199 (DB sample id: 2061, crown rot score 3.2). Max crown rot score: 8
cd $input/targets
echo "2013" > sbc_samples.txt
echo "2061" >> sbc_samples.txt
echo "2194" >> sbc_samples.txt
python $scripts/ananassa_genotypes_db.py sbc_samples.txt sbc_samples.out

#Output the genotypes in the VCF format with locations substituted  according to map positions relative to vesca 1.1. 
python $scripts/ananassa_genotypes_vcf.py sbc_samples.out istraw90_vesca_v1.1_snp_positions.gff3