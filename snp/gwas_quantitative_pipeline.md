# Using data with with continuous phenotype values
```
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/cr_gwas
cd $input
```
## Data preparation
Select sample ids of individuals to be included in the analysis and extract their
genotypes from the db. 
Goal - create 3 datasets:
A) istraw35 samples only
B) istraw90 samples only
C) joint analysis of istraw35 and istraw90 samples - use intersection of istraw35 and istraw90 markers
```
input_file=cr_list.txt
qsub $scripts/sub_ananassa_genotypes_db.sh $input_file ${input_file}.out
```
Separate GWAS dataset by plate type (datasets A) and B) above)
```
awk -F"\t" '$4 == "istraw90" { print $0 }' ${input_file}.out  OFS='\t' >${input_file}_istraw90.out
awk -F"\t" '$4 == "istraw35" { print $0 }' ${input_file}.out  OFS='\t' >${input_file}_istraw35.out
```
In cases where sample ids belong to the same cultivar (clone), will select the sample with the most genotypes.s
Save a copy of sample table for reference in scripts.
```
a="SELECT id, clone_id, file, path, type, batch FROM sample"
echo $a | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >sample.txt
```
```
for infile in ${input_file}.out ${input_file}_istraw35.out ${input_file}_istraw90.out
do
python $scripts/eliminate_duplicate_clones.py sample.txt $infile >${infile%.out}.lst
done
```
Now need to re-run the sub_ananassa_genotypes_db.sh script to obtain genotypes
```
for infile in ${input_file}.lst ${input_file}_istraw35.lst ${input_file}_istraw90.lst
do
qsub $scripts/sub_ananassa_genotypes_db.sh $infile ${infile%.lst}.out
done
```

Convert the resulting files to VCF format to be used in PLINK and TASSEL.
Have a choice of 3 GFF files for SNP assignment to the position on the chromosome.

A) vesca genome ver. 1.1 $input/istraw90_vesca_v1.1_snp_positions.gff3

B) vesca genome ver. 2.0 $input/istraw90_vesca_v2.0_snp_positions.gff3

C) ananassa genome $input/vesca2consensus_map_noambiguous_2017-08-22.gff --> assignment of positions of significantly fewer number of markers than in vesca.

Need to check back with Rob to get the latest version of those GFF files, as he's working
on improving them. 

Here, using A) as this version has been used in the past in Charlotte's QTL analysis.
```
gff_file=$input/istraw90_vesca_v1.1_snp_positions.gff3
for infile in ${input_file}.out ${input_file}_istraw35.out ${input_file}_istraw90.out
do
qsub $scripts/sub_ananassa_genotypes_vcf.sh $infile $gff_file
done
```
## Data analysis with PLINK
First, plink does not allow chromosome names to start with a letter, so fix that, and sort by coordinates.
```
for infile in ${input_file}.out.vcf ${input_file}_istraw35.out.vcf ${input_file}_istraw90.out.vcf
do
cat $infile | sed 's/LG//' | sed 's/Unknown/0/' | awk 'NR<3{print $0;next}{print $0| "sort -k1,2"}'  >${infile%.vcf}_fix.vcf
done
```
Convert VCF to plink format and retain only informative SNPs with minor allele freqeuncy (MAF) of at least 0.05.
Log messages output by plink are useful and captured here, as they give info such as genotyping rate and number of variants/individuals before and after filtering.
```
for infile in ${input_file}.out_fix.vcf ${input_file}_istraw35.out_fix.vcf ${input_file}_istraw90.out_fix.vcf
do
plink --vcf $infile --maf 0.05 --recode --out ${infile%.vcf}_min05 > ${infile%.vcf}_min05.log
done
```
Change sex to male from Unknown in the input files to make GWAS analysis possible.
```
for infile in ${input_file}.out_fix_min05.ped ${input_file}_istraw35.out_fix_min05.ped ${input_file}_istraw90.out_fix_min05.ped
do
awk '{$5 = "1"; print}' $infile > temp
mv temp $infile
done
```
Substitute the missing phenotype values for mean crown rot scores, or other continous phenotype values. Requires an input table with sample_id in the first column and phenotype value in the second column (see example file below)
```
phenotype_file=crown_rot_scores.txt
for infile in ${input_file}.out_fix_min05.ped ${input_file}_istraw35.out_fix_min05.ped ${input_file}_istraw90.out_fix_min05.ped
do
python $scripts/add_phenotype_ped.py $infile $phenotype_file >${infile%.ped}_pheno.ped 
cp ${infile%.ped}.map ${infile%.ped}_pheno.map
done
```
Identification of individuals with elevated missing data rates or outlying heterozygosity rate
```
for infile in ${input_file}.out_fix_min05_pheno ${input_file}_istraw35.out_fix_min05_pheno  ${input_file}_istraw90.out_fix_min05_pheno
do
plink --file $infile --missing --out ${infile%.ped}
plink --file $infile --het --out ${infile%.ped}
done
```

This reates files with extension ".het", in which the third column denotes the observed number of homozygous genotypes [O(Hom)] and the fifth column denotes the number of nonmissing genotypes [N(NM)] per individual.
The script below calculates the observed heterozygosity rate per individual using the formula (N(NM) − O(Hom))/N(NM). It then creates a graph in which the observed heterozygosity rate per individual is plotted on the x axis and the proportion of missing SNPs per individuals is plotted on the y axis.
*Note* This script requires "geneplotter" and "tools" R libraries to be installed in your account.
```
for infile in ${input_file}.out_fix_min05_pheno ${input_file}_istraw35.out_fix_min05_pheno  ${input_file}_istraw90.out_fix_min05_pheno
do
Rscript --vanilla $scripts/imiss-vs-het.Rscript $infile
done
```
The graph (extension _het.pdf) suggests one outlier sample with very low heterozygosity, which is suspicious and needs to be removed. In order to find the sample (IID) and family (FID) id of that sample (both ids identical in our analyses), we need to open the table used to make the figure (extension _het.txt) and inspect the values in the columns (F_MISS - ie. Proportion of missing genotypes, meanHet - ie. Heterozygosity rate)

Here, adding individual 880 to the list of individuals to be removed.
NB, it is the F. chiloensis sample so not surprising.

Saving sample (IID) and family id (FID) of that sample in a list, to be removed later - one sample per line.
```
echo '880 880' >>to_remove.txt
```
In some cases, may be useful to calculate pairwise identity-by-descent (IBS - "DST" in the table output below) and PI_HAT (measure of identity-by-descent, but estimates only reliable using a big sample of individuals - not reliable here). We may then want to eliminate samples which are too closely related. As we have few samples and a lot of cultivars are inherently derived from a limited genetic pool, ignoring the results of this step here.
```
for infile in ${input_file}.out_fix_min05_pheno ${input_file}_istraw35.out_fix_min05_pheno  ${input_file}_istraw90.out_fix_min05_pheno
do
plink --file $infile --genome --out ${infile%.ped} > ${infile}.log
done
```
However, in certain cases we may want to filter our closely related individuals in a pair above a certain IBS threshold. 
The code below looks at the individual call rates and outputs the IDs of the individual with the lowest call rate for subsequent removal, if desired.
Example below - however, not going to filter out the indviduals output below in our analysis.
```
threshold=0.95
for infile in ${input_file}.out_fix_min05_pheno ${input_file}_istraw35.out_fix_min05_pheno  ${input_file}_istraw90.out_fix_min05_pheno
do
perl $scripts/run-IBS-QC.pl $infile $threshold > ${infile}_to_filter
done
```

To estimate population stratification, PLINK offers tools to cluster individuals into homogeneous subsets (which is achieved through complete linkage agglomerative clustering based on pair-wise IBS distance) and to perform classical MDS to visualize substructure and provide quantitative indices of population genetic variation that can be used as covariates in subsequent association analysis to control for stratification, instead of using discrete clusters.

```
for infile in ${input_file}.out_fix_min05_pheno ${input_file}_istraw35.out_fix_min05_pheno  ${input_file}_istraw90.out_fix_min05_pheno
do
plink --file $infile --read-genome ${infile}.genome --cluster --mds-plot 4 --silent --out ${infile}
done
```

Convert weird spacing between columns to tabs in PLINK output 
```
for infile in ${input_file}.out_fix_min05_pheno ${input_file}_istraw35.out_fix_min05_pheno  ${input_file}_istraw90.out_fix_min05_pheno
do
cat ${infile}.mds | awk '{$1=$1;print}' OFS='\t' >temp
mv temp ${infile}.mds
done
```
MDS plot based on Dimensions 1 and 2. 
*Note*  Plotting requires R libraries "ggplot2", ""ggrepel" to be installed.
```
for infile in ${input_file}.out_fix_min05_pheno ${input_file}_istraw35.out_fix_min05_pheno  ${input_file}_istraw90.out_fix_min05_pheno
do
Rscript --vanilla $scripts/plot_plink_mds.R ${infile}.mds
done
```
EmxFe samples cluster together but keeping them in the analysis as we are so short for samples in this analysis.

Now, filter select individuals from the analysis (here only sample 880 - F. chiloensis, ignoring samples with high IBS).
For this step, need to convert to internal BAM file format used by PLINK.
```
for infile in ${input_file}.out_fix_min05_pheno ${input_file}_istraw35.out_fix_min05_pheno  ${input_file}_istraw90.out_fix_min05_pheno
do
plink --file $infile --remove to_remove.txt --make-bed --out ${infile}_filtered 
done
```

Calculate the missing genotype rate for each marker. The results of this analysis can be found in files with ".lmiss" extension.
```
for infile in ${input_file}.out_fix_min05_pheno_filtered ${input_file}_istraw35.out_fix_min05_pheno_filtered  ${input_file}_istraw90.out_fix_min05_pheno_filtered 
do
plink --bfile $infile --missing --out ${infile%.bam}
done
```

Plot a histogram of the missing genotype rate to identify a threshold for extreme genotype failure rate. This can be carried out using the data in column 5 of the .lmiss file.
But first, convert weird spacing between columns to tabs in PLINK output, as before.
```
for infile in ${input_file}.out_fix_min05_pheno_filtered ${input_file}_istraw35.out_fix_min05_pheno_filtered  ${input_file}_istraw90.out_fix_min05_pheno_filtered 
do
cat ${infile}.lmiss | awk '{$1=$1;print}' OFS='\t' >temp
mv temp ${infile}.lmiss 
Rscript --vanilla $scripts/plot_missing_genotypes_plink.R ${infile}.lmiss
done
```

Remove SNPs with more than 80% missing data ('relaxed' criterion)
```
for infile in ${input_file}.out_fix_min05_pheno_filtered ${input_file}_istraw35.out_fix_min05_pheno_filtered  ${input_file}_istraw90.out_fix_min05_pheno_filtered 
do
plink --bfile $infile --geno 0.8 --make-bed --out ${infile}_relaxed >${infile}_relaxed.log
done
```

Remove SNPs with more than 50% missing data ('stringent' criterion)
```
for infile in ${input_file}.out_fix_min05_pheno_filtered ${input_file}_istraw35.out_fix_min05_pheno_filtered  ${input_file}_istraw90.out_fix_min05_pheno_filtered 
do
plink --bfile $infile --geno 0.5 --make-bed --out ${infile}_stringent >${infile}_strinent.log
done
```
