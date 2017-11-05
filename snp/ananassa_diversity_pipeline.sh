#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/diversity
cd $input 
#Select sample ids of individuals to be included in the analysis and extract their
#genotypes from the db. 
#Goal - create 3 datasets:
#A) istraw35 samples only
#B) istraw90 samples only
#C) joint analysis of istraw35 and istraw90 samples - use intersection of istraw35 and istraw90 markers

input_file=master_list_191917.txt
qsub $scripts/sub_ananassa_genotypes_db.sh $input_file ${input_file}.out

#Separate GWAS dataset by plate type (datasets A) and B) above)
awk -F"\t" '$4 == "istraw90" { print $0 }' ${input_file}.out  OFS='\t' >${input_file}_istraw90.out
awk -F"\t" '$4 == "istraw35" { print $0 }' ${input_file}.out  OFS='\t' >${input_file}_istraw35.out

#Save a copy of sample table for reference in scripts.
a="SELECT id, clone_id, file, path, type, batch FROM sample"
echo $a | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >sample.txt

#In cases where sample ids belong to the same cultivar (clone), will select the sample with the most genotypes.
for infile in ${input_file}.out ${input_file}_istraw35.out ${input_file}_istraw90.out
do
python $scripts/eliminate_duplicate_clones.py sample.txt $infile >${infile%.out}.lst
done

#Now need to re-run the sub_ananassa_genotypes_db.sh script to obtain genotypes
for infile in ${input_file}.lst ${input_file}_istraw35.lst ${input_file}_istraw90.lst
do
qsub $scripts/sub_ananassa_genotypes_db.sh $infile ${infile%.lst}.out
done

#Convert the resulting files to VCF format to be used in PLINK and TASSEL.
#Have a choice of 3 GFF files for SNP assignment to the position on the chromosome.
#A) vesca genome ver. 1.1 $input/istraw90_vesca_v1.1_snp_positions.gff3
#B) vesca genome ver. 2.0 $input/istraw90_vesca_v2.0_snp_positions.gff3
#C) ananassa genome $input/vesca2consensus_map_noambiguous_2017-08-22.gff --> assignment of positions of significantly fewer number of markers than in vesca.

#Need to check back with Rob to get the latest version of those GFF files, as he's working
#on improving them. 

#Here, using B) and C)
gff_file=istraw90_vesca_v2.0_snp_positions.gff3
for infile in ${input_file}.out ${input_file}_istraw35.out ${input_file}_istraw90.out
do
qsub $scripts/sub_ananassa_genotypes_vcf.sh $infile $gff_file
done
mv *.vcf vesca2.0

gff_file=vesca2consensus_map_noambiguous_2017-08-22.gff
for infile in ${input_file}.out ${input_file}_istraw35.out ${input_file}_istraw90.out
do
qsub $scripts/sub_ananassa_genotypes_vcf.sh $infile $gff_file
done
mv *.vcf ananassa

#First, plink does not allow chromosome names to start with a letter, so fix that, and sort by coordinates.
for infile in vesca2.0/${input_file}.out.vcf vesca2.0/${input_file}_istraw35.out.vcf vesca2.0/${input_file}_istraw90.out.vcf ananassa/${input_file}.out.vcf ananassa/${input_file}_istraw35.out.vcf ananassa/${input_file}_istraw90.out.vcf
do
cat $infile | sed 's/LG//' | sed 's/Unknown/0/' | awk 'NR<3{print $0;next}{print $0| "sort -k1,2"}'  >${infile%.vcf}_fix.vcf
done

#Convert VCF to plink format and retain only informative SNPs with minor allele freqeuncy (MAF)of at least 0.05.
#Log messages output by plink are useful and captured here, as they give info such as genotyping rate and number of variants/individuals before and after filtering.
for infile in vesca2.0/${input_file}.out_fix.vcf vesca2.0/${input_file}_istraw35.out_fix.vcf vesca2.0/${input_file}_istraw90.out_fix.vcf ananassa/${input_file}.out_fix.vcf ananassa/${input_file}_istraw35.out_fix.vcf ananassa/${input_file}_istraw90.out_fix.vcf
do
plink --allow-extra-chr --vcf $infile --maf 0.05 --recode --out ${infile%.vcf}_min05 > ${infile%.vcf}_min05.log
done

#Calculate the missing genotype rate for each marker. The results of this analysis can be found in files with ".lmiss" extension.
for infile in vesca2.0/${input_file}.out_fix_min05 vesca2.0/${input_file}_istraw35.out_fix_min05  vesca2.0/${input_file}_istraw90.out_fix_min05 ananassa/${input_file}.out_fix_min05 ananassa/${input_file}_istraw35.out_fix_min05 ananassa/${input_file}_istraw90.out_fix_min05
do
plink --allow-extra-chr --file $infile --missing --out ${infile%.ped} >${infile%.ped}.log
done

#Plot a histogram of the missing genotype rate to identify a threshold for extreme genotype failure rate. This can be carried out using the data in column 5 of the .lmiss file.
#But first, convert weird spacing between columns to tabs in PLINK output, as before.
for infile in vesca2.0/${input_file}.out_fix_min05 vesca2.0/${input_file}_istraw35.out_fix_min05  vesca2.0/${input_file}_istraw90.out_fix_min05 ananassa/${input_file}.out_fix_min05 ananassa/${input_file}_istraw35.out_fix_min05 ananassa/${input_file}_istraw90.out_fix_min05
do
cat ${infile}.lmiss | awk '{$1=$1;print}' OFS='\t' >temp
mv temp ${infile}.lmiss 
Rscript --vanilla $scripts/plot_missing_genotypes_plink.R ${infile}.lmiss
done

#Remove SNPs with more than a given % of missing data. Here, 50% and 20%. Need to convert to PLINK's BAM format at this point to run the command.
for infile in vesca2.0/${input_file}.out_fix_min05 vesca2.0/${input_file}_istraw35.out_fix_min05  vesca2.0/${input_file}_istraw90.out_fix_min05 ananassa/${input_file}.out_fix_min05 ananassa/${input_file}_istraw35.out_fix_min05 ananassa/${input_file}_istraw90.out_fix_min05
do
for per_missing in 0 0.01 0.05 0.1
do
    plink --file $infile --allow-extra-chr --geno $per_missing --make-bed --out ${infile}_${per_missing} >${infile}_${per_missing}.log
    cp ${infile%.ped}.map ${infile%.ped}_${per_missing}.map
done
done

#*From this point onwards, using only the subsets filtered for markers with low genotyping rates.*
#Identification of individuals with elevated missing data rates or outlying heterozygosity rate
for infile in vesca2.0/${input_file}.out_fix_min05 vesca2.0/${input_file}_istraw35.out_fix_min05  vesca2.0/${input_file}_istraw90.out_fix_min05 ananassa/${input_file}.out_fix_min05 ananassa/${input_file}_istraw35.out_fix_min05 ananassa/${input_file}_istraw90.out_fix_min05
do
for per_missing in 0 0.01 0.05 0.1
do
    plink --bfile ${infile}_${per_missing} --missing --allow-extra-chr --out ${infile}_${per_missing}  >${infile}_${per_missing}_post_filtering.log
    plink --bfile ${infile}_${per_missing} --het --allow-extra-chr --out ${infile}_${per_missing}
done
done

#This reates files with extension ".het", in which the third column denotes the observed number of homozygous genotypes [O(Hom)] and the fifth column denotes the number of nonmissing genotypes [N(NM)] per individual.
#The script below calculates the observed heterozygosity rate per individual using the formula (N(NM) − O(Hom))/N(NM). It then creates a graph in which the observed heterozygosity rate per individual is plotted on the x axis and the proportion of missing SNPs per individuals is plotted on the y axis.
#*Note* This script requires "geneplotter" and "tools" R libraries to be installed in your account.
for infile in vesca2.0/${input_file}.out_fix_min05 vesca2.0/${input_file}_istraw35.out_fix_min05  vesca2.0/${input_file}_istraw90.out_fix_min05 ananassa/${input_file}.out_fix_min05 ananassa/${input_file}_istraw35.out_fix_min05 ananassa/${input_file}_istraw90.out_fix_min05
do
for per_missing in 0 0.01 0.05 0.1
do
    Rscript --vanilla $scripts/imiss-vs-het.Rscript ${infile}_${per_missing}
done
done

#The graph (extension _het.pdf) suggests one outlier sample with very low heterozygosity, which is suspicious and needs to be removed. In order to find the sample (IID) and family (FID) id of that sample (both ids identical in our analyses), we need to open the table used to make the figure (extension _het.txt) and inspect the values in the columns (F_MISS - ie. Proportion of missing genotypes, meanHet - ie. Heterozygosity rate)

#In some cases, may be useful to calculate pairwise identity-by-descent (IBS - "DST" in the table output below) and PI_HAT (measure of identity-by-descent, but estimates only reliable using a big sample of individuals - not reliable here). We may then want to eliminate samples which are too closely related. As we have few samples and a lot of cultivars are inherently derived from a limited genetic pool, ignoring the results of this step here.
for infile in vesca2.0/${input_file}.out_fix_min05 vesca2.0/${input_file}_istraw35.out_fix_min05  vesca2.0/${input_file}_istraw90.out_fix_min05 ananassa/${input_file}.out_fix_min05 ananassa/${input_file}_istraw35.out_fix_min05 ananassa/${input_file}_istraw90.out_fix_min05
do
for per_missing in 0 0.01 0.05 0.1
do
    plink --bfile ${infile}_${per_missing} --allow-extra-chr --genome --out ${infile}_${per_missing}
done
done

#However, in certain cases we may want to filter our closely related individuals in a pair above a certain IBS threshold. 
#The code below looks at the individual call rates and outputs the IDs of the individual with the lowest call rate for subsequent removal, if desired.
threshold=0.90
for infile in vesca2.0/${input_file}.out_fix_min05 vesca2.0/${input_file}_istraw35.out_fix_min05  vesca2.0/${input_file}_istraw90.out_fix_min05 ananassa/${input_file}.out_fix_min05 ananassa/${input_file}_istraw35.out_fix_min05 ananassa/${input_file}_istraw90.out_fix_min05
do
for per_missing in 0 0.01 0.05 0.1
do
    perl $scripts/run-IBS-QC.pl ${infile}_${per_missing} $threshold > ${infile}_${per_missing}_to_filter
done
done

#Now, filter select individuals from the analysis (here samples related above 90% level by IBS).
for infile in vesca2.0/${input_file}.out_fix_min05 vesca2.0/${input_file}_istraw35.out_fix_min05  vesca2.0/${input_file}_istraw90.out_fix_min05 ananassa/${input_file}.out_fix_min05 ananassa/${input_file}_istraw35.out_fix_min05 ananassa/${input_file}_istraw90.out_fix_min05
do
for per_missing in 0 0.01 0.05 0.1
do
    plink --bfile ${infile}_${per_missing} --allow-extra-chr --remove ${infile}_${per_missing}_to_filter --make-bed --out ${infile}_${per_missing}_filtered >${infile}_${per_missing}.log
done
done

#To estimate population stratification, PLINK offers tools to cluster individuals into homogeneous subsets (which is achieved through complete linkage agglomerative clustering based on pair-wise IBS distance) and to perform classical MDS to visualize substructure and provide quantitative indices of population genetic variation that can be used as covariates in subsequent association analysis to control for stratification, instead of using discrete clusters.
for infile in vesca2.0/${input_file}.out_fix_min05 vesca2.0/${input_file}_istraw35.out_fix_min05  vesca2.0/${input_file}_istraw90.out_fix_min05 ananassa/${input_file}.out_fix_min05 ananassa/${input_file}_istraw35.out_fix_min05 ananassa/${input_file}_istraw90.out_fix_min05
do
for per_missing in 0 0.01 0.05 0.1
do
    plink --bfile ${infile}_${per_missing}_filtered --allow-extra-chr --read-genome ${infile}_${per_missing}.genome --cluster --mds-plot 4 --silent --out ${infile}_${per_missing}_filtered
done
done

#Convert weird spacing between columns to tabs in PLINK output 
for infile in vesca2.0/${input_file}.out_fix_min05 vesca2.0/${input_file}_istraw35.out_fix_min05  vesca2.0/${input_file}_istraw90.out_fix_min05 ananassa/${input_file}.out_fix_min05 ananassa/${input_file}_istraw35.out_fix_min05 ananassa/${input_file}_istraw90.out_fix_min05
do
for per_missing in 0 0.01 0.05 0.1
do
    cat ${infile}_${per_missing}_filtered.mds | awk '{$1=$1;print}' OFS='\t' >temp
    mv temp ${infile}_${per_missing}_filtered.mds 
done
done

#MDS plot based on Dimensions 1 and 2. 
#*Note*  Plotting requires R libraries "ggplot2", ""ggrepel" to be installed.
for infile in vesca2.0/${input_file}.out_fix_min05 vesca2.0/${input_file}_istraw35.out_fix_min05  vesca2.0/${input_file}_istraw90.out_fix_min05 ananassa/${input_file}.out_fix_min05 ananassa/${input_file}_istraw35.out_fix_min05 ananassa/${input_file}_istraw90.out_fix_min05
do
for per_missing in 0 0.01 0.05 0.1
do
    Rscript --vanilla $scripts/plot_plink_mds.R ${infile}_${per_missing}_filtered.mds
done
done

#Convert the filtered input files used in Plink to VCF so that can be used in other programs.
for infile in vesca2.0/${input_file}.out_fix_min05 vesca2.0/${input_file}_istraw35.out_fix_min05  vesca2.0/${input_file}_istraw90.out_fix_min05 ananassa/${input_file}.out_fix_min05 ananassa/${input_file}_istraw35.out_fix_min05 ananassa/${input_file}_istraw90.out_fix_min05
do
for per_missing in 0 0.01 0.05 0.1
do
plink --bfile ${infile}_${per_missing}_filtered --allow-extra-chr --recode vcf-iid --out ${infile}_${per_missing}_filtered
done
done

#Convert all PDFs to PNG.
for my_pdf in vesca2.0/*.pdf ananassa/*.pdf
do
convert -verbose -density 500 "${my_pdf}" "${my_pdf%.*}.png"
done