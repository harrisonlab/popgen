#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/crown_rot_gwas/plate
cd $input

#Re-extract the GWAS dataset (without sample 880), as well as the SBC sample dataset.
python $scripts/ananassa_genotypes_db.py sample_ids_crown_rot_nodup.txt sample_ids_crown_rot_nodup.out
python $scripts/ananassa_genotypes_db.py sbc_candidates.txt sbc_candidates.out
a="SELECT id, clone_id, file, path, type, batch FROM sample WHERE id in (2194,2013,2061,2163,1943,1991,2087,2218,2185,1949,2047,2049,1907,2230,2053,2206,2208,1931,2027,2123,2171)"
echo $a | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login)> sbc_samples.txt
b="SELECT id, clone_id, file, path, type, batch FROM sample WHERE id in (1468,1341,1196,1203,1420,1255,859,1114,1113,2234,2208,2206,2202,2200,2198,2194,2192,2190,2171,2167,2165,2163,2161,2159,2157,2155,2153,2151,2123,2117,2111,2107,2105,2103,2067,2065,2063,2061,2059,2057,2019,2017,2015,2013,2011,2009,1971,1969,1967,1963,1931,1925,1923,1921,1919,1917,1915,1913,1881,1879,1875,1873,1871,1869,1867,1387,1218,1181,1180,1112,1111,1110,1108,1104,1103,1102,1101,1098,1096,1095,1094,1093,1092,1091,1089,1088,1087,894,892,878,869,860,851,850,842,831,808,677,675,644,601,595,587,577,568,554,551)"
echo $b | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login)> sample_ids_crown_rot_nodup_samples.txt

#Modify the path to simple plate name.
sed -i 's/\/home\/groups\/harrisonlab\/raw_data\/raw_celfiles\/istraw35_1plate_P160793_B1/istraw35_1plate_P160793_B1/' sbc_samples.txt
sed -i 's/\/home\/groups\/harrisonlab\/raw_data\/raw_celfiles\/istraw35_1plate_P160793_B1/istraw35_1plate_P160793_B1/' sample_ids_crown_rot_nodup_samples.txt 
sed -i 's/\/home\/groups\/harrisonlab\/raw_data\/raw_celfiles\/istraw90_2plates_P150761Q151138_20160722\/symlinks/istraw90_2plates_P150761Q151138_20160722/' sample_ids_crown_rot_nodup_samples.txt 
sed -i 's/\/home\/groups\/harrisonlab\/raw_data\/raw_celfiles\/istraw90_4plates_rosbreed\/symlinks\/known/istraw90_4plates_rosbreed/' sample_ids_crown_rot_nodup_samples.txt 
sed -i 's/\/home\/groups\/harrisonlab\/raw_data\/raw_celfiles\/istraw90_6plates_RGxHA_EMxFE_FLxCH_BSxEL_etc\/symlinks/istraw90_6plates_RGxHA_EMxFE_FLxCH_BSxEL_etc/' sample_ids_crown_rot_nodup_samples.txt 
#Modify the MDS dataframe to add plate name and plot
R CMD BATCH $scripts/plot_plink_mds2.R

#Separate GWAS dataset by plate type
awk -F"\t" '$4 == "istraw90" { print $0 }' sample_ids_crown_rot_nodup.out  OFS='\t' >sample_ids_crown_rot_nodup_istraw90.out
awk -F"\t" '$4 == "istraw35" { print $0 }' sample_ids_crown_rot_nodup.out  OFS='\t' >sample_ids_crown_rot_nodup_istraw35.out

python $scripts/ananassa_genotypes_vcf.py sample_ids_crown_rot_nodup_istraw90.out istraw90_vesca_v1.1_snp_positions.gff3
python $scripts/ananassa_genotypes_vcf.py sample_ids_crown_rot_nodup_istraw35.out istraw90_vesca_v1.1_snp_positions.gff3

for vcf in sample_ids_crown_rot_nodup_istraw35.out.vcf sample_ids_crown_rot_nodup_istraw90.out.vcf
do
cat $vcf | sed 's/LG//' | sed 's/Unknown/0/' | awk 'NR<3{print $0;next}{print $0| "sort -k1,2"}'  >${vcf%.vcf}_fix.vcf
done
#54 people (0 males, 0 females, 54 ambiguous) loaded from .fam.
#Total genotyping rate is 0.999999.
#5911 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#16385 variants and 54 people pass filters and QC.
plink --vcf sample_ids_crown_rot_nodup_istraw35.out_fix.vcf --maf 0.05 --recode --out sample_ids_crown_rot_nodup_istraw35.out_fix_min05
#50 people (0 males, 0 females, 50 ambiguous) loaded from .fam.
#Total genotyping rate is 0.746874.
#48961 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#45510 variants and 50 people pass filters and QC.
plink --vcf sample_ids_crown_rot_nodup_istraw90.out_fix.vcf --maf 0.05 --recode --out sample_ids_crown_rot_nodup_istraw90.out_fix_min05

#Change sex to male from Unknown
for p in 35 90
do
awk '{$5 = "1"; print}' sample_ids_crown_rot_nodup_istraw${p}.out_fix_min05.ped > temp
mv temp sample_ids_crown_rot_nodup_istraw${p}.out_fix_min05.ped
#Substitute the missing phenotype values for mean crown rot scores.
python $scripts/add_phenotype_ped.py sample_ids_crown_rot_nodup_istraw${p}.out_fix_min05.ped crown_rot_scores.txt >sample_ids_crown_rot_nodup_istraw${p}.out_fix_min05_pheno.ped
cp sample_ids_crown_rot_nodup_istraw${p}.out_fix_min05.map sample_ids_crown_rot_nodup_istraw${p}.out_fix_min05_pheno.map
done
#Follow the rest of the QC part of the GWAS pipeline in different subdirs
mkdir istraw35
mv *istraw35* ./istraw35

mkdir istraw90
mv *istraw90* ./istraw90

for p in 35 90
do
cd istraw${p}
input_file=sample_ids_crown_rot_nodup_istraw${p}.out_fix_min05_pheno
plink --file $input_file --missing --out raw-GWA-data
plink --file $input_file --het --out raw-GWA-data
R CMD BATCH $scripts/imiss-vs-het.Rscript
plink --file $input_file --genome --out raw-GWA-data
perl $scripts/run-IBD-QC.pl raw-GWA-data
plink --file $input_file --read-genome raw-GWA-data.genome --cluster --mds-plot 4 --silent --out mds
cat mds.mds | awk '{$1=$1;print}' OFS='\t' >temp
mv temp mds.mds
#MDS plot based on Dim 1 and 2
R CMD BATCH $scripts/plot_plink_mds.R
plink --file $input_file --missing --out clean-inds-GWA-data
cat clean-inds-GWA-data.lmiss | awk '{$1=$1;print}' OFS='\t' >temp
mv temp clean-inds-GWA-data.lmiss 
R CMD BATCH $scripts/plot_missing_genotypes_plink.R
cd ../
done

#Get the genotypes for the select 5 markers for SBC samples.
python $scripts/ananassa_genotypes_vcf.py sbc_candidates.out istraw90_vesca_v1.1_snp_positions.gff3
cat sbc_candidates.out.vcf | sed 's/LG//' | sed 's/Unknown/0/' | awk 'NR<3{print $0;next}{print $0| "sort -k1,2"}'  >sbc_candidates.out2.vcf 
#File given to the 'extract' option contains the list of markers we are interested in - here the 5 markers predictive of crown rot resistance.
plink --vcf sbc_candidates.out2.vcf --extract $input/istraw_35_outliers.txt --recode A --out sbc_samples_outliers
cat sbc_samples_outliers.raw | awk '{$1=$1;print}' OFS='\t' >temp
mv temp sbc_samples_outliers.raw

#Get the genotypes for the select 5 markers for all EM and EMR samples.
python $scripts/ananassa_genotypes_db.py EM_EMR.txt EM_EMR.out
python $scripts/ananassa_genotypes_vcf.py EM_EMR.out istraw90_vesca_v1.1_snp_positions.gff3
cat EM_EMR.out.vcf | sed 's/LG//' | sed 's/Unknown/0/' | awk 'NR<3{print $0;next}{print $0| "sort -k1,2"}'  >EM_EMR.out2.vcf 
#File given to the 'extract' option contains the list of markers we are interested in - here the 5 markers predictive of crown rot resistance.
plink --vcf EM_EMR.out2.vcf --extract $input/istraw_35_outliers.txt --recode A --out EM_EMR_outliers
cat EM_EMR_outliers.raw | awk '{$1=$1;print}' OFS='\t' >temp
mv temp EM_EMR_outliers.raw

#Get the genotypes for the select 5 markers for all samples in the master strawberry spreadsheet.
a=all_cultivars_ids
python $scripts/ananassa_genotypes_db.py $a.txt $a.out
python $scripts/ananassa_genotypes_vcf.py $a.out istraw90_vesca_v1.1_snp_positions.gff3
cat $a.out.vcf | sed 's/LG//' | sed 's/Unknown/0/' | awk 'NR<3{print $0;next}{print $0| "sort -k1,2"}'  >$a.out2.vcf 
#File given to the 'extract' option contains the list of markers we are interested in - here the 5 markers predictive of crown rot resistance.
plink --vcf $a.out2.vcf --extract $input/istraw_35_outliers.txt --recode A --out ${a}_outliers
cat ${a}_outliers.raw | awk '{$1=$1;print}' OFS='\t' >temp
mv temp ${a}_outliers.raw2

#Repeated after fixed the bug in the db table to VCF conversion script.

#Print the genotype QC table for all samples.
cd strawberry_db
python $scripts/db_qc.py alias genotype sample >qc_table
python $scripts/db_qc2.py alias 
python $scripts/db_qc3.py alias 
python $scripts/db_qc4.py alias 
