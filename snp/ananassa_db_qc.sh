#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/crown_rot_gwas/plate/strawberry_db
#Various checks of the ananassa snp chip database - QC. 
cd $input

#Dump all the relevant strawberry_sample tables to files to be analyzed.
a="SELECT id, pipeline_id, sample_id, alias_id, clustering, genotype FROM genotype"
echo $a | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >genotype
b="SELECT id, platform, probe_id, snp_id, marker_id FROM alias"
echo $b | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >alias
c="SELECT id, name, ref, alt, upstream, downstream FROM marker"
echo $c | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >marker
d="SELECT marker_id, tag FROM marker_tag"
echo $d | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >marker_tag
e="SELECT id, clone_id, file, path, type, date, batch FROM sample"
echo $e | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >sample
#Print the genotype QC table for all samples.
#Take a given marker_id and print the number of genotypes available for a given marker_id on istraw35 and istraw90, and alias_id, Affx codes associated with it on each chip.
python $scripts/db_qc.py alias genotype sample >qc_table
Rscript --vanilla $scripts/analyze_db.R
#Print my list of multiform and secondaryform markers.
#multiform means "Affx name with more than one marker id associated", secondaryform mean "Affx name with more than one marker id associated, *this* form does not match any istraw35 marker"
python $scripts/db_qc2.py alias 
#Distribution of number of probes per the same Affx ID in istraw35 and istraw90.
python $scripts/db_qc3.py alias 
#Output:
## istraw35
#1    30014
#2    4246

## istraw90
#1    67561
#2    25736
#4    4
#6    347
#12    1414

#For cases, where istraw90 and istraw35 both contain one unique Affx, check if they match - if they don't print them out. Also, check the total number of matches.
python $scripts/db_qc4.py alias 

#For istraw35/90 multiforms separately, need to check that no duplicates linking the same marker Affx id to the same alleles with the same flanking sequence (but MAYBE use different probes) and the same marker id - highlight them, and choose one with more genotypes available to proceed in each chip. 
python $scripts/db_qc5.py alias marker
#Positive control test using fake data to check if script works correctly, as obtained no dupliactes
python $scripts/db_qc5.py alias_test marker_test
