# Using data with with continous phenotype values
```
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/cr_gwas
```
## Data preparation
Select sample ids of individuals to be included in the analysis and extract their
genotypes from the db. 
```
input_file=cr_list.txt
qsub $scripts/sub_ananassa_genotypes_db.sh $input_file ${input_file}.out
```

In cases where sample ids belong to the same cultivar (clone), will select the sample with the most genotypes.
python $scripts/eliminate_duplicate_clones.py

Save a copy of sample table for reference in scripts.
```
a="SELECT id, clone_id, file, path, type, batch FROM sample"
echo $a | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login) >sample.txt
```

A) istraw35 samples only
B) istraw90 samples only
C) joint analysis of istraw35 and istraw90 samples - use intersection of istraw35 and istraw90 markers