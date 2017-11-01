# Using data with with continous phenotype values
```
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/cr_gwas
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
In cases where sample ids belong to the same cultivar (clone), will select the sample with the most genotypes.
```
for infile in ${input_file}.out ${input_file}_istraw35.out ${input_file}_istraw90.out
do
python $scripts/eliminate_duplicate_clones.py $infile
done
```
Now need to re-run the sub_ananassa_genotypes_db.sh s
```
for infile in ${input_file}.out ${input_file}_istraw35.out ${input_file}_istraw90.out
do
python $scripts/eliminate_duplicate_clones.py $infile
done
```
Now need to re-run the sub_ananassa_genotypes_db.sh s