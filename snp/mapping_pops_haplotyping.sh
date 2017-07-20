#!/bin/bash

#
# dump genotype data for haplotyping all samples
# datadir: /home/sobczm/popgen/snp/snp_chip/haplotypes/comparison/haplotyping4_mapping_pops
#

set -eu

input=/home/vicker/octoploid_mapping/haplotyping4
datadir=/home/sobczm/popgen/snp/snp_chip/haplotypes/comparison/haplotyping4_mapping_pops
scripts=/home/sobczm/bin/popgen/snp
export PATH=${scripts}:${PATH}
export PATH=/home/vicker/programs/shapeit.v2.r837.GLIBCv2.12.Linux.static/bin:${PATH}
export PATH=/home/vicker/programs/plink-1.90beta:${PATH}
source ${scripts}/haplotyping_all_samples_funcs.sh
export NODES="-N2,5,6,7,8,9,10"
export MAXJOBS="-L14"


#Retrieve file with all the FLxCH genotypes dumped.
mysql -B -u marias -h mongo -p$(cat /home/sobczm/bin/mysql_sample_database/login) -D strawberry_samples \
    < $scripts/extract_sql_flch.sql \
    | gzip \
    > genotypes-flch.tsv.gz

#Retrieve file with all the RGxHA genotypes dumped.
mysql -B -u marias -h mongo -p$(cat /home/sobczm/bin/mysql_sample_database/login) -D strawberry_samples \
    < $scripts/extract_sql_rgha.sql \
    | gzip \
    > genotypes-rgha.tsv.gz
#Retrieve file with all the EMxFE genotypes dumped.
mysql -B -u marias -h mongo -p$(cat /home/sobczm/bin/mysql_sample_database/login) -D strawberry_samples \
    < $scripts/extract_sql_emfe.sql \
    | gzip \
    > genotypes-emfe.tsv.gz

#select best pipeline per sample
mysql -B -u marias -h mongo -p$(cat /home/sobczm/bin/mysql_sample_database/login)  -D strawberry_samples \
    < $scripts/find_best_pipeline.sql \
    | gzip \
    > best_pipeline.tsv.gz

#filter out calls not from best pipeline, drop null calls
for genotypes in genotypes-flch.tsv.gz genotypes-rgha.tsv.gz genotypes-emfe.tsv.gz
do
zcat $genotypes \
    | python $scripts/filter_pipeline.py best_pipeline.tsv.gz \
    | gzip \
    > ${genotypes%.tsv.gz}_best_pipeline.tsv.gz
done 

#filter out calls from flagged markers (multiallelic, sexerror, crosslg)
for genotypes in *_best_pipeline.tsv.gz
do
zcat $genotypes \
    | python $scripts/filter_markers.py \
    | gzip \
    > ${genotypes%_best_pipeline.tsv.gz}_filtered_markers.tsv.gz
done

#convert to tablular csv format
for genotypes in *_filtered_markers.tsv.gz
do
python $scripts/convert_to_table.py $genotypes \
    | gzip \
    > ${genotypes%_filtered_markers.tsv.gz}_table_all_samples.csv.gz
done

#reformat ready for shapeit
for genotypes in *_table_all_samples.csv.gz
do
convert2shapeit \
/home/groups/harrisonlab/project_files/fragaria_x_ananassa/octoseq/maps/consensus/vesca2consensus_map_noduplicates_2017-05-17.csv \
$genotypes ${genotypes%_table_all_samples.csv.gz}_shapeit
done

for genotypes in *_shapeit
do
impute_haplotypes $genotypes
estimate_haploblocks $genotypes
extract_haplotypes2ped $genotypes
done

##Comparison of shapeit and map-based haplotypes.
input=/home/sobczm/popgen/snp/snp_chip/haplotypes/comparison
#For EMxFE
#conversion of map-based haplos to ped
mkdir $datadir/genotypes-emfe_shapeit/comparison 
cd $datadir/genotypes-emfe_shapeit/comparison 
cp $input/id_filename.tsv $input/popn_EMxFE/1A_haplotypes.csv $input/vesca2consensus_map_noduplicates_2017-05-17.csv ./
python $scripts/map_to_ped.py 1A_haplotypes.csv id_filename.tsv vesca2consensus_map_noduplicates_2017-05-17.csv
python $scripts/overall_haplotype_similarity.py 1A_phased.ped ../1A.ped.phased 1A.info ../1A.info 

#For FLxCH
#conversion of map-based haplos to ped
mkdir $datadir/genotypes-flch_shapeit/comparison
cd $datadir/genotypes-flch_shapeit/comparison 
cp $input/id_filename.tsv $input/popn_FLxCH/1A_haplotypes.csv $input/vesca2consensus_map_noduplicates_2017-05-17.csv ./
python $scripts/map_to_ped.py 1A_haplotypes.csv id_filename.tsv vesca2consensus_map_noduplicates_2017-05-17.csv
python $scripts/overall_haplotype_similarity.py 1A_phased.ped ../1A.ped.phased 1A.info ../1A.info 

#For RGxH
#A#conversion of map-based haplos to ped
mkdir $datadir/genotypes-rgha_shapeit/comparison
cd $datadir/genotypes-rgha_shapeit/comparison 
cp $input/id_filename.tsv $input/popn_RGxHA/1A_haplotypes.csv $input/vesca2consensus_map_noduplicates_2017-05-17.csv ./
python $scripts/map_to_ped.py 1A_haplotypes.csv id_filename.tsv vesca2consensus_map_noduplicates_2017-05-17.csv
python $scripts/overall_haplotype_similarity.py 1A_phased.ped ../1A.ped.phased 1A.info ../1A.info 
