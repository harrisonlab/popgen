#!/bin/bash
input=/home/sobczm/popgen/snp/snp_chip/haplotypes/comparison
scripts=/home/sobczm/bin/popgen/snp

#Copy tha map-based phasing results.

#Extract results from matching individuals - haplotyping.
for a in *_map.ped
do
python $scripts/overall_haplotype_similarity.py $a ${a%_map.ped}_impute.ped ${a%.ped}.info ${a%_map.ped}_impute.info
python $scripts/short_identical_haplotypes.py $a ${a%_map.ped}_impute.ped ${a%.ped}.info ${a%_map.ped}_impute.info
python $scripts/haplotype_length_distribution.py $a ${a%_map.ped}_impute.ped ${a%.ped}.info ${a%_map.ped}_impute.info
done

#Plot results using ggplot2. 
for a in *_overall.stat
do
Rscript --vanilla $scripts/overall_haplotype_similarity_fig.R $a 8 4
done

for a in *_matching_fragments.stat
do
Rscript --vanilla $scripts/short_identical_haplotypes_fig1.R $a 20 4
Rscript --vanilla $scripts/short_identical_haplotypes_fig2.R $a 8 4
done

for a in *_match_length.stat
do
Rscript --vanilla $scripts/haplotype_length_distribution_fig.R $a 20 4
done

#Copy haplotyping results from shapeit
haplotyping_shapeit1=/home/vicker/octoploid_mapping/haplotyping4/shapeit_all_samples
haplotyping_shapeit2=/home/vicker/octoploid_mapping/haplotyping4/shapeit_filtered_progeny
cp -r $haplotyping_shapeit1 $haplotyping_shapeit2 $input
#Copy haplotyping results from map
map=/home/vicker/octoploid_mapping/consensus_map4
cp -r $map/popn_EMxFE/map/haplotypes popn_EMxFE
cp -r $map/popn_RGxHA/map/haplotypes popn_RGxHA
cp -r $map/popn_FLxCH/map/haplotypes popn_FLxCH

###Convert from the custom map-based hapotype format to .ped and .info
#Obtain mapping of filenames to sample_id.
echo "SELECT id, file FROM sample" | mysql -u marias -h mongo -D strawberry_copy -p$(cat /home/sobczm/bin/mysql_sample_database/login)> id_filename.tsv

python $scripts/map_to_ped.py 1A_haplotypes.csv id_filename.tsv vesca2consensus_map_noduplicates_2017-05-17.csv
python $scripts/overall_haplotype_similarity.py 1A_phased.ped 1A.shapeit.ped.phased 1A.info 1A.shapeit.info
python $scripts/short_identical_haplotypes.py 1A_phased.ped 1A.shapeit.ped.phased 1A.info 1A.shapeit.info
python $scripts/haplotype_length_distribution.py 1A_phased.ped 1A.shapeit.ped.phased 1A.info 1A.shapeit.info