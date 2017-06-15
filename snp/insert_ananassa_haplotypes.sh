#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/haplotypes
###Insert a new table with haplotype info to the strawberry database,
#Insert a new entry into the pipeline_run table describing the haplotyping pipeline.
python $scripts/insert_into_pipeline.py
#Oops, I realised the entry for this pipeline is already in the db!
#Remove it now.
DELETE FROM pipeline_run WHERE id=4;
#Create a new table schema.
source /home/sobczm/bin/popgen/snp/haplotype.sql
##Prepare a table gathering info to be loaded into the database.
#Extract the marker info from the alias table.
echo "SELECT * FROM marker" | mysql -u marias -h mongo -D strawberry_copy -p$(cat /home/sobczm/bin/mysql_sample_database/login)> marker.tsv
#Extract the alias info from the alias table. 
echo "SELECT * FROM alias" | mysql -u marias -h mongo -D strawberry_copy -p$(cat /home/sobczm/bin/mysql_sample_database/login)> alias.tsv
#Extract the sample info from the sample table.
echo "SELECT * FROM sample" | mysql -u marias -h mongo -D strawberry_copy -p$(cat /home/sobczm/bin/mysql_sample_database/login)> sample.tsv
#Extract the genotype info from the genotype table. 
echo "SELECT * FROM genotype" | mysql -u marias -h mongo -D strawberry_copy -p$(cat /home/sobczm/bin/mysql_sample_database/login)> genotype.tsv
#Extract the markers to be analyzed - ie. the ones which map to just one subgenome.
cut -d"," -f1 vesca2consensus_map_noduplicates_2017-05-17.csv | sort | uniq -c | sort -rk1 >markers_subgenomes.txt
#The markers with more than 1 occurrence (193 of them) will be eliminated from
#inclusion in the database.
cut -d"," -f1 vesca2consensus_map_noduplicates_2017-05-17.csv | sort | uniq -c | sort -rk1 | awk '$1 == 2' | sed 's/^[ \t]*//' | cut -d " " -f2 >markers_to_eliminate.txt
#Print table with marker_name, genotype_id, phased_genotype
python $scripts/prepare_haplotypes_table.py marker.tsv sample.tsv genotype.tsv alias.tsv markers_to_eliminate.txt
#Need to add key indices, doing it hack-ish way here:
cat -n haplotype_table.txt >haplotype_table_keys.txt
#Load the data into the database.
#Get errors using those methods:
source /home/sobczm/bin/popgen/snp/load_haplotype_table.sql
mysqlimport --fields-terminated-by=\t --local -u root -p strawberry_copy haplotype_table_keys.txt
#Trying Rob's way
python $scripts/load_haplotype_table.py haplotype_table.txt