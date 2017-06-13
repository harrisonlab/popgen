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
#Extract the sample info from the sample table.
echo "SELECT * FROM sample" | mysql -u marias -h mongo -D strawberry_copy -p$(cat /home/sobczm/bin/mysql_sample_database/login)> sample.tsv
#Extract the genotype info from the genotype table. 
echo "SELECT * FROM genotype" | mysql -u marias -h mongo -D strawberry_copy -p$(cat /home/sobczm/bin/mysql_sample_database/login)> gentype.tsv
#Load the data into the database.

