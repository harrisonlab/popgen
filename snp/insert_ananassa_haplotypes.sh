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
mysql> DESCRIBE haplotype;
+-----------------+-------------+------+-----+---------+----------------+
| Field           | Type        | Null | Key | Default | Extra          |
+-----------------+-------------+------+-----+---------+----------------+
| haplo_id        | bigint(20)  | NO   | PRI | NULL    | auto_increment |
| genotype_id     | bigint(20)  | NO   |     | NULL    |                |
| pipeline_id     | int(11)     | NO   |     | NULL    |                |
| chromosome      | varchar(45) | YES  |     | NULL    |                |
| subgenome       | varchar(45) | YES  |     | NULL    |                |
| phased_genotype | varchar(45) | NO   |     | NULL    |                |
+-----------------+-------------+------+-----+---------+----------------+
6 rows in set (0.01 sec)

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
#Worked!

#Extract phased genotypes from the db and convert to VCF.
python $scripts/ananassa_haplotypes_db.py samples_to_analyze.txt samples_to_analyze.out

#Output the genotypes in the VCF format with locations substituted 
#according to map positions relative to vesca 1.1. 
python $scripts/ananassa_haplotypes_vcf.py samples_to_analyze.out istraw90_vesca_v1.1_snp_positions.gff3

python $scripts/substitute_sample_names.py samples_to_analyze.out.vcf cultivar_names.txt sample_clone_id.txt

#Utilise old scripts involving haplotype-based stats.