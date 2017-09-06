#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/final_samples
cd $input

#Dump the clone table
a="SELECT * FROM clone"
echo $a | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login)> all_clones.txt

#Dump the sample table
a="SELECT id, clone_id, file, path, type, date, batch FROM sample"
echo $a | mysql -u marias -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login)> all_samples.txt
