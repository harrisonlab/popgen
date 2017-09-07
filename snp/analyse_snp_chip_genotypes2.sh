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

python $scripts/match_up_crown_rot_ids.py crown_rot_data.txt all_clones.txt all_samples.txt >congruent_ids_crown_rot

#Add sample id to the master strawberry table. Only perfect case-blind matches accepted. 
python $scripts/add_sample_id.py master_strawberry_spreadsheet.txt all_samples.txt  >master_strawberry_spreadsheet_sample_id.txt
python $scripts/add_sample_id.py just_sample_names_spreadsheet.txt all_samples.txt  >just_sample_names_spreadsheet_sample_id.txt