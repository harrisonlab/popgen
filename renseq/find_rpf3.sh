#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq/rpf2
input=/home/sobczm/popgen/renseq/strawberry

#Discover the sequence of the 3 markers for Rpf2 in vesca 1.1 
mysql -u strawreadonly -h mongo -D strawberry_samples -p$(cat /home/sobczm/.my.cnf)
mysql -u strawreadonly -h mongo -D strawberry_samples -p$(cat /home/sobczm/bin/mysql_sample_database/login)
SELECT probe_id, snp_id, alleles, strand, context FROM marker WHERE probe_id = "AX-89882775" OR probe_id = "AX-89825497" OR probe_id = "AX-89785467";

---------------------------------+
| probe_id    | snp_id        | alleles | strand | context                                                       |
+-------------+---------------+---------+--------+---------------------------------------------------------------+
| AX-89785467 | Affx-88833022 | T/C     | f      | GAAATATGTGTCTAAATTTCCATTCATTAT/ATGCAGATCCCTTGCAATTAACTCTTGCGC |
| AX-89825497 | Affx-88833139 | A/G     | f      | GAGGATTGGTTCAGGGACATTGCTTGGGAC/ACTCACACTTATGTTGTTGGTAGTGGCACA |
| AX-89882775 | Affx-88833184 | T/G     | r      | CCTTTTGCACATCTGCTCTGTTATAGAATT/TTTTGAATAATCATCCAGACATGGATCATA |
+-------------+---------------+---------+--------+---------------------------------------------------------------+


 select * from sample join clone on sample.clone_id = clone.id where name like '%red%';
+------+----------+---------------------------------------+------------------------------------------------------------------------------------------------------+----------+------------+------------------------+------+-------------------------------+--------+--------+
| id   | clone_id | file                                  | path                                                                                                 | type     | date       | batch                  | id   | name                          | mat_id | pat_id |
+------+----------+---------------------------------------+------------------------------------------------------------------------------------------------------+----------+------------+------------------------+------+-------------------------------+--------+--------+
| 1024 |        1 | EMR1_6_E06_Redgauntlet.CEL            | /home/groups/harrisonlab/raw_data/raw_celfiles/istraw90_6plates_RGxHA_EMxFE_FLxCH_BSxEL_etc/symlinks | istraw90 | 2014-12-05 | istraw90_batch1_plate6 |    1 | Redgauntlet                   |   NULL |   NULL |
| 1113 |        1 | EMR2_P2_C03_Redgauntlet.CEL           | /home/groups/harrisonlab/raw_data/raw_celfiles/istraw90_2plates_P150761Q151138_20160722/symlinks     | istraw90 | 2016-07-22 | istraw90_batch2_plate2 |    1 | Redgauntlet                   |   NULL |   NULL |
| 1294 |     3228 | ROS_175_C06_false_Redgauntlet.CEL     | /home/groups/harrisonlab/raw_data/raw_celfiles/istraw90_4plates_rosbreed/symlinks/known              | istraw90 | 2013-03-19 | rosbreed_plate4        | 3228 | ROS_175_C06_false_Redgauntlet |   NULL |   NULL |
| 1377 |     3680 | ROS_173_B09_Redglow.CEL               | /home/groups/harrisonlab/raw_data/raw_celfiles/istraw90_4plates_rosbreed/symlinks/known              | istraw90 | 2013-03-19 | rosbreed_plate2        | 3680 | ROS_173_B09_Redglow           |   NULL |   NULL |
| 1513 |        1 | B09-P160793-B1_Q3-A05_Redgauntlet.CEL | /home/groups/harrisonlab/raw_data/raw_celfiles/istraw35_1plate_P160793_B1                            | istraw35 | 2017-02-16 | istraw35_batch1_plate1 |    1 | Redgauntlet                   |   NULL |   NULL |
+------+----------+---------------------------------------+------------------------------------------------------------------------------------------------------+----------+------------+------------------------+------+-------------------------------+--------+--------+
5 rows in set (0.01 sec)

#Check alias for SNPs
SELECT * FROM alias WHERE probe_id = "AX-89882775" OR probe_id = "AX-89825497" OR probe_id = "AX-89785467";
mysql> SELECT * FROM alias WHERE probe_id = "AX-89882775" OR probe_id = "AX-89825497" OR probe_id = "AX-89785467";
+--------+----------+-------------+---------------+-----------+
| id     | platform | probe_id    | snp_id        | marker_id |
+--------+----------+-------------+---------------+-----------+
|   7171 | istraw90 | AX-89785467 | Affx-88833022 |      7171 |
| 173094 | istraw35 | AX-89785467 | Affx-88833022 |      7171 |
|  47164 | istraw90 | AX-89825497 | Affx-88833139 |     47164 |
| 173101 | istraw35 | AX-89825497 | Affx-88833139 |     47164 |
|  98092 | istraw90 | AX-89882775 | Affx-88833184 |     98092 |
+--------+----------+-------------+---------------+-----------+

#What alleles are present in RG?
##istraw90 and istraw35
#First marker
select * FROM genotype WHERE (sample_id = 1024 OR sample_id = 1113 OR sample_id = 1513) AND (alias_id = 7171 OR alias_id = 173094);
+---------+-------------+-----------+----------+------------+----------+-----------+
| id      | pipeline_id | sample_id | alias_id | clustering | genotype | haplotype |
+---------+-------------+-----------+----------+------------+----------+-----------+
| 7527352 |           1 |      1513 |   173094 | NMH        | BB       | NULL      |
| 9961100 |           2 |      1024 |     7171 | NMH        | AB       | NULL      |
| 9961127 |           2 |      1113 |     7171 | NMH        | AB       | NULL      |
+---------+-------------+-----------+----------+------------+----------+-----------+
3 rows in set (40.34 sec)

#Second marker
select * FROM genotype WHERE (sample_id = 1024 OR sample_id = 1113 OR sample_id = 1513) AND (alias_id = 47164 OR alias_id = 173101);
+----------+-------------+-----------+----------+------------+----------+-----------+
| id       | pipeline_id | sample_id | alias_id | clustering | genotype | haplotype |
+----------+-------------+-----------+----------+------------+----------+-----------+
| 11608100 |           2 |      1024 |    47164 | NMH        | AB       | NULL      |
| 11608127 |           2 |      1113 |    47164 | NMH        | AB       | NULL      |
+----------+-------------+-----------+----------+------------+----------+-----------+
2 rows in set (40.32 sec)

#Third marker
select * FROM genotype WHERE (sample_id = 1024 OR sample_id = 1113 OR sample_id = 1513) AND (alias_id = 98092);
----+
| id       | pipeline_id | sample_id | alias_id | clustering | genotype | haplotype |
+----------+-------------+-----------+----------+------------+----------+-----------+
| 13368300 |           2 |      1024 |    98092 | NMH        | AB       | NULL      |
| 13368327 |           2 |      1113 |    98092 | NMH        | AB       | NULL      |
+----------+-------------+-----------+----------+------------+----------+-----------+
2 rows in set (34.56 sec)

#What haplotype is it?

