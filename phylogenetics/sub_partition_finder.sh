#$ -S /bin/bash
#$ -cwd 
#$ -pe smp 16
#$ -l h_vmem=1G 
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace|blacklace12.blacklace


#### Submit a Partition Finder job in a folder specified 

input=$1

python=/home/sobczm/bin/anaconda2/bin/python 
partition_finder=/home/sobczm/bin/PartitionFinder1.1.1/PartitionFinder.py
 
 $python $partition_finder $input