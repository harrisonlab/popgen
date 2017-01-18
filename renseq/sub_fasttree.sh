
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 5
#$ -l virtual_free=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace

input=$1
output=${input%.fa*}.nwk
FastTreeMP < $input > $output
