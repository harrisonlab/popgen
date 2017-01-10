#$ -S /bin/bash
#$ -cwd
#$ -pe smp 2
#$ -l h_vmem=7G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

input=$1
output="{$input%.vcf}_haplo.vcf.gz"
beagle=/home/sobczm/bin/beagle/beagle.27Jul16.86a.jar

java -Xmx12g -jar $beagle gt=$input out=$output nthreads=20 window=1000 overlap=100` gprobs=false
