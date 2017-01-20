#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l h_vmem=3G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

#Script will phase diploid genotypes in a given VCF file. 
#Output: VCF file with phased genotypes with the suffix "_haplo"

input=$1
output=${input%.vcf}_haplo
beagle=/home/sobczm/bin/beagle/beagle.27Jul16.86a.jar

java -Xmx8g -Xms1g -jar $beagle gt=$input out=$output nthreads=20 window=1000 overlap=100 gprobs=false
