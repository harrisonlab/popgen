#$ -S /bin/bash
#$ -cwd
#$ -pe smp 8
#$ -l h_vmem=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace|blacklace12.blacklace
###

bam=$1
splitters=$2
discordant=$3

cpath=$PWD

temp_dir="$TMPDIR"
mkdir -p $temp_dir

cp -r $bam $splitters $discordant $temp_dir
bam_in=$(basename "$bam")
splitters_in=$(basename "$splitters")
discordant_in=$(basename "$discordant")

cd $temp_dir

lumpy=/home/sobczm/bin/lumpy-sv/bin
output=${bam_in%.*}.gt

$lumpy/lumpyexpress -B $bam_in -S $splitters_in -D $discordant_in -o $output 

samtools index $bam_in
samtools index $splitters_in
samtools index $discordant_in

svtyper=/home/sobczm/bin/svtyper
$svtyper/svtyper -B $bam_in -S $splitters_in -i $output >${bam_in%.*}.vcf

rm $bam_in $splitters_in $discordant_in
cp -r * $cpath