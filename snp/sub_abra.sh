#$ -S /bin/bash
#$ -cwd
#$ -pe smp 8
#$ -l h_vmem=2G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace|blacklace12.blacklace

### ABRA for realignment of reads around indels

input_bam=$1
reference=$2
filename=$(basename "$input_bam")
reference_b=$(basename "$input_bam")
name="${filename%.*}"
output_dir="${name%.*}_$(date +%F_%R)"

### Output folder
cpath=$PWD
temp_dir="$TMPDIR"

### Temp folder prep
mkdir -p $temp_dir
cp $input_bam $temp_dir
cp $reference $temp_dir
cd $temp_dir

### Create BED file
btb=/home/armita/prog/bedtools2/bin/bamToBed
$btb -i $filename >$name.bed

abra=/home/sobczm/bin/abra/bin/abra-0.97-SNAPSHOT-jar-with-dependencies.jar
java -Xmx4G -jar $abra --in $filename --out $name\_abra.bam --ref $reference_b \
--targets $name.bed --threads 8 --working $output_dir > $name.log 2>&1

samtools sort $name\_abra.bam  $name\_abra_sorted
samtools index $name\_abra_s.bam

picard=/home/sobczm/bin/picard-tools-2.5.0/picard.jar
java -jar $picard FixMateInformation \
I=$name\_abra_s.bam O=$name\_abra_sf.bam

mv * $cpath
