#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace


### BUSCO analysis to identify single copy genes conserved in Fungi in all genomes in the study. Sample submission script for one genome.
### Note: use the genome mode for the transcriptome contigs, otherwise get all contigs re-named as Transcript 1 in the output!

### Do not forget to input the path to the assembly
assembly=$1

### Output folder
filename=$(basename "$assembly")
name="${filename%.*}_$(date +%F_%R)" #current assembly filename + timestamp for output dir

### Setting variables
cpath=$PWD
temp_dir="$TMPDIR"

busco=/home/sobczm/bin/BUSCO_v1.22/BUSCO_v1.22.py
db=/home/sobczm/bin/BUSCO_v1.22/fungi

### Prep
mkdir -p $temp_dir
cd $temp_dir
cp $assembly $temp_dir

### Execute
python $busco -o $name -in $assembly -l $db -m genome

### Cleanup
mv $temp_dir/$filename $temp_dir/run_$name
mv $temp_dir/run_$name $cpath
rm -r $temp_dir
