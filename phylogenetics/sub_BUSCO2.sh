#$ -S /bin/bash
#$ -cwd
#$ -pe smp 5
#$ -l h_vmem=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace


### BUSCO analysis to identify single copy genes conserved in Fungi in all genomes in the study. Sample submission script for one genome.
### Note: use the genome mode for the transcriptome contigs, otherwise get all contigs re-named as Transcript 1 in the output!

### Do not forget to input the path to the assembly
### And select the database (options: Eukaryotic, Fungal, Plant, Bacteria)
### NOTE that if you are working with Fungi or Bacteria, you can download lower taxonomic level
### databases to use (e.g. just for Sordariomyceta)

# NOTE2 to easily prepare a figure with a comparison of complete, fragmented, suplicated
# and missing BUSCO genes in each genome, use the script BUSCO_plot.py in the BUSCO folder.
# Instructions in the user guide BUSCO_v2.0_userguide.pdf
###

assembly=$1
db_opt=$2

if [ "$db_opt" = "Eukaryotic" ]
then
    db=/home/sobczm/bin/BUSCO_v2/eukaryota_odb9
elif [ "$db_opt" = "Fungal" ]
then
    db=/home/sobczm/bin/BUSCO_v2/fungi_odb9
elif [ "$db_opt" = "Plant" ]
then
    db=/home/sobczm/bin/BUSCO_v2/embryophyta_odb9
elif [ "$db_opt" = "Bacteria" ]
then
    db=/home/sobczm/bin/BUSCO_v2/bacteria_odb9
else
    exit
fi


### Output folder
filename=$(basename "$assembly")
name="${filename%.*}_$(date +%F_%R)" #current assembly filename + timestamp for output dir

### Setting variables
cpath=$PWD
temp_dir="$TMPDIR"

busco=/home/sobczm/bin/BUSCO_v2/BUSCO.py

### Prep
mkdir -p $temp_dir
cp $assembly $temp_dir
cd $temp_dir

### Execute
python $busco -o $name -i $filename -l $db -m geno -c 5

### Cleanup
cp $temp_dir/$filename $temp_dir/run_$name
cp -r $temp_dir/run_$name $cpath
rm -r $temp_dir
