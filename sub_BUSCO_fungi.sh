#!/bin/bash
#$ -S /bin/bash
#$ -cwd 
#$ -pe smp 1
#$ -l h_vmem=1G 

### BUSCO analysis to identify single copy genes conserved in Fungi in all genomes in the study. Sample submission script for one genome:
### Note: use the genome mode for the transcriptome contigs, otherwise get all contigs re-named as Transcript 1 in the output!

busco=/home/sobczm/bin/BUSCO_v1.22/BUSCO_v1.22.py
db=/home/sobczm/bin/BUSCO_v1.22/fungi
name=${PWD##*/} #current directory name

## To change:
assembly=PATH/TO/ASSEMBLY/INVESTIGATED

python $busco -o $name -in $assembly -l $db -m genome



