#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#RepeatModeler and Masker submission script based on Andy's.

# This script uses repeatmodeler and repeatmasker to mask Interspersed repeats
# and low complexity regions within the genome. Firstly, repeatmodeler identifies
# repeat element boundaries and relationships within repeat families. The repeats
# identified within the genome are provided to repeatmasker, which uses this data
# along with it's own repeat libraries to identify these repetitive regions and
# perform masking. Masking is done at 3 levels:
# Hardmasking = repetitive sequence is replaced with N's.
# Softmasking = repetitive sequence is converted to lower case.
# Ignoring low-complexity regions = only interspersed repetitive elements are masked.

assembly=$1
filename=$(basename "$assembly")
name=${filename%.*}
echo $name

repeat_masker=/home/armita/prog/RepeatMasker
repeat_modeler=/home/armita/prog/RepeatModeler
$repeat_modeler/BuildDatabase -name $name $assembly
$repeat_modeler/RepeatModeler -database $name

#Hardmask
$repeat_masker/RepeatMasker -gff -pa 4 -lib ./lsRM_*.*/consensi.fa.classified $assembly
#mv $filename.cat.gz "$name"_hardmasked.fasta.cat.gz
#mv $filename.masked "$name"_hardmasked.fasta
#mv $filename.out "$name"_hardmasked.out
#mv $filename.tbl "$name"_hardmasked.tbl
#grep -v '#' $filename.out.gff > "$name"_hardmasked.gff

#softmask
# repeat_masker/RepeatMasker -xsmall -gff -pa 4 -lib RM_*.*/consensi.fa.classified $assembly
#mv $filename.cat.gz "$name"_softmasked.fasta.cat.gz
#mv $filename.masked "$name"_softmasked.fasta
#mv $filename.out "$name"_softmasked.out
#mv $filename.tbl "$name"_softmasked.tbl
#grep -v '#' $filename.out.gff > "$name"_softmasked.gff
