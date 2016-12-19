#!/bin/bash
input=/home/sobczm/popgen/clock/pep_genomes
scripts=/home/sobczm/bin/popgen/clock/motif_discovery

##Determine the orthogroup where the Neurospora protein ID of interest ( taken from
## Orthology_analysis.xlsx) is placed.

##Frequency
head -1 $input/OrthoFinder2/Results_Oct26/Orthogroups.csv >>freq_orthogroup
grep "ESA42013" $input/OrthoFinder2/Results_Oct26/Orthogroups.csv >>freq_orthogroup
#WC1
head -1 $input/OrthoFinder2/Results_Oct26/Orthogroups.csv >>wc1_orthogroup
grep "ESA41977" $input/OrthoFinder2/Results_Oct26/Orthogroups.csv >>wc1_orthogroup
#Frequency
head -1 $input/OrthoFinder2/Results_Oct26/Orthogroups.csv >>wc2_orthogroup
grep "EAA34583" $input/OrthoFinder2/Results_Oct26/Orthogroups.csv >>wc2_orthogroup
#frh
head -1 $input/OrthoFinder2/Results_Oct26/Orthogroups.csv >>frh_orthogroup
grep "EAA27062" $input/OrthoFinder2/Results_Oct26/Orthogroups.csv >>frh_orthogroup
#Fwd-1
head -1 $input/OrthoFinder2/Results_Oct26/Orthogroups.csv >>fwd1_orthogroup
grep "EAA26744" $input/OrthoFinder2/Results_Oct26/Orthogroups.csv >>fwd1_orthogroup
#vvd
head -1 $input/OrthoFinder2/Results_Oct26/Orthogroups.csv >>vvd_orthogroup
grep "EAA28370" $input/OrthoFinder2/Results_Oct26/Orthogroups.csv >>vvd_orthogroup

##Retrieve the protein sequences in question and save them to individual files per gene
##per species
##(File can contain more than 1 protein sequences if multiple orthologs present)
python $scripts/extract_protein.py freq_orthogroup freq
python $scripts/extract_protein.py wc1_orthogroup wc1
python $scripts/extract_protein.py wc2_orthogroup wc2
python $scripts/extract_protein.py frh_orthogroup frh
python $scripts/extract_protein.py fwd1_orthogroup fwd1
python $scripts/extract_protein.py vvd_orthogroup vvd

#Run InterProScan on each gene in the folder:
for file in freq*.fa wc1*.fa wc2*.fa frh*.fa fwd1*.fa vvd*.fa
do
    Jobs=$(qstat | grep 'sub_interp' | wc -l)
    while [ $Jobs -gt 100 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_interp' | wc -l)
    done
qsub /home/sobczm/bin/popgen/renseq/sub_interproscan.sh $file
done
