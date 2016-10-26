#!/bin/bash
input=/home/sobczm/popgen/clock
scripts=/home/sobczm/bin/popgen/clock/motif_discovery

#Identification of orthologs amongst selected genomes using the OrthoFinder.

##DATASET A)
#Get whole-genome protein sequences and rename the files
cd $input/pep_genomes
mv Botrytis_cinerea.ASM83294v1.pep.all.fa Botrytis_cinerea.pep.fa
mv Fus2_final_genes_combined.pep.fasta Fus2.pep.fa
mv Fusarium_graminearum.RR.pep.all.fa Fusarium_graminearum.pep.fa
mv Fusarium_solani.v2.0.pep.all.fa Fusarium_solani.pep.fa
mv Magnaporthe_oryzae.MG8.pep.all.fa Magnaporthe_oryzae.pep.fa
mv Neonectria_ditissima.R0905_v2.0.pep.all.fa Neonectria_ditissima.pep.fa
mv Neurospora_crassa.NC12.pep.all.fa Neurospora_crassa.pep.fa
mv Podospora_anserina_s_mat_.ASM22654v1.pep.all.fa Podospora_anserina.pep.fa
mv Sordaria_macrospora.ASM18280v2.pep.all.fa Sordaria_macrospora.pep.fa
mv Trichoderma_reesei.GCA_000167675.2.pep.all.fa Trichoderma_reesei.pep.fa
mv Verticillium_alfalfae_vamsstat_102.ASM15082v1.pep.all.fa Verticillium_alfalfae.pep.fa
mv Verticillium_dahliae.ASM15067v2.pep.all.fa Verticillium_dahliae.pep.fa
#Create a directory for OrthoFinder run, copy input FASTA files and run orthofinder
mkdir OrthoFinder
for a in *fa; do cp $a ./OrthoFinder; done
qsub $scripts/sub_orthofinder.sh OrthoFinder

#It turns out that the pipeline terminates for no reason at different random
#points when executing through SGE. Therefore, will run it on a worker node.

# Running a dummy script to reserve resources on blacklace01
qsub /home/sobczm/bin/popgen/other/reserve_cluster.sh
#Running orthofinder on the directory OrthoFinder containing all input files
orthofinder=/home/sobczm/bin/OrthoFinder-1.0.7/orthofinder
anaconda=/home/sobczm/bin/anaconda2/bin/python
dir=OrthoFinder
$anaconda $orthofinder/orthofinder.py -f $dir -t 16

##DATASET B)
