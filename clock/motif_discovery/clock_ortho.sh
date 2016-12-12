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
#points when executing through SGE. Therefore, will run it on a worker node in a screen.

##Running a dummy script to reserve resources on blacklace01 (16 threads)
qsub /home/sobczm/bin/popgen/other/reserve_cluster.sh
#Running orthofinder on the directory OrthoFinder containing all input files
sh ./run_orthofinder.sh $input/$dir

##DATASET B)
#Get whole-genome protein sequences and rename the files
cd $input/pep_genomes
mv Chaetomium_globosum_cbs_148_51.ASM14336v1.pep.all.fa Chaetomium_globosum.pep.fa
mv Chaetomium_thermophilum_var_thermophilum_dsm_1495.CTHT_3.0.pep.all.fa Chaetomium_thermophilum.pep.fa
mv Fusarium_fujikuroi.EF1.pep.all.fa Fusarium_fujikuroi.pep.fa
mv Fusarium_langsethiae.ASM129263v1.pep.all.fa Fusarium_langsethiae.pep.fa
mv Fusarium_pseudograminearum.GCA_000303195.1.pep.all.fa Fusarium_pseudograminearum.pep.fa
mv Fusarium_verticillioides.ASM14955v1.pep.all.fa Fusarium_verticillioides.pep.fa
mv Gaeumannomyces_graminis.Gae_graminis_V2.pep.all.fa Gaeumannomyces_graminis.pep.fa
mv Magnaporthe_poae.Mag_poae_ATCC_64411_V1.pep.all.fa Magnaporthe_poae.pep.fa
mv Neurospora_tetrasperma_fgsc_2508.v2.0.pep.all.fa Neurospora_tetrasperma.pep.fa
mv Sclerotinia_borealis_f_4157.SBOR_1.pep.all.fa Sclerotinia_borealis.pep.fa
mv Sclerotinia_sclerotiorum.ASM14694v1.pep.all.fa Sclerotinia_sclerotiorum.pep.fa
mv Thielavia_terrestris_nrrl_8126.ASM22611v1.pep.all.fa Thielavia_terrestris.pep.fa
mv Trichoderma_atroviride_imi_206040.TRIAT_v2_0.pep.all.fa Trichoderma_atroviride.pep.fa
mv Trichoderma_gamsii.ASM148177v1.pep.all.fa Trichoderma_gamsii.pep.fa
mv Trichoderma_harzianum.ASM98886v1.pep.all.fa Trichoderma_harzianum.pep.fa
mv Trichoderma_virens.ASM17099v1.pep.all.fa Trichoderma_virens.pep.fa
mv Verticillium_longisporum_gca_001268165.vl2.denovo.v1.pep.all.fa Verticillium_longisporum.pep.fa

#Create a directory for OrthoFinder run, copy input FASTA files and run orthofinder
mkdir OrthoFinder2
for a in *fa; do cp $a ./OrthoFinder2; done
# Running a dummy script to reserve resources on blacklace01
qsub /home/sobczm/bin/popgen/other/reserve_cluster.sh
#Running orthofinder on the directory OrthoFinder containing all input files
dir=OrthoFinder2
sh ./run_orthofinder.sh $input/$dir
