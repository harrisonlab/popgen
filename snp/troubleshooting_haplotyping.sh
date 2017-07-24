#!/bin/bash

#
# dump genotype data for haplotyping all samples
# datadir: /home/sobczm/popgen/snp/snp_chip/haplotypes/comparison/haplotyping4_mapping_pops
#

set -eu

input=/home/vicker/octoploid_mapping/haplotyping4
datadir=/home/sobczm/popgen/snp/snp_chip/haplotypes/comparison/haplotyping4_mapping_pops
scripts=/home/sobczm/bin/popgen/snp
export PATH=${scripts}:${PATH}
export PATH=/home/vicker/programs/shapeit.v2.r837.GLIBCv2.12.Linux.static/bin:${PATH}
export PATH=/home/vicker/programs/plink-1.90beta:${PATH}
source ${scripts}/haplotyping_all_samples_funcs.sh
export NODES="-N2,5,6,7,8,9,10"
export MAXJOBS="-L14"

cd $datadir/genotypes-flch_shapeit
shapeit -convert --input-haps 1A.phased --output-vcf 1A.phased.haps.vcf
#Subtitute B for C
cat 1A.phased.haps.vcf | awk -v OFS='\t' '{ if ($4 == "B") {$4 = "C"}; print }' | awk -v OFS='\t' '{ if ($5 == "B") {$5 = "C"}; print }' >1A.phased.haps.noB.vcf
#Drop lines where alt==ref
cat 1A.phased.haps.noB.vcf | awk -v OFS='\t' '$4!=$5 {print }' >1A.phased.haps.noB_nohoms.vcf
#Phase with Beagle
qsub /home/sobczm/bin/popgen/snp/sub_beagle.sh 1A.phased.haps.noB_nohoms.vcf
#No improvement observed regarding flipping.

#Compare input and output of shapeit
#Need to substitute B with Cs in the unphased file
sed 's/B/C/g' 1A.ped >1A_C.ped
python $scripts/overall_haplotype_similarity.py 1A.ped.phased 1A_C.ped 1A.info 1A.info 
#All randmom matches as before, so shapeit actually does something to the input files

#Carry out iterative haplotyping - haplotype the output of shapeit.
#Use the VCF conversion from above as input.
grid_run -Jshapeit ${MAXJOBS} ${NODES} -M4 -C4 \
shapeit --input-vcf 1A.phased.haps.noB.vcf \
--output-max 1A.phased.haps.noB.haps 1A.phased.haps.noB.sample \
--thread 4 
#And extract to ped
extract_haplotypes_to_ped.py 1A_C.ped 1A.phased.haps.noB.haps 1A.phased.haps.noB.ped
sed 's/B/C/g' 1A.ped.phased >1A_C.ped.phased
#Compare with map-based haplotypes.
python $scripts/overall_haplotype_similarity.py 1A.phased.haps.noB.ped ./comparison/1A_phased_C.ped 1A.info ./comparison/1A.info 

#Shapeit3 
grid_run -Jshapeit3.r884.1 ${MAXJOBS} ${NODES} -M4 -C4 shapeit3.r884.1 -G 1A.phased -M 1A.pmap -O 1A.shapeit3 --cluster-size 4000